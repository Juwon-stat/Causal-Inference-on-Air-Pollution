###############################################################################
# Title   : CBPS_GLOBAL_FINAL.R
# Purpose : Block-wise Covariate Balancing Propensity Score (CBPS)
#           + block-specific ATT
#           + Hajek aggregation
#           + conditional SE
#           + block bootstrap
#           + unit-level multiplier bootstrap
#           + SMD diagnostics for before-vs-after balance plots
#
# ATT direction:
#   ATT = weighted control mean - treated mean
#   Positive ATT means PM concentration was lower during the lockdown period.
#
# Design:
#   - Observational unit: district-date-hour record
#   - Treatment-date block: one block per treated calendar date
#   - Control pool: non-2020 observations within ±d calendar days
#   - CBPS model is fitted within each block only
#   - CBPS weights are treated as fixed for inference
#
# Output:
#   - cbps_all_att_summary.csv
#   - cbps_all_delta_table.csv
#   - cbps_all_block_att_summary.csv
#   - cbps_all_hajek_att_summary.csv
#   - cbps_all_obs_contrib.csv
#   - cbps_all_smd_variable_summary.csv
#   - cbps_all_smd_plot_summary.csv
#   - CBPS_SMD_values_master_rt_weighted_mean.csv
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(CBPS)
})

# =============================================================================
# 0. USER CONFIGURATION
# =============================================================================

cfg <- list(
  root_in  = "C:/path/to/Seoul_RData",
  root_out = "C:/path/to/CBPS_results",
  
  q_values = c(3, 4, 5),
  day_windows = c(3, 7, 14),
  region_vars = c("dist", "area"),
  
  treatment_start = as.Date("2020-01-24"),
  treatment_end   = as.Date("2020-02-09"),
  treatment_year  = 2020L,
  
  outcomes = c("PM10", "PM25"),
  
  # CBPS settings
  CBPS_METHOD = "exact",
  ATT_TRIM_EPS = 1e-6,
  
  # Stability settings:
  # TRUE means sparse / separated covariate levels are screened block-by-block.
  SCREEN_COVARIATES = TRUE,
  
  # Bootstrap
  B = 1000L,
  seed = 20260507L,
  
  # Save options
  save_full_rds = FALSE,
  save_weighted_rows_csv = FALSE
)

dir.create(cfg$root_out, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. BASIC HELPERS
# =============================================================================

safe_as_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  as.Date(x)
}

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

make_month_day_window <- function(date_i, d) {
  format(seq(date_i - d, date_i + d, by = "day"), "%m-%d")
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  mu <- weighted_mean(x[ok], w[ok])
  sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok])
}

get_all_covariates <- function(region_var) {
  c(
    region_var,
    "wday",
    "time_of_day",
    "CO_quantile",
    "O3_quantile",
    "wind_sp_ind",
    "wind_dir_8",
    "rain_binary",
    "temp_round"
  )
}

fit_cbps_capture <- function(formula, data, cfg) {
  warning_msgs <- character(0)
  
  fit <- withCallingHandlers(
    tryCatch(
      CBPS::CBPS(
        formula = formula,
        data = data,
        ATT = 1,
        method = cfg$CBPS_METHOD
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warning_msgs <<- c(warning_msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  list(
    fit = fit,
    warnings = unique(warning_msgs),
    ok = !inherits(fit, "error")
  )
}

extract_cbps_ps <- function(fit, newdata) {
  ps <- NULL
  
  ps <- tryCatch(
    as.numeric(predict(fit, newdata = newdata, type = "response")),
    error = function(e) NULL
  )
  
  if (!is.null(ps) && length(ps) == nrow(newdata)) {
    return(ps)
  }
  
  if (!is.null(fit$fitted.values)) {
    ps <- as.numeric(fit$fitted.values)
    if (length(ps) == nrow(newdata)) return(ps)
  }
  
  if (!is.null(fit$ps)) {
    ps <- as.numeric(fit$ps)
    if (length(ps) == nrow(newdata)) return(ps)
  }
  
  if (!is.null(fit$propensity.score)) {
    ps <- as.numeric(fit$propensity.score)
    if (length(ps) == nrow(newdata)) return(ps)
  }
  
  stop("Could not extract propensity scores from CBPS object.", call. = FALSE)
}

screen_cbps_covariates <- function(pool, candidate_covariates) {
  keep <- character(0)
  dropped <- list()
  
  for (v in candidate_covariates) {
    x <- pool[[v]]
    z <- pool$treated
    
    if (length(unique(x[!is.na(x)])) < 2) {
      dropped[[length(dropped) + 1L]] <- data.table(
        covariate = v,
        reason = "single_level_in_block"
      )
      next
    }
    
    tab <- table(z, x, useNA = "no")
    
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      dropped[[length(dropped) + 1L]] <- data.table(
        covariate = v,
        reason = "insufficient_2way_table"
      )
      next
    }
    
    level_has_both_groups <- apply(tab, 2, function(col) all(col > 0))
    
    if (!all(level_has_both_groups)) {
      dropped[[length(dropped) + 1L]] <- data.table(
        covariate = v,
        reason = "level_not_shared_by_treated_and_control"
      )
      next
    }
    
    keep <- c(keep, v)
  }
  
  drop_dt <- if (length(dropped) > 0) {
    rbindlist(dropped, use.names = TRUE, fill = TRUE)
  } else {
    data.table()
  }
  
  list(keep = unique(keep), dropped = drop_dt)
}

# =============================================================================
# 2. DATA PREPARATION
# =============================================================================

load_q_data <- function(q, root_in) {
  file_path <- file.path(root_in, sprintf("seoul_q%d.RData", q))
  
  if (!file.exists(file_path)) {
    stop(sprintf("Input file not found: %s", file_path), call. = FALSE)
  }
  
  env <- new.env()
  load(file_path, envir = env)
  
  obj_name <- sprintf("seoul_q%d", q)
  
  if (!exists(obj_name, envir = env)) {
    stop(
      sprintf("Object `%s` not found inside %s", obj_name, file_path),
      call. = FALSE
    )
  }
  
  as.data.table(get(obj_name, envir = env))
}

prepare_data <- function(DT, q, cfg) {
  DT <- copy(DT)
  
  DT[, date := safe_as_date(date)]
  
  if (!("year" %in% names(DT))) {
    DT[, year := as.integer(format(date, "%Y"))]
  }
  
  setorder(DT, date)
  DT[, row_id := .I]
  DT[, month_day := format(date, "%m-%d")]
  
  # Temperature covariate
  DT[, temp_round := round(temp)]
  
  # Observation-level ID for reuse-aware multiplier bootstrap.
  # Same district-date-hour observation receives same multiplier when reused.
  DT[, obs_unit_id := paste(dist, as.character(date), hour, sep = "_")]
  
  for (outcome in cfg$outcomes) {
    DT[, (outcome) := to_num(get(outcome))]
  }
  
  DT
}

# =============================================================================
# 3. SMD HELPERS
# =============================================================================

smd_binary_weighted <- function(x_t, w_t, x_c, w_c) {
  mt <- weighted_mean(x_t, w_t)
  mc <- weighted_mean(x_c, w_c)
  
  vt <- weighted_var(x_t, w_t)
  vc <- weighted_var(x_c, w_c)
  
  sp <- sqrt((vt + vc) / 2)
  
  if (!is.finite(sp) || sp == 0) {
    if (is.finite(mt) && is.finite(mc) &&
        abs(mt - mc) < .Machine$double.eps) {
      return(0)
    }
    return(NA_real_)
  }
  
  abs(mt - mc) / sp
}

compute_smd_categorical <- function(treated_dt,
                                    control_dt,
                                    covariates,
                                    treated_weight_col = "smd_weight",
                                    control_weight_col = "smd_weight",
                                    label = "after") {
  out <- list()
  
  for (v in covariates) {
    xt <- as.character(treated_dt[[v]])
    xc <- as.character(control_dt[[v]])
    
    wt <- treated_dt[[treated_weight_col]]
    wc <- control_dt[[control_weight_col]]
    
    lvls <- sort(unique(c(xt, xc)))
    lvls <- lvls[!is.na(lvls)]
    
    if (length(lvls) == 0) next
    
    tmp <- lapply(lvls, function(lv) {
      zt <- as.numeric(xt == lv)
      zc <- as.numeric(xc == lv)
      
      data.table(
        sample = label,
        covariate = v,
        level = lv,
        smd = smd_binary_weighted(zt, wt, zc, wc),
        treated_mean = weighted_mean(zt, wt),
        control_mean = weighted_mean(zc, wc)
      )
    })
    
    out[[v]] <- rbindlist(tmp, use.names = TRUE, fill = TRUE)
  }
  
  if (length(out) == 0) return(data.table())
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

summarize_smd_by_variable <- function(smd_level_dt) {
  if (nrow(smd_level_dt) == 0) return(data.table())
  
  smd_level_dt[
    ,
    .(
      max_abs_smd = suppressWarnings(max(smd, na.rm = TRUE)),
      mean_abs_smd = suppressWarnings(mean(smd, na.rm = TRUE)),
      n_levels = .N
    ),
    by = .(sample, covariate)
  ][
    !is.infinite(max_abs_smd)
  ]
}

compute_block_smd <- function(pool_before,
                              weighted_after,
                              covariates,
                              block_date,
                              block_weight_before,
                              block_weight_after) {
  # BEFORE: eligible sample in the block before weighting
  before_treated <- pool_before[treated == 1L, c(covariates), with = FALSE]
  before_control <- pool_before[treated == 0L, c(covariates), with = FALSE]
  
  before_treated[, smd_weight := 1]
  before_control[, smd_weight := 1]
  
  smd_before_level <- compute_smd_categorical(
    treated_dt = before_treated,
    control_dt = before_control,
    covariates = covariates,
    label = "before"
  )
  
  smd_before_var <- summarize_smd_by_variable(smd_before_level)
  
  if (nrow(smd_before_var) > 0) {
    smd_before_var[
      ,
      `:=`(
        block_date = block_date,
        block_weight = block_weight_before
      )
    ]
  }
  
  # AFTER: CBPS weighted sample
  after_treated <- weighted_after[treated == 1L, c(covariates, "weights"), with = FALSE]
  after_control <- weighted_after[treated == 0L, c(covariates, "weights"), with = FALSE]
  
  setnames(after_treated, "weights", "smd_weight")
  setnames(after_control, "weights", "smd_weight")
  
  smd_after_level <- compute_smd_categorical(
    treated_dt = after_treated,
    control_dt = after_control,
    covariates = covariates,
    label = "after"
  )
  
  smd_after_var <- summarize_smd_by_variable(smd_after_level)
  
  if (nrow(smd_after_var) > 0) {
    smd_after_var[
      ,
      `:=`(
        block_date = block_date,
        block_weight = block_weight_after
      )
    ]
  }
  
  rbindlist(
    list(smd_before_var, smd_after_var),
    use.names = TRUE,
    fill = TRUE
  )
}

aggregate_smd_blocks_weighted_mean <- function(smd_block_variable) {
  if (nrow(smd_block_variable) == 0) return(data.table())
  
  smd_block_variable[
    is.finite(max_abs_smd) &
      is.finite(block_weight) &
      block_weight > 0,
    .(
      max_abs_smd = weighted.mean(max_abs_smd, w = block_weight, na.rm = TRUE),
      mean_abs_smd = weighted.mean(mean_abs_smd, w = block_weight, na.rm = TRUE),
      n_blocks = uniqueN(block_date)
    ),
    by = .(sample, covariate)
  ]
}

make_smd_master_for_panel <- function(smd_variable,
                                      q,
                                      region_var,
                                      unit = "rt",
                                      d_win) {
  dt <- copy(smd_variable)
  
  if (nrow(dt) == 0) return(data.table())
  
  dt[
    ,
    `:=`(
      variable = as.character(covariate),
      stage = tolower(as.character(sample)),
      region = region_var,
      unit = unit,
      q = paste0("q", q),
      d_win = as.integer(d_win),
      smd_value = as.numeric(max_abs_smd)
    )
  ]
  
  max_na <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    max(x)
  }
  
  wide <- dcast(
    dt,
    variable + region + unit + q + d_win ~ stage,
    value.var = "smd_value",
    fun.aggregate = max_na
  )
  
  if (!("before" %in% names(wide))) wide[, before := NA_real_]
  if (!("after" %in% names(wide))) wide[, after := NA_real_]
  
  wide[
    ,
    .(
      variable,
      before,
      after,
      region,
      unit,
      q,
      d_win
    )
  ]
}

# =============================================================================
# 4. ATT, CONDITIONAL SE, BLOCK BOOTSTRAP, MULTIPLIER BOOTSTRAP
# =============================================================================

make_delta_table <- function(rows, outcomes, weight_col = "weights") {
  # Delta_bu = weighted control mean in block b - treated outcome for u.
  out <- list()
  
  for (outcome_name in outcomes) {
    tmp <- rows[
      ,
      {
        tr <- .SD[treated == 1L]
        co <- .SD[treated == 0L]
        
        wt <- tr[[weight_col]]
        wc <- co[[weight_col]]
        
        yt <- tr[[outcome_name]]
        yc <- co[[outcome_name]]
        
        mu_c <- weighted_mean(yc, wc)
        mu_t <- weighted_mean(yt, wt)
        
        Wb <- sum(wt[is.finite(wt) & wt > 0], na.rm = TRUE)
        tau_b <- mu_c - mu_t
        
        data.table(
          outcome = outcome_name,
          obs_unit_id = tr$obs_unit_id,
          block_date = tr$block_date,
          omega_bu = wt,
          Y_treated = yt,
          control_mean_b = mu_c,
          Delta_bu = mu_c - yt,
          Wb = Wb,
          tau_b = tau_b,
          numerator_b = Wb * tau_b,
          n_treated = nrow(tr),
          n_control = nrow(co)
        )
      },
      by = block_date
    ]
    
    out[[outcome_name]] <- tmp
  }
  
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

compute_block_att_from_delta <- function(delta_dt) {
  unique(
    delta_dt[
      ,
      .(
        Wb = Wb[1],
        tau_b = tau_b[1],
        numerator = numerator_b[1],
        n_treated = n_treated[1],
        n_control = n_control[1]
      ),
      by = .(outcome, block_date)
    ]
  )
}

hajek_from_block_table <- function(block_table) {
  block_table[
    is.finite(Wb) & Wb > 0 & is.finite(numerator),
    .(
      att = sum(numerator, na.rm = TRUE) / sum(Wb, na.rm = TRUE),
      W_total = sum(Wb, na.rm = TRUE),
      n_blocks = uniqueN(block_date)
    ),
    by = outcome
  ]
}

conditional_se_from_delta <- function(delta_dt) {
  # Var_cond(tau_hat) =
  # sum_b a_b^2 * (1/W_b^2) * sum_u omega_bu^2 (Delta_bu - tau_b)^2.
  block_var <- delta_dt[
    is.finite(omega_bu) & omega_bu > 0 &
      is.finite(Delta_bu) &
      is.finite(Wb) & Wb > 0,
    .(
      tau_b = tau_b[1],
      Wb = Wb[1],
      var_delta_b = sum(omega_bu^2 * (Delta_bu - tau_b[1])^2, na.rm = TRUE) /
        (Wb[1]^2),
      n_treated = .N
    ),
    by = .(outcome, block_date)
  ]
  
  Wtot <- block_var[
    ,
    .(W_total = sum(Wb, na.rm = TRUE)),
    by = outcome
  ]
  
  block_var <- merge(block_var, Wtot, by = "outcome", all.x = TRUE)
  block_var[, a_b := Wb / W_total]
  
  out <- block_var[
    ,
    {
      var_total <- sum((a_b^2) * var_delta_b, na.rm = TRUE)
      se <- sqrt(var_total)
      att0 <- sum(Wb * tau_b, na.rm = TRUE) / sum(Wb, na.rm = TRUE)
      
      .(
        estimator = "conditional",
        att = att0,
        std_error = se,
        ci_low = att0 - 1.96 * se,
        ci_high = att0 + 1.96 * se,
        n = sum(n_treated),
        n_blocks = uniqueN(block_date)
      )
    },
    by = outcome
  ]
  
  out[]
}

block_bootstrap_from_block_table <- function(block_table, outcomes, B_boot, seed) {
  # Resample treatment-date blocks with replacement.
  set.seed(seed)
  
  out <- list()
  
  block_dates_all <- sort(unique(block_table$block_date))
  K <- length(block_dates_all)
  
  for (outcome_name in outcomes) {
    bt <- block_table[outcome == outcome_name]
    
    boot_stat <- numeric(B_boot)
    
    for (b in seq_len(B_boot)) {
      sampled_blocks <- sample(block_dates_all, size = K, replace = TRUE)
      
      kappa <- data.table(block_date = sampled_blocks)[, .N, by = block_date]
      setnames(kappa, "N", "kappa")
      
      boot_bt <- merge(bt, kappa, by = "block_date", all.x = FALSE)
      
      denom <- boot_bt[, sum(kappa * Wb, na.rm = TRUE)]
      
      boot_stat[b] <- if (is.finite(denom) && denom > 0) {
        boot_bt[, sum(kappa * numerator, na.rm = TRUE) / denom]
      } else {
        NA_real_
      }
    }
    
    att0 <- hajek_from_block_table(bt)$att[1]
    
    out[[outcome_name]] <- data.table(
      estimator = "block_bootstrap",
      outcome = outcome_name,
      att = att0,
      std_error = sd(boot_stat, na.rm = TRUE),
      ci_low = as.numeric(quantile(boot_stat, 0.025, na.rm = TRUE)),
      ci_high = as.numeric(quantile(boot_stat, 0.975, na.rm = TRUE)),
      n = B_boot,
      n_blocks = K
    )
  }
  
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

make_observation_contribution_table <- function(rows,
                                                block_table,
                                                outcomes,
                                                weight_col = "weights") {
  # Contribution representation for observation-level multiplier bootstrap:
  # tau_hat = sum_k A_k / sum_k B_k.
  out <- list()
  
  for (outcome_name in outcomes) {
    bt <- block_table[
      outcome == outcome_name,
      .(block_date, outcome, Wb)
    ]
    
    tmp <- merge(
      rows,
      bt,
      by = "block_date",
      allow.cartesian = TRUE
    )
    
    tmp <- tmp[is.finite(Wb) & Wb > 0]
    
    denom <- tmp[
      ,
      .(
        Wt = sum(get(weight_col)[treated == 1L], na.rm = TRUE),
        Wc = sum(get(weight_col)[treated == 0L], na.rm = TRUE)
      ),
      by = block_date
    ]
    
    tmp <- merge(tmp, denom, by = "block_date", all.x = TRUE)
    tmp <- tmp[is.finite(Wt) & Wt > 0 & is.finite(Wc) & Wc > 0]
    
    y <- tmp[[outcome_name]]
    w <- tmp[[weight_col]]
    
    tmp[
      ,
      `:=`(
        coef_y = fifelse(
          treated == 1L,
          -Wb * w / Wt,
          Wb * w / Wc
        ),
        coef_den = fifelse(
          treated == 1L,
          Wb * w / Wt,
          0
        ),
        y_value = y
      )
    ]
    
    contrib <- tmp[
      ,
      .(
        A = sum(coef_y * y_value, na.rm = TRUE),
        B = sum(coef_den, na.rm = TRUE),
        n_uses = .N
      ),
      by = .(outcome, obs_unit_id)
    ]
    
    out[[outcome_name]] <- contrib
  }
  
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

unit_multiplier_bootstrap_from_contrib <- function(contrib_dt,
                                                   outcomes,
                                                   B_boot,
                                                   seed) {
  set.seed(seed)
  
  out <- list()
  
  for (outcome_name in outcomes) {
    dt <- contrib_dt[outcome == outcome_name]
    N_star <- nrow(dt)
    
    boot_stat <- numeric(B_boot)
    
    for (b in seq_len(B_boot)) {
      xi <- rexp(N_star, rate = 1)
      
      denom <- sum(xi * dt$B, na.rm = TRUE)
      
      boot_stat[b] <- if (is.finite(denom) && denom > 0) {
        sum(xi * dt$A, na.rm = TRUE) / denom
      } else {
        NA_real_
      }
    }
    
    denom0 <- sum(dt$B, na.rm = TRUE)
    att0 <- sum(dt$A, na.rm = TRUE) / denom0
    
    out[[outcome_name]] <- data.table(
      estimator = "unit_multiplier_bootstrap",
      outcome = outcome_name,
      att = att0,
      std_error = sd(boot_stat, na.rm = TRUE),
      ci_low = as.numeric(quantile(boot_stat, 0.025, na.rm = TRUE)),
      ci_high = as.numeric(quantile(boot_stat, 0.975, na.rm = TRUE)),
      n = B_boot,
      n_units = N_star
    )
  }
  
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

# =============================================================================
# 5. ONE CBPS RUN
# =============================================================================

run_cbps_one <- function(DT, q, region_var, d, cfg) {
  message(
    sprintf(
      "\n[CBPS] q=%d | region=%s | temp=round(temp) | window=±%d days",
      q, region_var, d
    )
  )
  
  all_covariates <- get_all_covariates(region_var)
  
  work <- copy(DT)
  
  work[
    ,
    treated := as.integer(
      date >= cfg$treatment_start &
        date <= cfg$treatment_end &
        year == cfg$treatment_year
    )
  ]
  
  for (v in all_covariates) {
    work[, (v) := as.character(get(v))]
  }
  
  treatment_dates <- seq(cfg$treatment_start, cfg$treatment_end, by = "day")
  
  weighted_rows_list <- list()
  smd_block_list <- list()
  drop_log <- list()
  
  for (td in treatment_dates) {
    td <- as.Date(td)
    mdays <- make_month_day_window(td, d)
    
    Tb <- work[
      treated == 1L &
        date == td
    ]
    
    Cb <- work[
      treated == 0L &
        year != cfg$treatment_year &
        month_day %in% mdays
    ]
    
    if (nrow(Tb) == 0) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "no_treated"
      )
      next
    }
    
    if (nrow(Cb) == 0) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "no_controls"
      )
      next
    }
    
    pool <- rbindlist(
      list(Tb, Cb),
      use.names = TRUE,
      fill = TRUE
    )
    
    cc_vars <- unique(c("treated", cfg$outcomes, all_covariates))
    cc <- complete.cases(pool[, ..cc_vars])
    
    if (!all(cc)) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "complete_case_drop",
        n_drop = sum(!cc)
      )
    }
    
    pool <- pool[cc]
    
    if (nrow(pool) == 0 || length(unique(pool$treated)) < 2) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "insufficient_groups_after_complete_case"
      )
      next
    }
    
    for (v in all_covariates) {
      pool[, (v) := as.factor(get(v))]
    }
    
    # -------------------------------------------------------------------------
    # Block-specific CBPS model
    # -------------------------------------------------------------------------
    
    if (isTRUE(cfg$SCREEN_COVARIATES)) {
      screened <- screen_cbps_covariates(pool, all_covariates)
      cbps_covariates <- screened$keep
      
      if (nrow(screened$dropped) > 0) {
        tmp_drop <- copy(screened$dropped)
        tmp_drop[
          ,
          `:=`(
            block_date = td,
            reason = paste0("cbps_covariate_dropped_", reason)
          )
        ]
        drop_log[[length(drop_log) + 1L]] <- tmp_drop
      }
    } else {
      cbps_covariates <- all_covariates
    }
    
    if (length(cbps_covariates) == 0) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "no_usable_cbps_covariates"
      )
      next
    }
    
    cbps_formula <- as.formula(
      paste("treated ~", paste(cbps_covariates, collapse = " + "))
    )
    
    pool_df <- as.data.frame(pool)
    
    cbps_fit_obj <- fit_cbps_capture(
      formula = cbps_formula,
      data = pool_df,
      cfg = cfg
    )
    
    cbps_fit <- cbps_fit_obj$fit
    
    if (!isTRUE(cbps_fit_obj$ok)) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "cbps_fail",
        message = conditionMessage(cbps_fit),
        cbps_formula = deparse(cbps_formula)
      )
      next
    }
    
    if (length(cbps_fit_obj$warnings) > 0) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "cbps_warning",
        message = paste(cbps_fit_obj$warnings, collapse = " | "),
        cbps_formula = deparse(cbps_formula)
      )
    }
    
    ps <- tryCatch(
      extract_cbps_ps(cbps_fit, pool_df),
      error = function(e) e
    )
    
    if (inherits(ps, "error")) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "cbps_ps_extract_fail",
        message = conditionMessage(ps),
        cbps_formula = deparse(cbps_formula)
      )
      next
    }
    
    ps <- pmin(pmax(as.numeric(ps), cfg$ATT_TRIM_EPS), 1 - cfg$ATT_TRIM_EPS)
    
    pool[, ps := ps]
    
    # ATT-type CBPS weights:
    # treated: 1
    # control: e(X)/(1-e(X))
    pool[
      ,
      weights := fifelse(
        treated == 1L,
        1,
        ps / (1 - ps)
      )
    ]
    
    pool <- pool[is.finite(weights) & weights > 0]
    
    if (nrow(pool) == 0 || length(unique(pool$treated)) < 2) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "no_valid_cbps_weighted_rows"
      )
      next
    }
    
    pool[
      ,
      `:=`(
        block_date = td,
        q = q,
        region_var = region_var,
        temp_var = "temp_round",
        day_window = d
      )
    ]
    
    weighted_rows_list[[length(weighted_rows_list) + 1L]] <- pool
    
    block_weight_before <- nrow(pool[treated == 1L])
    block_weight_after <- sum(pool[treated == 1L]$weights, na.rm = TRUE)
    
    smd_block <- compute_block_smd(
      pool_before = rbindlist(list(Tb, Cb), use.names = TRUE, fill = TRUE),
      weighted_after = pool,
      covariates = all_covariates,
      block_date = td,
      block_weight_before = block_weight_before,
      block_weight_after = block_weight_after
    )
    
    if (nrow(smd_block) > 0) {
      smd_block_list[[length(smd_block_list) + 1L]] <- smd_block
    }
  }
  
  if (length(weighted_rows_list) == 0) {
    warning(sprintf("No valid CBPS weights: q=%d, region=%s, d=%d", q, region_var, d))
    return(NULL)
  }
  
  weighted_rows <- rbindlist(weighted_rows_list, use.names = TRUE, fill = TRUE)
  
  # ---------------------------------------------------------------------------
  # ATT and inference
  # ---------------------------------------------------------------------------
  
  delta_dt <- make_delta_table(
    rows = weighted_rows,
    outcomes = cfg$outcomes,
    weight_col = "weights"
  )
  
  block_att <- compute_block_att_from_delta(delta_dt)
  hajek_att <- hajek_from_block_table(block_att)
  
  cond_table <- conditional_se_from_delta(delta_dt)
  
  block_boot_table <- block_bootstrap_from_block_table(
    block_table = block_att,
    outcomes = cfg$outcomes,
    B_boot = cfg$B,
    seed = cfg$seed + q * 10000L + d * 100L +
      match(region_var, cfg$region_vars) * 10L + 1L
  )
  
  obs_contrib <- make_observation_contribution_table(
    rows = weighted_rows,
    block_table = block_att,
    outcomes = cfg$outcomes,
    weight_col = "weights"
  )
  
  mult_table <- unit_multiplier_bootstrap_from_contrib(
    contrib_dt = obs_contrib,
    outcomes = cfg$outcomes,
    B_boot = cfg$B,
    seed = cfg$seed + q * 10000L + d * 100L +
      match(region_var, cfg$region_vars) * 10L + 2L
  )
  
  att_summary <- rbindlist(
    list(cond_table, block_boot_table, mult_table),
    use.names = TRUE,
    fill = TRUE
  )
  
  att_summary[
    ,
    `:=`(
      method = "CBPS",
      q = q,
      region_var = region_var,
      temp_var = "temp_round",
      day_window = d,
      CBPS_METHOD = cfg$CBPS_METHOD,
      ATT_TRIM_EPS = cfg$ATT_TRIM_EPS,
      SCREEN_COVARIATES = cfg$SCREEN_COVARIATES,
      n_weighted_rows = nrow(weighted_rows),
      n_unique_treated_obs = uniqueN(weighted_rows[treated == 1L]$obs_unit_id),
      n_unique_control_obs = uniqueN(weighted_rows[treated == 0L]$obs_unit_id),
      n_treatment_date_blocks = uniqueN(weighted_rows$block_date)
    )
  ]
  
  setcolorder(
    att_summary,
    c(
      "method",
      "q",
      "region_var",
      "temp_var",
      "day_window",
      "estimator",
      "outcome",
      "att",
      "std_error",
      "ci_low",
      "ci_high",
      "n",
      "CBPS_METHOD",
      "ATT_TRIM_EPS",
      "SCREEN_COVARIATES",
      "n_weighted_rows",
      "n_unique_treated_obs",
      "n_unique_control_obs",
      "n_treatment_date_blocks"
    )
  )
  
  # ---------------------------------------------------------------------------
  # SMD aggregation
  # ---------------------------------------------------------------------------
  
  smd_block_variable <- if (length(smd_block_list) > 0) {
    rbindlist(smd_block_list, use.names = TRUE, fill = TRUE)
  } else {
    data.table()
  }
  
  smd_variable <- aggregate_smd_blocks_weighted_mean(smd_block_variable)
  
  if (nrow(smd_block_variable) > 0) {
    smd_block_variable[
      ,
      `:=`(
        method = "CBPS",
        q = q,
        region_var = region_var,
        temp_var = "temp_round",
        day_window = d
      )
    ]
  }
  
  if (nrow(smd_variable) > 0) {
    smd_variable[
      ,
      `:=`(
        method = "CBPS",
        q = q,
        region_var = region_var,
        temp_var = "temp_round",
        day_window = d,
        date_agg = "weighted_mean"
      )
    ]
  }
  
  smd_master_panel <- make_smd_master_for_panel(
    smd_variable = smd_variable,
    q = q,
    region_var = region_var,
    unit = "rt",
    d_win = d
  )
  
  smd_plot_summary <- copy(smd_variable)
  if (nrow(smd_plot_summary) > 0) {
    setnames(smd_plot_summary, "sample", "stage")
    setnames(smd_plot_summary, "max_abs_smd", "SMD")
  }
  
  drop_log_dt <- if (length(drop_log) > 0) {
    rbindlist(drop_log, use.names = TRUE, fill = TRUE)
  } else {
    data.table()
  }
  
  list(
    config = cfg,
    q = q,
    region_var = region_var,
    temp_var = "temp_round",
    day_window = d,
    all_covariates = all_covariates,
    weighted_rows = weighted_rows,
    delta_dt = delta_dt,
    block_att = block_att,
    hajek_att = hajek_att,
    obs_contrib = obs_contrib,
    att_summary = att_summary,
    smd_block_variable = smd_block_variable,
    smd_variable = smd_variable,
    smd_plot_summary = smd_plot_summary,
    smd_master_panel = smd_master_panel,
    drop_log = drop_log_dt
  )
}

# =============================================================================
# 6. SAVE FUNCTION
# =============================================================================

save_cbps_result <- function(res, cfg) {
  if (is.null(res)) return(invisible(NULL))
  
  out_dir <- file.path(
    cfg$root_out,
    sprintf("region_%s", res$region_var),
    sprintf("q%d", res$q)
  )
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  tag <- sprintf(
    "cbps_%s_rt_q%d_d%d",
    res$region_var,
    res$q,
    res$day_window
  )
  
  if (isTRUE(cfg$save_full_rds)) {
    saveRDS(
      res,
      file = file.path(out_dir, paste0(tag, "_full_result.rds"))
    )
  } else {
    light_res <- res
    light_res$weighted_rows <- NULL
    saveRDS(
      light_res,
      file = file.path(out_dir, paste0(tag, "_result_no_weighted_rows.rds"))
    )
  }
  
  if (isTRUE(cfg$save_weighted_rows_csv)) {
    fwrite(
      res$weighted_rows,
      file = file.path(out_dir, paste0(tag, "_weighted_rows.csv"))
    )
  }
  
  fwrite(res$delta_dt, file.path(out_dir, paste0(tag, "_delta_table.csv")))
  fwrite(res$block_att, file.path(out_dir, paste0(tag, "_block_att.csv")))
  fwrite(res$hajek_att, file.path(out_dir, paste0(tag, "_hajek_att.csv")))
  fwrite(res$obs_contrib, file.path(out_dir, paste0(tag, "_obs_contrib.csv")))
  fwrite(res$att_summary, file.path(out_dir, paste0(tag, "_att_summary.csv")))
  fwrite(res$smd_block_variable, file.path(out_dir, paste0(tag, "_smd_block_variable.csv")))
  fwrite(res$smd_variable, file.path(out_dir, paste0(tag, "_smd_variable_weighted_mean.csv")))
  fwrite(res$smd_plot_summary, file.path(out_dir, paste0(tag, "_smd_plot_summary.csv")))
  fwrite(res$smd_master_panel, file.path(out_dir, paste0(tag, "_smd_master_panel.csv")))
  fwrite(res$drop_log, file.path(out_dir, paste0(tag, "_drop_log.csv")))
  
  message(sprintf("[Saved] %s", out_dir))
}

# =============================================================================
# 7. MAIN LOOP
# =============================================================================

all_att <- list()
all_delta <- list()
all_block_att <- list()
all_hajek_att <- list()
all_obs_contrib <- list()
all_smd_block <- list()
all_smd_variable <- list()
all_smd_plot <- list()
all_smd_panel <- list()
all_drop_log <- list()

for (q in cfg$q_values) {
  message(sprintf("\n================ Loading q%d data ================", q))
  
  DT_raw <- load_q_data(q, cfg$root_in)
  DT <- prepare_data(DT_raw, q, cfg)
  
  for (region_var in cfg$region_vars) {
    if (!(region_var %in% names(DT))) {
      warning(sprintf("Skipping region_var=%s because it is not in the data.", region_var))
      next
    }
    
    for (d in cfg$day_windows) {
      res <- run_cbps_one(
        DT = DT,
        q = q,
        region_var = region_var,
        d = d,
        cfg = cfg
      )
      
      if (!is.null(res)) {
        save_cbps_result(res, cfg)
        
        key <- sprintf("%s_rt_q%d_d%d", region_var, q, d)
        
        all_att[[key]] <- res$att_summary
        all_delta[[key]] <- res$delta_dt
        all_block_att[[key]] <- res$block_att
        all_hajek_att[[key]] <- res$hajek_att
        all_obs_contrib[[key]] <- res$obs_contrib
        all_smd_block[[key]] <- res$smd_block_variable
        all_smd_variable[[key]] <- res$smd_variable
        all_smd_plot[[key]] <- res$smd_plot_summary
        all_smd_panel[[key]] <- res$smd_master_panel
        
        if (nrow(res$drop_log) > 0) {
          dl <- copy(res$drop_log)
          dl[
            ,
            `:=`(
              q = q,
              region_var = region_var,
              temp_var = "temp_round",
              day_window = d
            )
          ]
          all_drop_log[[key]] <- dl
        }
      }
    }
  }
}

# =============================================================================
# 8. GLOBAL SUMMARY FILES
# =============================================================================

if (length(all_att) > 0) {
  fwrite(
    rbindlist(all_att, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "cbps_all_att_summary.csv")
  )
}

if (length(all_delta) > 0) {
  fwrite(
    rbindlist(all_delta, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "cbps_all_delta_table.csv")
  )
}

if (length(all_block_att) > 0) {
  fwrite(
    rbindlist(all_block_att, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "cbps_all_block_att_summary.csv")
  )
}

if (length(all_hajek_att) > 0) {
  fwrite(
    rbindlist(all_hajek_att, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "cbps_all_hajek_att_summary.csv")
  )
}

if (length(all_obs_contrib) > 0) {
  fwrite(
    rbindlist(all_obs_contrib, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "cbps_all_obs_contrib.csv")
  )
}

if (length(all_smd_block) > 0) {
  fwrite(
    rbindlist(all_smd_block, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "cbps_all_smd_block_variable_summary.csv")
  )
}

if (length(all_smd_variable) > 0) {
  fwrite(
    rbindlist(all_smd_variable, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "cbps_all_smd_variable_summary.csv")
  )
}

if (length(all_smd_plot) > 0) {
  fwrite(
    rbindlist(all_smd_plot, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "cbps_all_smd_plot_summary.csv")
  )
}

if (length(all_smd_panel) > 0) {
  fwrite(
    rbindlist(all_smd_panel, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "CBPS_SMD_values_master_rt_weighted_mean.csv")
  )
}

if (length(all_drop_log) > 0) {
  fwrite(
    rbindlist(all_drop_log, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "cbps_all_drop_log.csv")
  )
}

message("\nCBPS weighting + block-wise ATT + bootstrap SE + SMD diagnostics completed.")
