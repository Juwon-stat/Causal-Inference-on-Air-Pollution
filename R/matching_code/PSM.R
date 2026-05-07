###############################################################################
# Title   : PSM_GLOBAL_FINAL.R
# Purpose : Block-wise Propensity Score Matching (PSM)
#           + block-specific ATT
#           + Hajek aggregation
#           + conditional SE
#           + block bootstrap
#           + unit-level multiplier bootstrap
#           + SMD diagnostics for before-vs-after balance plots
#
# ATT direction:
#   ATT = matched control mean - treated mean
#   Positive ATT means PM concentration was lower during the lockdown period.
#
# Design:
#   - Observational unit: district-date-hour record
#   - Treatment-date block: one block per treated calendar date
#   - Control pool: non-2020 observations within ±d calendar days
#   - PS model is fitted within each block only
#   - Matching weights are treated as fixed for inference
#
# Output:
#   - psm_all_att_summary.csv
#   - psm_all_block_att_summary.csv
#   - psm_all_smd_variable_summary.csv
#   - psm_all_smd_plot_summary.csv
#   - PSM_SMD_values_master_rt_weighted_mean.csv
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(MatchIt)
})

# =============================================================================
# 0. USER CONFIGURATION
# =============================================================================

cfg <- list(
  root_in  = "C:/path/to/Seoul_RData",
  root_out = "C:/path/to/PSM_results",
  
  q_values = c(3, 4, 5),
  day_windows = c(3, 7, 14),
  region_vars = c("dist", "area"),
  
  treatment_start = as.Date("2020-01-24"),
  treatment_end   = as.Date("2020-02-09"),
  treatment_year  = 2020L,
  
  outcomes = c("PM10", "PM25"),
  
  # PSM settings
  K_MAX = 3L,
  CALIPER_SD = 0.2,
  REPLACE_CTL = TRUE,
  
  # Stabilization settings
  # These reflect the previous PSM implementation logic:
  # exact matching on spatial/weekday/time/temperature structure,
  # and PS modeling on the remaining covariates to reduce separation.
  USE_EXACT_IN_PSM = TRUE,
  exact_base_vars = c("wday", "time_of_day", "temp_round"),
  USE_REDUCED_PS_MODEL = TRUE,
  
  # Bootstrap
  B = 1000L,
  seed = 20260507L,
  
  # Save options
  save_full_rds = TRUE,
  save_matched_rows_csv = FALSE
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

logit_safe <- function(p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  qlogis(p)
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

get_exact_covariates <- function(region_var, cfg) {
  if (!isTRUE(cfg$USE_EXACT_IN_PSM)) return(character(0))
  unique(c(region_var, cfg$exact_base_vars))
}

fit_glm_capture <- function(formula, data) {
  warning_msgs <- character(0)
  
  fit <- withCallingHandlers(
    tryCatch(
      glm(
        formula,
        data = data,
        family = binomial(),
        control = glm.control(maxit = 50)
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
    converged = if (!inherits(fit, "error")) isTRUE(fit$converged) else FALSE
  )
}

screen_ps_covariates <- function(pool, candidate_covariates) {
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
  
  # Temperature matching/PS covariate
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
                              matched_after,
                              covariates,
                              block_date,
                              block_weight_before,
                              block_weight_after) {
  # BEFORE: eligible sample in the block before matching
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
  
  # AFTER: matched weighted sample from MatchIt
  after_treated <- matched_after[treated == 1L, c(covariates, "weights"), with = FALSE]
  after_control <- matched_after[treated == 0L, c(covariates, "weights"), with = FALSE]
  
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
  # This table implements the paper's Delta_bu structure:
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
  # Implements:
  # Var_cond(tau_hat) = sum_b a_b^2 * (1/W_b^2) *
  #                     sum_u omega_bu^2 (Delta_bu - tau_b)^2.
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
  # Resample treatment-date blocks, as described in the user's block-wise design.
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
    
    # Need within-block treated/control denominators.
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
# 5. ONE PSM RUN
# =============================================================================

run_psm_one <- function(DT, q, region_var, d, cfg) {
  message(
    sprintf(
      "\n[PSM] q=%d | region=%s | temp=round(temp) | window=±%d days",
      q, region_var, d
    )
  )
  
  all_covariates <- get_all_covariates(region_var)
  exact_covariates <- get_exact_covariates(region_var, cfg)
  
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
  
  matched_rows_list <- list()
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
    # Block-specific propensity score model
    # -------------------------------------------------------------------------
    
    if (isTRUE(cfg$USE_REDUCED_PS_MODEL)) {
      ps_model_candidates <- setdiff(all_covariates, exact_covariates)
    } else {
      ps_model_candidates <- all_covariates
    }
    
    screened <- screen_ps_covariates(pool, ps_model_candidates)
    usable_covariates <- screened$keep
    
    if (nrow(screened$dropped) > 0) {
      tmp_drop <- copy(screened$dropped)
      tmp_drop[
        ,
        `:=`(
          block_date = td,
          reason = paste0("ps_covariate_dropped_", reason)
        )
      ]
      drop_log[[length(drop_log) + 1L]] <- tmp_drop
    }
    
    if (length(usable_covariates) == 0) {
      ps_formula <- treated ~ 1
    } else {
      ps_formula <- as.formula(
        paste("treated ~", paste(usable_covariates, collapse = " + "))
      )
    }
    
    ps_fit_obj <- fit_glm_capture(ps_formula, data = pool)
    ps_fit <- ps_fit_obj$fit
    
    if (inherits(ps_fit, "error")) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "glm_fail",
        message = conditionMessage(ps_fit),
        ps_formula = deparse(ps_formula)
      )
      next
    }
    
    if (length(ps_fit_obj$warnings) > 0 || !isTRUE(ps_fit_obj$converged)) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "glm_warning_or_nonconvergence",
        message = paste(ps_fit_obj$warnings, collapse = " | "),
        converged = ps_fit_obj$converged,
        ps_formula = deparse(ps_formula)
      )
    }
    
    pool[, ps := pmin(pmax(as.numeric(fitted(ps_fit)), 1e-6), 1 - 1e-6)]
    pool[, lp := logit_safe(ps)]
    
    sd_lp <- sd(pool$lp, na.rm = TRUE)
    
    if (!is.finite(sd_lp) || sd_lp <= 0) {
      pool[, lp := 0]
      caliper_value <- NULL
      
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "constant_or_invalid_logit_ps_use_no_caliper",
        ps_formula = deparse(ps_formula)
      )
    } else {
      caliper_value <- cfg$CALIPER_SD * sd_lp
    }
    
    # -------------------------------------------------------------------------
    # MatchIt nearest-neighbor PSM within the block
    # -------------------------------------------------------------------------
    
    pool_df <- as.data.frame(pool)
    rownames(pool_df) <- as.character(pool_df$row_id)
    
    exact_use <- intersect(exact_covariates, names(pool_df))
    exact_use <- exact_use[
      vapply(
        exact_use,
        function(v) length(unique(pool_df[[v]][!is.na(pool_df[[v]])])) > 1,
        logical(1)
      )
    ]
    
    m <- tryCatch(
      {
        if (isTRUE(cfg$USE_EXACT_IN_PSM) && length(exact_use) > 0) {
          matchit(
            formula = treated ~ 1,
            data = pool_df,
            method = "nearest",
            distance = pool_df$lp,
            replace = cfg$REPLACE_CTL,
            ratio = cfg$K_MAX,
            caliper = caliper_value,
            std.caliper = FALSE,
            exact = exact_use,
            estimand = "ATT"
          )
        } else {
          matchit(
            formula = treated ~ 1,
            data = pool_df,
            method = "nearest",
            distance = pool_df$lp,
            replace = cfg$REPLACE_CTL,
            ratio = cfg$K_MAX,
            caliper = caliper_value,
            std.caliper = FALSE,
            estimand = "ATT"
          )
        }
      },
      error = function(e) e
    )
    
    if (inherits(m, "error")) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "matchit_fail",
        message = conditionMessage(m),
        ps_formula = deparse(ps_formula)
      )
      next
    }
    
    md <- as.data.table(match.data(m, data = pool_df, drop.unmatched = TRUE))
    
    if (!("weights" %in% names(md))) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "weights_missing"
      )
      next
    }
    
    md <- md[is.finite(weights) & weights > 0]
    
    if (nrow(md) == 0 || length(unique(md$treated)) < 2) {
      drop_log[[length(drop_log) + 1L]] <- data.table(
        block_date = td,
        reason = "no_valid_matched_rows"
      )
      next
    }
    
    md[
      ,
      `:=`(
        block_date = td,
        q = q,
        region_var = region_var,
        temp_var = "temp_round",
        day_window = d
      )
    ]
    
    matched_rows_list[[length(matched_rows_list) + 1L]] <- md
    
    block_weight_before <- nrow(pool[treated == 1L])
    block_weight_after <- sum(md[treated == 1L]$weights, na.rm = TRUE)
    
    smd_block <- compute_block_smd(
      pool_before = pool,
      matched_after = md,
      covariates = all_covariates,
      block_date = td,
      block_weight_before = block_weight_before,
      block_weight_after = block_weight_after
    )
    
    if (nrow(smd_block) > 0) {
      smd_block_list[[length(smd_block_list) + 1L]] <- smd_block
    }
  }
  
  if (length(matched_rows_list) == 0) {
    warning(sprintf("No valid PSM matches: q=%d, region=%s, d=%d", q, region_var, d))
    return(NULL)
  }
  
  matched_rows <- rbindlist(matched_rows_list, use.names = TRUE, fill = TRUE)
  
  # ---------------------------------------------------------------------------
  # ATT and inference
  # ---------------------------------------------------------------------------
  
  delta_dt <- make_delta_table(
    rows = matched_rows,
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
    rows = matched_rows,
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
      method = "PSM",
      q = q,
      region_var = region_var,
      temp_var = "temp_round",
      day_window = d,
      K_MAX = cfg$K_MAX,
      CALIPER_SD = cfg$CALIPER_SD,
      REPLACE_CTL = cfg$REPLACE_CTL,
      USE_EXACT_IN_PSM = cfg$USE_EXACT_IN_PSM,
      USE_REDUCED_PS_MODEL = cfg$USE_REDUCED_PS_MODEL,
      n_weighted_rows = nrow(matched_rows),
      n_unique_treated_obs = uniqueN(matched_rows[treated == 1L]$obs_unit_id),
      n_unique_control_obs = uniqueN(matched_rows[treated == 0L]$obs_unit_id),
      n_treatment_date_blocks = uniqueN(matched_rows$block_date)
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
      "K_MAX",
      "CALIPER_SD",
      "REPLACE_CTL",
      "USE_EXACT_IN_PSM",
      "USE_REDUCED_PS_MODEL",
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
        method = "PSM",
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
        method = "PSM",
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
    exact_covariates = exact_covariates,
    matched_rows = matched_rows,
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

save_psm_result <- function(res, cfg) {
  if (is.null(res)) return(invisible(NULL))
  
  out_dir <- file.path(
    cfg$root_out,
    sprintf("region_%s", res$region_var),
    sprintf("q%d", res$q)
  )
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  tag <- sprintf(
    "psm_%s_rt_q%d_d%d",
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
    light_res$matched_rows <- NULL
    saveRDS(
      light_res,
      file = file.path(out_dir, paste0(tag, "_result_no_matched_rows.rds"))
    )
  }
  
  if (isTRUE(cfg$save_matched_rows_csv)) {
    fwrite(
      res$matched_rows,
      file = file.path(out_dir, paste0(tag, "_matched_weighted_rows.csv"))
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
      res <- run_psm_one(
        DT = DT,
        q = q,
        region_var = region_var,
        d = d,
        cfg = cfg
      )
      
      if (!is.null(res)) {
        save_psm_result(res, cfg)
        
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
    file.path(cfg$root_out, "psm_all_att_summary.csv")
  )
}

if (length(all_delta) > 0) {
  fwrite(
    rbindlist(all_delta, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "psm_all_delta_table.csv")
  )
}

if (length(all_block_att) > 0) {
  fwrite(
    rbindlist(all_block_att, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "psm_all_block_att_summary.csv")
  )
}

if (length(all_hajek_att) > 0) {
  fwrite(
    rbindlist(all_hajek_att, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "psm_all_hajek_att_summary.csv")
  )
}

if (length(all_obs_contrib) > 0) {
  fwrite(
    rbindlist(all_obs_contrib, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "psm_all_obs_contrib.csv")
  )
}

if (length(all_smd_block) > 0) {
  fwrite(
    rbindlist(all_smd_block, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "psm_all_smd_block_variable_summary.csv")
  )
}

if (length(all_smd_variable) > 0) {
  fwrite(
    rbindlist(all_smd_variable, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "psm_all_smd_variable_summary.csv")
  )
}

if (length(all_smd_plot) > 0) {
  fwrite(
    rbindlist(all_smd_plot, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "psm_all_smd_plot_summary.csv")
  )
}

if (length(all_smd_panel) > 0) {
  fwrite(
    rbindlist(all_smd_panel, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "PSM_SMD_values_master_rt_weighted_mean.csv")
  )
}

if (length(all_drop_log) > 0) {
  fwrite(
    rbindlist(all_drop_log, use.names = TRUE, fill = TRUE),
    file.path(cfg$root_out, "psm_all_drop_log.csv")
  )
}

message("\nPSM matching + block-wise ATT + bootstrap SE + SMD diagnostics completed.")
