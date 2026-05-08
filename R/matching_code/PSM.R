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
  # Folder containing seoul_q3.RData, seoul_q4.RData, seoul_q5.RData
  root_in  = "C:/path/to/Seoul_RData",
  
  # Output folder
  root_out = "C:/path/to/PSM_results",
  q_values = c(3, 4, 5),
  day_windows = c(3, 7, 14),
  region_vars = c("dist", "area"),
  treatment_start = as.Date("2020-01-24"),
  treatment_end   = as.Date("2020-02-09"),
  treatment_year  = 2020L,
  outcomes = c("PM10", "PM25"),
  K_MAX = 3L,
  CALIPER_SD = 0.2,
  REPLACE_CTL = TRUE,
  exact_base_vars = c("wday", "time_of_day", "temp_round"),
  USE_REDUCED_PS_MODEL = TRUE,
  B = 1000L,
  seed = 20260508L,
  save_full_rds = FALSE,
  save_matched_rows_csv = FALSE
)

dir.create(cfg$root_out, recursive = TRUE, showWarnings = FALSE)

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

wmean <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

wvar <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  mu <- wmean(x[ok], w[ok])
  sum(w[ok] * (x[ok] - mu)^2) / sum(w[ok])
}

smd_bin <- function(x_t, w_t, x_c, w_c) {
  mt <- wmean(x_t, w_t)
  mc <- wmean(x_c, w_c)
  vt <- wvar(x_t, w_t)
  vc <- wvar(x_c, w_c)
  sp <- sqrt((vt + vc) / 2)
  
  if (!is.finite(sp) || sp == 0) {
    if (is.finite(mt) && is.finite(mc) && abs(mt - mc) < .Machine$double.eps) return(0)
    return(NA_real_)
  }
  
  abs(mt - mc) / sp
}

logit_safe <- function(p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  qlogis(p)
}

date_window_md <- function(date_i, d) {
  format(seq(date_i - d, date_i + d, by = "day"), "%m-%d")
}

get_covs <- function(region_var) {
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

load_q_data <- function(q, root_in) {
  f <- file.path(root_in, sprintf("seoul_q%d.RData", q))
  e <- new.env()
  load(f, envir = e)
  as.data.table(get(sprintf("seoul_q%d", q), envir = e))
}

prep_data <- function(DT) {
  DT <- copy(DT)
  DT[, date := as.Date(date)]
  if (!("year" %in% names(DT))) DT[, year := as.integer(format(date, "%Y"))]
  DT[, row_id := .I]
  DT[, month_day := format(date, "%m-%d")]
  DT[, temp_round := round(temp)]
  DT[, obs_unit_id := paste(dist, as.character(date), hour, sep = "_")]
  for (y in cfg$outcomes) DT[, (y) := to_num(get(y))]
  DT[]
}

screen_covs <- function(pool, covs) {
  keep <- character(0)
  drop <- list()
  
  for (v in covs) {
    x <- pool[[v]]
    z <- pool$treated
    
    if (length(unique(x[!is.na(x)])) < 2) {
      drop[[length(drop) + 1L]] <- data.table(covariate = v, reason = "single_level")
      next
    }
    
    tab <- table(z, x, useNA = "no")
    
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      drop[[length(drop) + 1L]] <- data.table(covariate = v, reason = "insufficient_table")
      next
    }
    
    both <- apply(tab, 2, function(a) all(a > 0))
    
    if (!all(both)) {
      drop[[length(drop) + 1L]] <- data.table(covariate = v, reason = "level_not_shared")
      next
    }
    
    keep <- c(keep, v)
  }
  
  list(
    keep = unique(keep),
    drop = if (length(drop)) rbindlist(drop, fill = TRUE) else data.table()
  )
}

fit_glm_capture <- function(formula, data) {
  warns <- character(0)
  
  fit <- withCallingHandlers(
    tryCatch(
      glm(formula, data = data, family = binomial(), control = glm.control(maxit = 50)),
      error = function(e) e
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  list(
    fit = fit,
    warnings = unique(warns),
    converged = if (!inherits(fit, "error")) isTRUE(fit$converged) else FALSE
  )
}

make_smd_pop <- function(rows, covs, weight_col, stage) {
  if (!nrow(rows)) return(data.table())
  
  dt <- copy(rows[, c("block_date", "treated", covs, weight_col), with = FALSE])
  setnames(dt, weight_col, "local_weight")
  dt[, local_weight := as.numeric(local_weight)]
  dt <- dt[is.finite(local_weight) & local_weight > 0]
  
  den <- dt[, .(Wg = sum(local_weight)), by = .(block_date, treated)]
  bw  <- den[treated == 1L, .(block_date, Wb = Wg)]
  
  dt <- merge(dt, den, by = c("block_date", "treated"), all.x = TRUE)
  dt <- merge(dt, bw, by = "block_date", all.x = TRUE)
  
  dt <- dt[is.finite(Wg) & Wg > 0 & is.finite(Wb) & Wb > 0]
  dt[, smd_weight := Wb * local_weight / Wg]
  dt[, stage := stage]
  dt[, c("local_weight", "Wg", "Wb") := NULL]
  dt[]
}

pooled_smd <- function(before_rows, after_rows, covs) {
  before_pop <- make_smd_pop(before_rows, covs, "before_weight", "before")
  after_pop  <- make_smd_pop(after_rows,  covs, "weights",        "after")
  
  pop <- rbindlist(list(before_pop, after_pop), fill = TRUE)
  
  lev_out <- list()
  
  for (st in c("before", "after")) {
    dts <- pop[stage == st]
    tr <- dts[treated == 1L]
    co <- dts[treated == 0L]
    
    if (!nrow(tr) || !nrow(co)) next
    
    for (v in covs) {
      xt <- as.character(tr[[v]])
      xc <- as.character(co[[v]])
      wt <- tr$smd_weight
      wc <- co$smd_weight
      lvls <- sort(unique(c(xt, xc)))
      lvls <- lvls[!is.na(lvls)]
      
      if (!length(lvls)) next
      
      tmp <- rbindlist(lapply(lvls, function(lv) {
        zt <- as.numeric(xt == lv)
        zc <- as.numeric(xc == lv)
        
        data.table(
          stage = st,
          covariate = v,
          level = lv,
          smd = smd_bin(zt, wt, zc, wc),
          treated_mean = wmean(zt, wt),
          control_mean = wmean(zc, wc),
          treated_weight_sum = sum(wt[is.finite(wt) & wt > 0]),
          control_weight_sum = sum(wc[is.finite(wc) & wc > 0])
        )
      }), fill = TRUE)
      
      lev_out[[paste(st, v, sep = "_")]] <- tmp
    }
  }
  
  level_dt <- if (length(lev_out)) rbindlist(lev_out, fill = TRUE) else data.table()
  
  var_dt <- level_dt[
    ,
    .(
      max_abs_smd = max(smd, na.rm = TRUE),
      mean_abs_smd = mean(smd, na.rm = TRUE),
      n_levels = .N
    ),
    by = .(stage, covariate)
  ]
  
  var_dt[is.infinite(max_abs_smd), max_abs_smd := NA_real_]
  var_dt[is.infinite(mean_abs_smd), mean_abs_smd := NA_real_]
  
  list(population = pop, level = level_dt, variable = var_dt)
}

make_panel_master <- function(smd_var, q, region_var, d) {
  if (!nrow(smd_var)) return(data.table())
  
  dt <- copy(smd_var)
  dt[, variable := covariate]
  dt[, q := paste0("q", q)]
  dt[, region := region_var]
  dt[, unit := "rt"]
  dt[, d_win := as.integer(d)]
  dt[, val := max_abs_smd]
  
  wide <- dcast(
    dt,
    variable + region + unit + q + d_win ~ stage,
    value.var = "val",
    fun.aggregate = function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      max(x)
    }
  )
  
  if (!("before" %in% names(wide))) wide[, before := NA_real_]
  if (!("after" %in% names(wide))) wide[, after := NA_real_]
  
  wide[, .(variable, before, after, region, unit, q, d_win)]
}

make_delta <- function(rows, outcomes) {
  out <- list()
  
  for (y in outcomes) {
    tmp <- rows[
      ,
      {
        tr <- .SD[treated == 1L]
        co <- .SD[treated == 0L]
        
        wt <- tr$weights
        wc <- co$weights
        
        yt <- tr[[y]]
        yc <- co[[y]]
        
        mu_c <- wmean(yc, wc)
        mu_t <- wmean(yt, wt)
        
        Wb <- sum(wt[is.finite(wt) & wt > 0])
        tau_b <- mu_c - mu_t
        
        data.table(
          outcome = y,
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
    
    out[[y]] <- tmp
  }
  
  rbindlist(out, fill = TRUE)
}

block_att_from_delta <- function(delta) {
  delta[
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
}

hajek_att <- function(block_att) {
  block_att[
    is.finite(Wb) & Wb > 0 & is.finite(numerator),
    .(
      att = sum(numerator) / sum(Wb),
      W_total = sum(Wb),
      n_blocks = uniqueN(block_date)
    ),
    by = outcome
  ]
}

conditional_se <- function(delta) {
  bv <- delta[
    is.finite(omega_bu) & omega_bu > 0 & is.finite(Delta_bu) & is.finite(Wb) & Wb > 0,
    .(
      tau_b = tau_b[1],
      Wb = Wb[1],
      var_delta_b = sum(omega_bu^2 * (Delta_bu - tau_b[1])^2) / Wb[1]^2,
      n_treated = .N
    ),
    by = .(outcome, block_date)
  ]
  
  wt <- bv[, .(W_total = sum(Wb)), by = outcome]
  bv <- merge(bv, wt, by = "outcome")
  bv[, a_b := Wb / W_total]
  
  bv[
    ,
    {
      att0 <- sum(Wb * tau_b) / sum(Wb)
      se0 <- sqrt(sum(a_b^2 * var_delta_b))
      .(
        estimator = "conditional",
        att = att0,
        std_error = se0,
        ci_low = att0 - 1.96 * se0,
        ci_high = att0 + 1.96 * se0,
        n = sum(n_treated),
        n_blocks = uniqueN(block_date)
      )
    },
    by = outcome
  ]
}

block_boot <- function(block_att, outcomes, B, seed) {
  set.seed(seed)
  
  blocks <- sort(unique(block_att$block_date))
  K <- length(blocks)
  out <- list()
  
  for (y in outcomes) {
    bt <- block_att[outcome == y]
    stat <- numeric(B)
    
    for (b in seq_len(B)) {
      samp <- sample(blocks, K, replace = TRUE)
      kap <- data.table(block_date = samp)[, .N, by = block_date]
      setnames(kap, "N", "kappa")
      bb <- merge(bt, kap, by = "block_date")
      den <- bb[, sum(kappa * Wb)]
      
      stat[b] <- if (is.finite(den) && den > 0) {
        bb[, sum(kappa * numerator) / den]
      } else {
        NA_real_
      }
    }
    
    att0 <- hajek_att(bt)$att[1]
    
    out[[y]] <- data.table(
      estimator = "block_bootstrap",
      outcome = y,
      att = att0,
      std_error = sd(stat, na.rm = TRUE),
      ci_low = as.numeric(quantile(stat, 0.025, na.rm = TRUE)),
      ci_high = as.numeric(quantile(stat, 0.975, na.rm = TRUE)),
      n = B,
      n_blocks = K
    )
  }
  
  rbindlist(out, fill = TRUE)
}

obs_contrib_table <- function(rows, block_att, outcomes) {
  out <- list()
  
  for (y in outcomes) {
    bt <- block_att[outcome == y, .(block_date, outcome, Wb)]
    tmp <- merge(rows, bt, by = "block_date", allow.cartesian = TRUE)
    
    den <- tmp[
      ,
      .(
        Wt = sum(weights[treated == 1L], na.rm = TRUE),
        Wc = sum(weights[treated == 0L], na.rm = TRUE)
      ),
      by = block_date
    ]
    
    tmp <- merge(tmp, den, by = "block_date")
    tmp <- tmp[is.finite(Wt) & Wt > 0 & is.finite(Wc) & Wc > 0]
    
    yy <- tmp[[y]]
    ww <- tmp$weights
    
    tmp[
      ,
      `:=`(
        coef_y = fifelse(treated == 1L, -Wb * ww / Wt, Wb * ww / Wc),
        coef_den = fifelse(treated == 1L, Wb * ww / Wt, 0),
        y_value = yy
      )
    ]
    
    out[[y]] <- tmp[
      ,
      .(
        A = sum(coef_y * y_value, na.rm = TRUE),
        B = sum(coef_den, na.rm = TRUE),
        n_uses = .N
      ),
      by = .(outcome, obs_unit_id)
    ]
  }
  
  rbindlist(out, fill = TRUE)
}

unit_multiplier <- function(contrib, outcomes, B, seed) {
  set.seed(seed)
  
  out <- list()
  
  for (y in outcomes) {
    dt <- contrib[outcome == y]
    n <- nrow(dt)
    stat <- numeric(B)
    
    for (b in seq_len(B)) {
      xi <- rexp(n)
      den <- sum(xi * dt$B)
      
      stat[b] <- if (is.finite(den) && den > 0) {
        sum(xi * dt$A) / den
      } else {
        NA_real_
      }
    }
    
    den0 <- sum(dt$B)
    att0 <- sum(dt$A) / den0
    
    out[[y]] <- data.table(
      estimator = "unit_multiplier_bootstrap",
      outcome = y,
      att = att0,
      std_error = sd(stat, na.rm = TRUE),
      ci_low = as.numeric(quantile(stat, 0.025, na.rm = TRUE)),
      ci_high = as.numeric(quantile(stat, 0.975, na.rm = TRUE)),
      n = B,
      n_units = n
    )
  }
  
  rbindlist(out, fill = TRUE)
}

run_psm_one <- function(DT, q, region_var, d) {
  message(sprintf("[PSM] q=%d | region=%s | d=%d", q, region_var, d))
  
  covs <- get_covs(region_var)
  exact_covs <- unique(c(region_var, cfg$exact_base_vars))
  
  work <- copy(DT)
  work[
    ,
    treated := as.integer(
      date >= cfg$treatment_start &
        date <= cfg$treatment_end &
        year == cfg$treatment_year
    )
  ]
  
  for (v in covs) work[, (v) := as.character(get(v))]
  
  dates <- seq(cfg$treatment_start, cfg$treatment_end, by = "day")
  
  before_list <- list()
  match_list <- list()
  drop_list <- list()
  
  for (td in dates) {
    td <- as.Date(td, origin = "1970-01-01")
    mdays <- date_window_md(td, d)
    
    Tb <- work[treated == 1L & date == td]
    Cb <- work[treated == 0L & year != cfg$treatment_year & month_day %in% mdays]
    
    if (!nrow(Tb) || !nrow(Cb)) {
      drop_list[[length(drop_list) + 1L]] <- data.table(block_date = td, reason = "no_treated_or_control")
      next
    }
    
    pool <- rbindlist(list(Tb, Cb), fill = TRUE)
    
    cc_vars <- unique(c("treated", cfg$outcomes, covs))
    pool <- pool[complete.cases(pool[, ..cc_vars])]
    
    if (!nrow(pool) || length(unique(pool$treated)) < 2) {
      drop_list[[length(drop_list) + 1L]] <- data.table(block_date = td, reason = "insufficient_after_cc")
      next
    }
    
    for (v in covs) pool[, (v) := factor(get(v))]
    
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
    
    ps_covs <- if (cfg$USE_REDUCED_PS_MODEL) setdiff(covs, exact_covs) else covs
    sc <- screen_covs(pool, ps_covs)
    use_covs <- sc$keep
    
    if (nrow(sc$drop)) {
      dd <- copy(sc$drop)
      dd[, block_date := td]
      drop_list[[length(drop_list) + 1L]] <- dd
    }
    
    fml <- if (length(use_covs)) {
      as.formula(paste("treated ~", paste(use_covs, collapse = " + ")))
    } else {
      treated ~ 1
    }
    
    fit_obj <- fit_glm_capture(fml, pool)
    fit <- fit_obj$fit
    
    if (inherits(fit, "error")) {
      drop_list[[length(drop_list) + 1L]] <- data.table(
        block_date = td,
        reason = "glm_fail",
        message = conditionMessage(fit)
      )
      next
    }
    
    pool[, ps := pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6)]
    pool[, lp := logit_safe(ps)]
    
    sdlp <- sd(pool$lp, na.rm = TRUE)
    cal <- if (is.finite(sdlp) && sdlp > 0) cfg$CALIPER_SD * sdlp else NULL
    
    pool_df <- as.data.frame(pool)
    rownames(pool_df) <- as.character(pool_df$row_id)
    
    exact_use <- intersect(exact_covs, names(pool_df))
    exact_use <- exact_use[
      vapply(
        exact_use,
        function(v) length(unique(pool_df[[v]][!is.na(pool_df[[v]])])) > 1,
        logical(1)
      )
    ]
    
    m <- tryCatch(
      matchit(
        formula = treated ~ 1,
        data = pool_df,
        method = "nearest",
        distance = pool_df$lp,
        replace = cfg$REPLACE_CTL,
        ratio = cfg$K_MAX,
        caliper = cal,
        std.caliper = FALSE,
        exact = exact_use,
        estimand = "ATT"
      ),
      error = function(e) e
    )
    
    if (inherits(m, "error")) {
      drop_list[[length(drop_list) + 1L]] <- data.table(
        block_date = td,
        reason = "matchit_fail",
        message = conditionMessage(m)
      )
      next
    }
    
    md <- as.data.table(match.data(m, data = pool_df, drop.unmatched = TRUE))
    md <- md[is.finite(weights) & weights > 0]
    
    if (!nrow(md) || length(unique(md$treated)) < 2) {
      drop_list[[length(drop_list) + 1L]] <- data.table(block_date = td, reason = "empty_matched")
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
    
    pool[, before_weight := 1]
    
    before_list[[length(before_list) + 1L]] <- pool
    match_list[[length(match_list) + 1L]] <- md
  }
  
  if (!length(match_list)) return(NULL)
  
  before_rows <- rbindlist(before_list, fill = TRUE)
  matched_rows <- rbindlist(match_list, fill = TRUE)
  
  delta <- make_delta(matched_rows, cfg$outcomes)
  batt <- block_att_from_delta(delta)
  hajek <- hajek_att(batt)
  cond <- conditional_se(delta)
  bboot <- block_boot(
    batt,
    cfg$outcomes,
    cfg$B,
    cfg$seed + q * 10000L + d * 100L + match(region_var, cfg$region_vars) * 10L + 1L
  )
  
  contrib <- obs_contrib_table(matched_rows, batt, cfg$outcomes)
  
  mult <- unit_multiplier(
    contrib,
    cfg$outcomes,
    cfg$B,
    cfg$seed + q * 10000L + d * 100L + match(region_var, cfg$region_vars) * 10L + 2L
  )
  
  att <- rbindlist(list(cond, bboot, mult), fill = TRUE)
  
  att[
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
      exact_covariates = paste(exact_covs, collapse = " + "),
      n_weighted_rows = nrow(matched_rows),
      n_unique_treated_obs = uniqueN(matched_rows[treated == 1L]$obs_unit_id),
      n_unique_control_obs = uniqueN(matched_rows[treated == 0L]$obs_unit_id),
      n_treatment_date_blocks = uniqueN(matched_rows$block_date)
    )
  ]
  
  smd <- pooled_smd(before_rows, matched_rows, covs)
  
  for (nm in c("level", "variable")) {
    smd[[nm]][
      ,
      `:=`(
        method = "PSM",
        q = q,
        region_var = region_var,
        temp_var = "temp_round",
        day_window = d,
        smd_type = "pooled_hajek"
      )
    ]
  }
  
  smd_plot <- copy(smd$variable)
  if (nrow(smd_plot)) setnames(smd_plot, "max_abs_smd", "SMD")
  
  panel <- make_panel_master(smd$variable, q, region_var, d)
  
  drop <- if (length(drop_list)) rbindlist(drop_list, fill = TRUE) else data.table()
  
  list(
    before_rows = before_rows,
    matched_rows = matched_rows,
    delta = delta,
    block_att = batt,
    hajek = hajek,
    contrib = contrib,
    att = att,
    smd_level = smd$level,
    smd_variable = smd$variable,
    smd_plot = smd_plot,
    smd_panel = panel,
    drop = drop,
    q = q,
    region_var = region_var,
    d = d
  )
}

save_one <- function(res) {
  if (is.null(res)) return(invisible(NULL))
  
  od <- file.path(
    cfg$root_out,
    sprintf("region_%s", res$region_var),
    sprintf("q%d", res$q)
  )
  
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  
  tag <- sprintf("psm_%s_rt_q%d_d%d", res$region_var, res$q, res$d)
  
  if (cfg$save_full_rds) {
    tmp <- res
    if (!cfg$save_matched_rows_csv) {
      tmp$before_rows <- NULL
      tmp$matched_rows <- NULL
    }
    saveRDS(tmp, file.path(od, paste0(tag, "_full_result.rds")))
  }
  
  if (cfg$save_matched_rows_csv) {
    fwrite(res$before_rows, file.path(od, paste0(tag, "_before_rows.csv")))
    fwrite(res$matched_rows, file.path(od, paste0(tag, "_matched_rows.csv")))
  }
  
  fwrite(res$att, file.path(od, paste0(tag, "_att_summary.csv")))
  fwrite(res$delta, file.path(od, paste0(tag, "_delta_table.csv")))
  fwrite(res$block_att, file.path(od, paste0(tag, "_block_att.csv")))
  fwrite(res$contrib, file.path(od, paste0(tag, "_obs_contrib.csv")))
  fwrite(res$smd_level, file.path(od, paste0(tag, "_smd_pooled_hajek_level.csv")))
  fwrite(res$smd_variable, file.path(od, paste0(tag, "_smd_pooled_hajek_variable.csv")))
  fwrite(res$smd_plot, file.path(od, paste0(tag, "_smd_plot_summary.csv")))
  fwrite(res$smd_panel, file.path(od, paste0(tag, "_smd_master_panel.csv")))
  fwrite(res$drop, file.path(od, paste0(tag, "_drop_log.csv")))
}

all_att <- list()
all_delta <- list()
all_block <- list()
all_hajek <- list()
all_contrib <- list()
all_smd_level <- list()
all_smd_var <- list()
all_smd_plot <- list()
all_panel <- list()
all_drop <- list()

for (q in cfg$q_values) {
  DT <- prep_data(load_q_data(q, cfg$root_in))
  
  for (rv in cfg$region_vars) {
    if (!(rv %in% names(DT))) next
    
    for (d in cfg$day_windows) {
      res <- run_psm_one(DT, q, rv, d)
      
      if (is.null(res)) next
      
      save_one(res)
      
      key <- sprintf("%s_q%d_d%d", rv, q, d)
      
      all_att[[key]] <- res$att
      all_delta[[key]] <- res$delta
      all_block[[key]] <- res$block_att
      all_hajek[[key]] <- res$hajek
      all_contrib[[key]] <- res$contrib
      all_smd_level[[key]] <- res$smd_level
      all_smd_var[[key]] <- res$smd_variable
      all_smd_plot[[key]] <- res$smd_plot
      all_panel[[key]] <- res$smd_panel
      
      if (nrow(res$drop)) {
        dd <- copy(res$drop)
        dd[, `:=`(q = q, region_var = rv, temp_var = "temp_round", day_window = d)]
        all_drop[[key]] <- dd
      }
    }
  }
}

if (length(all_att)) {
  fwrite(rbindlist(all_att, fill = TRUE), file.path(cfg$root_out, "psm_all_att_summary.csv"))
}

if (length(all_delta)) {
  fwrite(rbindlist(all_delta, fill = TRUE), file.path(cfg$root_out, "psm_all_delta_table.csv"))
}

if (length(all_block)) {
  fwrite(rbindlist(all_block, fill = TRUE), file.path(cfg$root_out, "psm_all_block_att_summary.csv"))
}

if (length(all_hajek)) {
  fwrite(rbindlist(all_hajek, fill = TRUE), file.path(cfg$root_out, "psm_all_hajek_att_summary.csv"))
}

if (length(all_contrib)) {
  fwrite(rbindlist(all_contrib, fill = TRUE), file.path(cfg$root_out, "psm_all_obs_contrib.csv"))
}

if (length(all_smd_level)) {
  fwrite(rbindlist(all_smd_level, fill = TRUE), file.path(cfg$root_out, "psm_all_smd_pooled_hajek_level_summary.csv"))
}

if (length(all_smd_var)) {
  smdv <- rbindlist(all_smd_var, fill = TRUE)
  fwrite(smdv, file.path(cfg$root_out, "psm_all_smd_pooled_hajek_variable_summary.csv"))
  fwrite(smdv, file.path(cfg$root_out, "psm_all_smd_variable_summary.csv"))
}

if (length(all_smd_plot)) {
  fwrite(rbindlist(all_smd_plot, fill = TRUE), file.path(cfg$root_out, "psm_all_smd_plot_summary.csv"))
}

if (length(all_panel)) {
  fwrite(rbindlist(all_panel, fill = TRUE), file.path(cfg$root_out, "PSM_SMD_values_master_rt_weighted_mean.csv"))
}

if (length(all_drop)) {
  fwrite(rbindlist(all_drop, fill = TRUE), file.path(cfg$root_out, "psm_all_drop_log.csv"))
}

message("DONE")
