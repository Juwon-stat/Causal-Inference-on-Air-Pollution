###############################################################################
# Title   : CEM_GLOBAL.R
# Purpose : Block-wise Coarsened Exact Matching (CEM)
#           + treated-level ATT
#           + conventional bootstrap
#           + SMD diagnostics
#
# Design  : Matches Algorithm 1 in the manuscript:
#           for q in {3,4,5}, d in {3,7,14},
#           for each treated observation i,
#           find controls j with |t_j - t_i| <= d and year(j) != 2020,
#           then exact match on coarsened covariates X_i^(q) = X_j^(q).
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# 0. USER CONFIGURATION
# =============================================================================

cfg <- list(
  # Folder containing seoul_q3.RData, seoul_q4.RData, seoul_q5.RData
  root_in  = "C:/path/to/Seoul_RData",

  # Output folder
  root_out = "C:/path/to/CEM_results",

  # O3 and CO quantile specifications
  q_values = c(3, 4, 5),

  # Time windows in days
  day_windows = c(3, 7, 14),

  # Spatial specifications
  # Main analysis: "dist"
  # Residential-zone sensitivity: "area"
  region_vars = c("dist", "area"),

  # Treatment period
  treatment_start = as.Date("2020-01-24"),
  treatment_end   = as.Date("2020-02-09"),
  treatment_year  = 2020L,

  # Bootstrap
  B = 1000L,
  seed = 20260507L,

  # Outcomes
  outcomes = c("PM10", "PM25"),

  # Save full matched object? Recommend: FALSE
  save_full_rds = FALSE
)

dir.create(cfg$root_out, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. HELPER FUNCTIONS
# =============================================================================

stop_if_missing <- function(DT, vars, data_name = "data") {
  missing_vars <- setdiff(vars, names(DT))

  if (length(missing_vars) > 0) {
    stop(
      sprintf(
        "[%s] Missing required variables: %s",
        data_name,
        paste(missing_vars, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

safe_as_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  as.Date(x)
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
    } else {
      return(NA_real_)
    }
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
    if (!(v %in% names(treated_dt)) || !(v %in% names(control_dt))) next

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

conventional_bootstrap_att <- function(delta_dt, outcomes, B, seed) {
  set.seed(seed)

  n <- nrow(delta_dt)

  if (n == 0) {
    stop("No matched treated observations available for bootstrap.", call. = FALSE)
  }

  boot_out <- vector("list", length(outcomes))
  names(boot_out) <- outcomes

  for (outcome in outcomes) {
    delta_col <- paste0(outcome, "_diff")
    delta <- delta_dt[[delta_col]]

    if (all(is.na(delta))) {
      boot_out[[outcome]] <- data.table(
        outcome = outcome,
        att = NA_real_,
        boot_se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        n_boot = B
      )
      next
    }

    boot_stat <- replicate(
      B,
      mean(delta[sample.int(n, size = n, replace = TRUE)], na.rm = TRUE)
    )

    boot_out[[outcome]] <- data.table(
      outcome = outcome,
      att = mean(delta, na.rm = TRUE),
      boot_se = sd(boot_stat, na.rm = TRUE),
      ci_low = as.numeric(quantile(boot_stat, 0.025, na.rm = TRUE)),
      ci_high = as.numeric(quantile(boot_stat, 0.975, na.rm = TRUE)),
      n_boot = B
    )
  }

  rbindlist(boot_out, use.names = TRUE, fill = TRUE)
}

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

  # Variables such as year, wday, time_of_day, rain_binary, wind_sp_ind,
  # wind_dir_8, CO_quantile, and O3_quantile are assumed to already exist.

  # Date conversion only
  DT[, date := safe_as_date(date)]

  # Stable observation ID
  setorder(DT, date)
  DT[, obs_id := .I]

  # Month-day key for cross-year calendar-window matching
  DT[, month_day := format(date, "%m-%d")]

  # Temperature is always matched using rounded temperature.
  # This is the only derived covariate created inside the CEM code.
  DT[, temp_round := round(temp)]

  DT
}

get_match_vars <- function(region_var) {
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

make_treated_window_table <- function(treated_key, d) {
  by_cols <- names(treated_key)

  treated_key[
    ,
    .(month_day = make_month_day_window(treated_date[1], d)),
    by = by_cols
  ]
}

# =============================================================================
# 2. ONE CEM RUN
# =============================================================================

run_cem_one <- function(DT, q, region_var, d, cfg) {
  message(
    sprintf(
      "\n[CEM] q=%d | region=%s | temp=round(temp) | window=±%d days",
      q, region_var, d
    )
  )

  match_vars <- get_match_vars(region_var)

  required_vars <- unique(c(
    "obs_id",
    "date",
    "year",
    "month_day",
    cfg$outcomes,
    match_vars
  ))

  stop_if_missing(DT, required_vars, data_name = sprintf("q%d data", q))

  work <- copy(DT)

  for (v in match_vars) {
    work[, (v) := as.character(get(v))]
  }

  # ---------------------------------------------------------------------------
  # 2.1 Define treated and control pools
  # ---------------------------------------------------------------------------

  treated <- work[
    date >= cfg$treatment_start &
      date <= cfg$treatment_end &
      year == cfg$treatment_year
  ]

  controls <- work[
    year != cfg$treatment_year
  ]

  if (nrow(treated) == 0) {
    stop("No treated observations found.", call. = FALSE)
  }

  if (nrow(controls) == 0) {
    stop("No control observations found.", call. = FALSE)
  }

  treated[, target_index := .I]

  treated_key_cols <- c(
    "target_index",
    "obs_id",
    "date",
    match_vars,
    cfg$outcomes
  )

  treated_key <- treated[, ..treated_key_cols]

  setnames(treated_key, "obs_id", "treated_obs_id")
  setnames(treated_key, "date", "treated_date")
  setnames(
    treated_key,
    old = cfg$outcomes,
    new = paste0(cfg$outcomes, "_treated")
  )

  treated_window <- make_treated_window_table(treated_key, d)

  control_cols <- unique(c(
    "obs_id",
    "date",
    "year",
    "month_day",
    match_vars,
    cfg$outcomes
  ))

  controls_key <- controls[, ..control_cols]

  setnames(controls_key, "obs_id", "control_obs_id")
  setnames(controls_key, "date", "control_date")

  # ---------------------------------------------------------------------------
  # 2.2 CEM exact matching
  # ---------------------------------------------------------------------------

  join_vars <- c("month_day", match_vars)

  matched_controls <- controls_key[
    treated_window,
    on = join_vars,
    allow.cartesian = TRUE,
    nomatch = 0
  ]

  if (nrow(matched_controls) == 0) {
    warning(
      sprintf(
        "No matches: q=%d, region=%s, d=%d",
        q, region_var, d
      )
    )
    return(NULL)
  }

  matched_controls <- unique(
    matched_controls,
    by = c("target_index", "control_obs_id")
  )

  n_control_by_target <- matched_controls[
    ,
    .(
      n_controls = .N,
      n_unique_controls = uniqueN(control_obs_id)
    ),
    by = target_index
  ]

  matched_controls <- merge(
    matched_controls,
    n_control_by_target,
    by = "target_index",
    all.x = TRUE
  )

  # ---------------------------------------------------------------------------
  # 2.3 Treated-level ATT contribution
  # ---------------------------------------------------------------------------

  control_means <- matched_controls[
    ,
    lapply(.SD, mean, na.rm = TRUE),
    by = target_index,
    .SDcols = cfg$outcomes
  ]

  setnames(
    control_means,
    old = cfg$outcomes,
    new = paste0(cfg$outcomes, "_control_mean")
  )

  treated_cols <- c(
    "target_index",
    "treated_obs_id",
    "treated_date",
    match_vars,
    paste0(cfg$outcomes, "_treated"),
    "n_controls",
    "n_unique_controls"
  )

  treated_matched <- unique(
    matched_controls[, ..treated_cols],
    by = "target_index"
  )

  treated_contrib <- merge(
    treated_matched,
    control_means,
    by = "target_index",
    all.x = TRUE
  )

  for (outcome in cfg$outcomes) {
    treated_col <- paste0(outcome, "_treated")
    control_col <- paste0(outcome, "_control_mean")
    diff_col <- paste0(outcome, "_diff")

    treated_contrib[
      ,
      (diff_col) := get(control_col) - get(treated_col)
    ]
  }

  # ---------------------------------------------------------------------------
  # 2.4 Conventional bootstrap over matched treated-level contrasts
  # ---------------------------------------------------------------------------

  boot_seed <- cfg$seed +
    q * 10000L +
    d * 100L +
    match(region_var, cfg$region_vars) * 10L

  boot_table <- conventional_bootstrap_att(
    delta_dt = treated_contrib,
    outcomes = cfg$outcomes,
    B = cfg$B,
    seed = boot_seed
  )

  boot_table[
    ,
    `:=`(
      method = "CEM",
      q = q,
      region_var = region_var,
      temp_var = "temp_round",
      day_window = d,
      n_treated_total = nrow(treated),
      n_treated_matched = nrow(treated_contrib),
      n_control_rows = nrow(matched_controls),
      n_control_unique = uniqueN(matched_controls$control_obs_id)
    )
  ]

  setcolorder(
    boot_table,
    c(
      "method",
      "q",
      "region_var",
      "temp_var",
      "day_window",
      "outcome",
      "att",
      "boot_se",
      "ci_low",
      "ci_high",
      "n_boot",
      "n_treated_total",
      "n_treated_matched",
      "n_control_rows",
      "n_control_unique"
    )
  )

  # ---------------------------------------------------------------------------
  # 2.5 SMD before and after matching
  # ---------------------------------------------------------------------------

  smd_covariates <- match_vars

  # BEFORE:
  # Compare treated observations against the union of controls eligible
  # under the same ±d calendar window, before exact covariate matching.
  # We only need unique eligible month-days to avoid a large cartesian join.
  eligible_month_days <- unique(treated_window$month_day)

  eligible_controls <- controls_key[
    month_day %in% eligible_month_days
  ]

  eligible_controls <- unique(
    eligible_controls,
    by = "control_obs_id"
  )

  before_treated <- treated[, c(smd_covariates), with = FALSE]
  before_treated[, smd_weight := 1]

  before_control <- eligible_controls[, c(smd_covariates), with = FALSE]
  before_control[, smd_weight := 1]

  smd_before_level <- compute_smd_categorical(
    treated_dt = before_treated,
    control_dt = before_control,
    covariates = smd_covariates,
    label = "before"
  )

  # AFTER:
  # ATT uses treated-level normalized matched-control means.
  # SMD uses matched-row-level diagnostic weights for the love plot.
  after_treated <- treated_contrib[, c(smd_covariates), with = FALSE]
  after_treated[, smd_weight := 1]

  after_control <- matched_controls[, c(smd_covariates), with = FALSE]
  after_control[, smd_weight := 1]

  smd_after_level <- compute_smd_categorical(
    treated_dt = after_treated,
    control_dt = after_control,
    covariates = smd_covariates,
    label = "after"
  )

  smd_level <- rbindlist(
    list(smd_before_level, smd_after_level),
    use.names = TRUE,
    fill = TRUE
  )

  smd_variable <- summarize_smd_by_variable(smd_level)

  if (nrow(smd_level) > 0) {
    smd_level[
      ,
      `:=`(
        method = "CEM",
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
        method = "CEM",
        q = q,
        region_var = region_var,
        temp_var = "temp_round",
        day_window = d
      )
    ]
  }

  list(
    config = cfg,
    q = q,
    region_var = region_var,
    temp_var = "temp_round",
    day_window = d,
    match_vars = match_vars,
    matched_controls = matched_controls,
    treated_contrib = treated_contrib,
    att_bootstrap = boot_table,
    smd_level = smd_level,
    smd_variable = smd_variable
  )
}

# =============================================================================
# 3. SAVE FUNCTION
# =============================================================================

save_cem_result <- function(res, cfg) {
  if (is.null(res)) return(invisible(NULL))

  out_dir <- file.path(
    cfg$root_out,
    sprintf("region_%s", res$region_var),
    sprintf("q%d", res$q)
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tag <- sprintf(
    "cem_%s_rt_q%d_d%d",
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
    light_res$matched_controls <- NULL
    saveRDS(
      light_res,
      file = file.path(out_dir, paste0(tag, "_result_no_matched_controls.rds"))
    )
  }

  fwrite(
    res$treated_contrib,
    file = file.path(out_dir, paste0(tag, "_treated_contrib.csv"))
  )

  fwrite(
    res$att_bootstrap,
    file = file.path(out_dir, paste0(tag, "_att_bootstrap.csv"))
  )

  fwrite(
    res$smd_level,
    file = file.path(out_dir, paste0(tag, "_smd_by_level.csv"))
  )

  fwrite(
    res$smd_variable,
    file = file.path(out_dir, paste0(tag, "_smd_by_variable.csv"))
  )

  message(sprintf("[Saved] %s", out_dir))
}

# =============================================================================
# 4. MAIN LOOP
# =============================================================================

all_att <- list()
all_smd_variable <- list()

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
      res <- run_cem_one(
        DT = DT,
        q = q,
        region_var = region_var,
        d = d,
        cfg = cfg
      )

      if (!is.null(res)) {
        save_cem_result(res, cfg)

        key <- sprintf("%s_rt_q%d_d%d", region_var, q, d)
        all_att[[key]] <- res$att_bootstrap
        all_smd_variable[[key]] <- res$smd_variable
      }
    }
  }
}

# =============================================================================
# 5. GLOBAL SUMMARY FILES
# =============================================================================

if (length(all_att) > 0) {
  att_all <- rbindlist(all_att, use.names = TRUE, fill = TRUE)

  fwrite(
    att_all,
    file.path(cfg$root_out, "cem_all_att_bootstrap_summary.csv")
  )

  message(
    sprintf(
      "[Saved] %s",
      file.path(cfg$root_out, "cem_all_att_bootstrap_summary.csv")
    )
  )
}

if (length(all_smd_variable) > 0) {
  smd_all <- rbindlist(all_smd_variable, use.names = TRUE, fill = TRUE)

  fwrite(
    smd_all,
    file.path(cfg$root_out, "cem_all_smd_variable_summary.csv")
  )

  smd_plot_all <- copy(smd_all)
  setnames(smd_plot_all, "sample", "stage")
  setnames(smd_plot_all, "max_abs_smd", "SMD")

  smd_plot_all <- smd_plot_all[
    ,
    .(
      method,
      q,
      region_var,
      temp_var,
      day_window,
      covariate,
      stage,
      SMD,
      mean_abs_smd,
      n_levels
    )
  ]

  fwrite(
    smd_plot_all,
    file.path(cfg$root_out, "cem_all_smd_plot_summary.csv")
  )

  message(
    sprintf(
      "[Saved] %s",
      file.path(cfg$root_out, "cem_all_smd_variable_summary.csv")
    )
  )

  message(
    sprintf(
      "[Saved] %s",
      file.path(cfg$root_out, "cem_all_smd_plot_summary.csv")
    )
  )
}

message("\nCEM matching + conventional bootstrap completed.")
