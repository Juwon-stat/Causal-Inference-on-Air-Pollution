###############################################################################
# Title   : Block-wise CBPS with unit-level Multiplier Bootstrap (ATT) + Diagnostics
# Author  : Juwon Jung
# Version : 1.0
# Date    : 2025-08-19
#
# Purpose :
#   - For each treatment date (block), fit CBPS(ATT) within a symmetric ±k-day
#     window and the same spatial units, then aggregate across blocks.
#   - Compute ATT for PM10 and PM2.5 using a global Hájek-style estimator.
#   - Perform unit-level multiplier bootstrap (Exp/Gamma with mean=1) to obtain
#     percentile 95% CI and bootstrap standard errors (se10, se25).
#   - Produce covariate balance diagnostics (SMD before/after weighting),
#     weight stability metrics (ESS, max_w_rel, etc.), and influence tables.
#
# Inputs  :
#   - Three RData files located in `root_in`, named:
#       seoul_q3.RData, seoul_q4.RData, seoul_q5.RData
#     Each file must contain an object named `seoul_q3` (or `_q4`, `_q5`) with at least:
#       date (Date or convertible), year (int), month, day, hour,
#       treated (0/1 for the calendar date in 2020),
#       PM10, PM25,
#       temp, temp_quantile,
#       rain_binary, wind_sp_binary, wind_dir_8, time_of_day,
#       CO_quantile, O3_quantile,
#       area (5-zone) and dist (25-district) identifiers,
#       wday (weekday label/numeric).
#
# Outputs :
#   - For each (q, unit_tag, day_window, region_var) combination:
#       * RData: "cbps_UNIT_ATT_boot_{q}_{unit_tag}_d{d}_{region_var}_WITH_SMD_DIAG.RData"
#                containing TALL, CALL, att_over, att10_B, att25_B, smd_table,
#                stab_table, drop_log, influence_table, influence_quant_dt,
#                se10, se25, n_boot_10, n_boot_25, B, region_var
#       * CSVs : SMD_*.csv, STAB_*.csv, DROP_*.csv, INFLUENCE_*.csv, INFLU_Q_*.csv
#   - A summary CSV over all combos:
#       "cbps_UNIT_ATT_summary_{region_var}.csv"
#
# Notes   :
#   - Control weights are scaled per block so that sum_w(control) = n_treated(block).
#   - unit-level multiplier bootstrap draws one multiplier per unique rid and reweights
#     both treated totals and control weighted totals consistently.
#   - Set `region_var <- "area"` or `"dist"` to switch the spatial matching unit.
#   - Set `unit_tag <- "tq"` for temp quantile matching or `"rt"` for rounded temp.
#
# Reproducibility :
#   - Use `future.seed=TRUE` and combine with Exp/Gamma multipliers (mean=1).
###############################################################################

suppressPackageStartupMessages({
  library(CBPS)
  library(future)
  library(future.apply)
  library(data.table)
})

## ============================== USER PATHS =================================
# Replace the following directories with your own.
root_in  <- "C:/path/to/Seoul_RData"                     # contains seoul_q{3,4,5}.RData
root_out <- "C:/path/to/Matching_Results/CBPS"           # outputs (RData/CSV)
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

## ============================= GLOBAL SETTINGS =============================
# Treatment period (inclusive). Adjust if needed.
start_date <- as.Date("2020-01-24")
end_date   <- as.Date("2020-02-09")
date_seq   <- as.Date(seq(start_date, end_date, by = "day"))

# Toggle spatial unit HERE: "area" (5 residential zones) or "dist" (25 districts)
region_var <- "dist"

# Temperature matching mode: "tq" (temp_quantile) or "rt" (rounded temp)
unit_tags   <- c("tq", "rt")

# Symmetric temporal windows in days for block-wise matching
day_windows <- c(3, 7, 14)

# Quantile-binning variants (input datasets)
q_list      <- c("q3", "q4", "q5")

# Parallelization toggle and bootstrap multiplier distribution
USE_PARALLEL <- TRUE        # TRUE = multisession; FALSE = sequential
MULT_ALPHA   <- 1           # 1: Exp(1); >1: Gamma(shape=α, rate=α) with mean 1

## ============================== UTILITIES ==================================
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  suppressWarnings(min(x, na.rm=TRUE))
}
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  suppressWarnings(max(x, na.rm=TRUE))
}

# Standardized Mean Difference for categorical covariates (before vs after weighting).
# pool_df must include: treated (0/1), covars as factors, and weight column w (controls only).
compute_smd_cat <- function(pool_df, covars, w_col="w"){
  out_list <- vector("list", length(covars))
  Trows <- pool_df$treated == 1L
  Crows <- !Trows
  nT <- sum(Trows, na.rm=TRUE); nC <- sum(Crows, na.rm=TRUE)
  wC_sum <- if (w_col %in% names(pool_df)) sum(as.double(pool_df[[w_col]][Crows]), na.rm=TRUE) else NA_real_
  
  for (j in seq_along(covars)){
    v <- covars[j]
    if (!v %in% names(pool_df)) next
    levs <- sort(unique(as.character(pool_df[[v]])))
    facT <- factor(as.character(pool_df[[v]][Trows]), levels = levs)
    facC <- factor(as.character(pool_df[[v]][Crows]), levels = levs)
    
    Nt <- as.numeric(table(facT)); Nc <- as.numeric(table(facC))
    pT_b <- if (nT>0) Nt/nT else rep(NA_real_, length(levs))
    pC_b <- if (nC>0) Nc/nC else rep(NA_real_, length(levs))
    den_b <- sqrt((pT_b*(1-pT_b) + pC_b*(1-pC_b))/2)
    smd_b <- (pC_b - pT_b)/den_b
    
    if (!is.na(wC_sum) && is.finite(wC_sum) && wC_sum>0){
      w_by_lev <- as.numeric(xtabs(pool_df[[w_col]][Crows] ~ facC, drop = FALSE))
      pC_a <- w_by_lev / wC_sum
    } else {
      pC_a <- rep(NA_real_, length(levs))
    }
    den_a <- sqrt((pT_b*(1-pT_b) + pC_a*(1-pC_a))/2)
    smd_a <- (pC_a - pT_b)/den_a
    
    smd_b[!is.finite(smd_b)] <- NA_real_
    smd_a[!is.finite(smd_a)] <- NA_real_
    
    out_list[[j]] <- data.table(
      variable   = v,
      level      = levs,
      pT_before  = pT_b,
      pC_before  = pC_b,
      pC_after   = pC_a,
      smd_before = smd_b,
      smd_after  = smd_a,
      abs_before = abs(smd_b),
      abs_after  = abs(smd_a),
      nT         = nT,
      nC         = nC,
      wC         = wC_sum
    )
  }
  rbindlist(out_list, use.names=TRUE, fill=TRUE)
}

# Weight stability metrics (controls only) BEFORE block-wise scaling.
compute_stability <- function(pool_df, w_col="w"){
  Crows <- pool_df$treated==0L
  if (!any(Crows)) {
    return(data.table(
      ESS_C=NA_real_, max_w=NA_real_, max_w_rel=NA_real_,
      sum_w=NA_real_, n_control=0L, n_treated=sum(pool_df$treated==1L, na.rm=TRUE)
    ))
  }
  w <- as.double(pool_df[[w_col]][Crows])
  sumw <- sum(w, na.rm=TRUE); sumw2 <- sum(w*w, na.rm=TRUE)
  ess <- if (sumw2>0) (sumw*sumw)/sumw2 else NA_real_
  data.table(
    ESS_C=ess,
    max_w=safe_max(w),
    max_w_rel=if (sumw>0) safe_max(w)/sumw else NA_real_,
    sum_w=sumw,
    n_control=sum(Crows, na.rm=TRUE),
    n_treated=sum(pool_df$treated==1L, na.rm=TRUE)
  )
}

## ============================ EXPERIMENT GRID ==============================
combos <- data.table::CJ(q=q_list, unit_tag=unit_tags, day_window=day_windows)

## ============================== CORE ROUTINE ===============================
run_combo <- function(idx){
  q        <- combos$q[idx]
  unit_tag <- combos$unit_tag[idx]
  day_w    <- combos$day_window[idx]
  
  # Tag used in all filenames: includes region_var for clarity.
  fname_tag <- sprintf("%s_%s_d%d_%s", q, unit_tag, day_w, region_var)
  
  # ---- Load data for this quantile scheme ----
  load(file.path(root_in, sprintf("seoul_%s.RData", q)))
  seoul <- get(sprintf("seoul_%s", q))
  seoul <- as.data.frame(seoul, stringsAsFactors = FALSE)
  seoul$date <- as.Date(seoul$date)
  
  # Unique observation unit (rid)
  if (!("rid" %in% names(seoul))) seoul$rid <- seq_len(nrow(seoul))
  
  # Basic numeric coercion (defensive)
  num_cols0 <- intersect(c("PM10","PM25","temp","wind_sp_binary","hour","day","month","year"),
                         names(seoul))
  for (nm in num_cols0) seoul[[nm]] <- to_num(seoul[[nm]])
  
  # Temperature matching variable
  if (unit_tag=="tq") {
    seoul$temp_match <- as.character(seoul$temp_quantile)
  } else {
    seoul$temp_match <- as.character(round(seoul$temp))  # e.g., round(temp/2)*2 for coarser bins
  }
  
  # Covariates used in CBPS (all treated as factors)
  covars <- c("wday","temp_match","rain_binary","CO_quantile","O3_quantile",
              "wind_sp_binary","wind_dir_8","time_of_day", region_var)
  
  ps_form <- as.formula(
    paste("treated ~", paste(c("wday","temp_match","rain_binary",
                               "CO_quantile","O3_quantile",
                               "wind_sp_binary","wind_dir_8",
                               "time_of_day", region_var), collapse=" + "))
  )
  
  # Containers for global (across blocks) aggregation
  TALL_list <- list()   # treated totals per rid per block (then summed)
  CALL_list <- list()   # control weighted totals per rid per block (then summed)
  
  # Diagnostics
  smd_list      <- list()
  stab_list     <- list()
  drop_log_list <- list()
  infl_raw_list <- list()  # per (rid, date) scaled weight summaries
  
  # =========================== BLOCK LOOP (by date) =========================
  for (td in date_seq){
    td  <- as.Date(td)
    win <- as.Date(td) + (-day_w:day_w)
    
    # Treated pool = 2020 observations on the target date td
    trt <- seoul[seoul$date==td & seoul$year==2020, , drop=FALSE]
    if (!nrow(trt)) { drop_log_list[[length(drop_log_list)+1L]] <- data.table(date=td, reason="no_treated"); next }
    trt$treated <- 1L
    
    # Controls = same calendar month-day across 2017–2019 & 2021–2023 within ±k days
    md <- format(win, "%m-%d")
    ctrl <- seoul[format(seoul$date,"%m-%d") %in% md & seoul$year!=2020, , drop=FALSE]
    # Restrict controls to the same spatial units present among treated
    ctrl <- ctrl[ ctrl[[region_var]] %in% unique(trt[[region_var]]), , drop=FALSE]
    if (!nrow(ctrl)) { drop_log_list[[length(drop_log_list)+1L]] <- data.table(date=td, reason="no_controls"); next }
    ctrl$treated <- 0L
    
    pool <- rbind(trt, ctrl)
    
    # Drop rows with missing covariates
    cc_idx <- stats::complete.cases(pool[, covars])
    if (!all(cc_idx)) {
      drop_log_list[[length(drop_log_list)+1L]] <- data.table(date=td, reason="covariate_missing", n_drop=sum(!cc_idx))
    }
    pool <- pool[cc_idx, , drop=FALSE]
    if (!nrow(pool)) { drop_log_list[[length(drop_log_list)+1L]] <- data.table(date=td, reason="all_dropped_cc"); next }
    
    # Factorize covariates; ensure outcomes numeric; treated as integer
    for (nm in covars) pool[[nm]] <- as.factor(pool[[nm]])
    pool$PM10 <- to_num(pool$PM10); pool$PM25 <- to_num(pool$PM25)
    pool$treated <- as.integer(pool$treated)
    
    # --------------------------- CBPS (ATT) fit -----------------------------
    fit <- tryCatch(CBPS(ps_form, data=pool, ATT=TRUE), error=function(e) NULL)
    if (is.null(fit) || anyNA(fit$weights)) {
      drop_log_list[[length(drop_log_list)+1L]] <- data.table(date=td, reason="CBPS_fail_or_NAweights")
      next
    }
    pool$w <- as.double(fit$weights)
    
    # ------------------ SMD & stability BEFORE scaling ----------------------
    smd_dt  <- compute_smd_cat(pool, covars, w_col="w"); smd_dt[, date := td]
    stab_dt <- compute_stability(pool, w_col="w");      stab_dt[, date := td]
    
    # ----------- Scale control weights so sum_w(control) = n_treated --------
    Tr <- pool$treated==1L; Cr <- !Tr
    nT  <- sum(Tr, na.rm=TRUE)
    wC  <- sum(as.double(pool$w[Cr]), na.rm=TRUE)
    if (!is.finite(nT) || nT<=0 || !is.finite(wC) || wC<=0) {
      drop_log_list[[length(drop_log_list)+1L]] <- data.table(date=td, reason="nonpositive_nT_or_wC", nT=nT, wC=wC)
      next
    }
    s_day <- nT / wC                         # per-block scaling factor
    w_scaled_vec <- s_day * as.double(pool$w[Cr])
    
    # ---------- Stability AFTER scaling (reference diagnostics) -------------
    sumw_after  <- sum(w_scaled_vec, na.rm=TRUE)
    sumw2_after <- sum(w_scaled_vec*w_scaled_vec, na.rm=TRUE)
    ess_after   <- if (sumw2_after>0) (sumw_after^2)/sumw2_after else NA_real_
    max_w_after <- safe_max(w_scaled_vec)
    max_w_rel_after <- if (sumw_after>0) max_w_after/sumw_after else NA_real_
    stab_dt[, `:=`(
      ESS_C_after      = ess_after,
      max_w_after      = max_w_after,
      max_w_rel_after  = max_w_rel_after,
      sum_w_after      = sumw_after
    )]
    
    smd_list[[length(smd_list)+1L]]   <- smd_dt
    stab_list[[length(stab_list)+1L]] <- stab_dt
    
    # ----------- Immediate per-block rid aggregation (memory-safe) ----------
    if (any(Tr)){
      dt_t <- data.table(rid = pool$rid[Tr],
                         Y10 = as.double(pool$PM10[Tr]),
                         Y25 = as.double(pool$PM25[Tr]))
      dt_t <- dt_t[, .(Y10 = sum(Y10), Y25 = sum(Y25), N = .N), by = .(rid)]
      TALL_list[[length(TALL_list)+1L]] <- dt_t
    }
    if (any(Cr)){
      dt_c <- data.table(rid  = pool$rid[Cr],
                         wY10 = w_scaled_vec * as.double(pool$PM10[Cr]),
                         wY25 = w_scaled_vec * as.double(pool$PM25[Cr]),
                         w    = w_scaled_vec)
      dt_c <- dt_c[, .(wY10 = sum(wY10), wY25 = sum(wY25), w = sum(w)), by = .(rid)]
      CALL_list[[length(CALL_list)+1L]] <- dt_c
      
      # Influence raw: per (rid, date) scaled w
      infl_raw_list[[length(infl_raw_list)+1L]] <- dt_c[, .(rid, date = td, w = w)]
    }
  } # end date loop
  
  # If no valid blocks, save drop log and return sentinel row
  if (!length(TALL_list) || !length(CALL_list)) {
    drop_log <- if (length(drop_log_list)) rbindlist(drop_log_list, fill=TRUE) else data.table()
    if (nrow(drop_log)) fwrite(drop_log, file.path(root_out, sprintf("DROP_%s.csv", fname_tag)))
    return(list(
      q=q, unit_tag=unit_tag, window=day_w, region_key=region_var,
      PM10=NA_real_, PM25=NA_real_,
      lo10=NA_real_, hi10=NA_real_, lo25=NA_real_, hi25=NA_real_,
      se10=NA_real_, se25=NA_real_, n_boot_10=NA_integer_, n_boot_25=NA_integer_,
      error=TRUE, message="no valid blocks"
    ))
  }
  
  # ====================== GLOBAL AGGREGATION (over blocks) ==================
  TALL <- rbindlist(TALL_list)[, .(Y10 = sum(Y10), Y25 = sum(Y25), N = sum(N)), by = .(rid)]
  CALL <- rbindlist(CALL_list)[, .(wY10 = sum(wY10), wY25 = sum(wY25), w = sum(w)), by = .(rid)]
  
  # Global Hájek-style point estimates (control weighted mean - treated mean)
  point10 <- (sum(CALL$wY10, na.rm=TRUE)/sum(CALL$w, na.rm=TRUE)) -
    (sum(TALL$Y10,  na.rm=TRUE)/sum(TALL$N, na.rm=TRUE))
  point25 <- (sum(CALL$wY25, na.rm=TRUE)/sum(CALL$w, na.rm=TRUE)) -
    (sum(TALL$Y25,  na.rm=TRUE)/sum(TALL$N, na.rm=TRUE))
  
  # =============================== INFLUENCE ================================
  influence_table <- NULL
  influence_quant_dt <- NULL
  if (length(infl_raw_list)){
    infl_raw <- rbindlist(infl_raw_list)
    setDT(infl_raw)
    influence_table <- infl_raw[, .(w_sum = sum(w),
                                    days  = uniqueN(date)), by = .(rid)]
    qs <- c(.5, .9, .95, .99)
    influence_quant_dt <- data.table(
      prob = qs,
      quantile = as.numeric(quantile(influence_table$w_sum, probs = qs, na.rm=TRUE))
    )
  }
  
  # ============================ BOOTSTRAP (unit) =============================
  B <- 1000L
  rids <- sort(unique(c(TALL$rid, CALL$rid)))
  key_t <- match(TALL$rid, rids)
  key_c <- match(CALL$rid, rids)
  
  att10_B <- numeric(B)
  att25_B <- numeric(B)
  
  draw_mult <- function(n){
    if (MULT_ALPHA == 1) rexp(n, rate=1) else rgamma(n, shape=MULT_ALPHA, rate=MULT_ALPHA) # mean=1
  }
  
  for (b in seq_len(B)){
    m_id <- draw_mult(length(rids))
    
    # Treated totals
    m_t   <- m_id[key_t]
    Y10Tb <- sum(m_t * TALL$Y10)
    Y25Tb <- sum(m_t * TALL$Y25)
    NTb   <- sum(m_t * TALL$N)
    
    # Control weighted totals
    m_c    <- m_id[key_c]
    wY10Cb <- sum(m_c * CALL$wY10)
    wY25Cb <- sum(m_c * CALL$wY25)
    wCb    <- sum(m_c * CALL$w)
    
    if (is.finite(NTb) && NTb>0 && is.finite(wCb) && wCb>0){
      att10_B[b] <- (wY10Cb / wCb) - (Y10Tb / NTb)
      att25_B[b] <- (wY25Cb / wCb) - (Y25Tb / NTb)
    } else {
      att10_B[b] <- NA_real_; att25_B[b] <- NA_real_
    }
  }
  
  ci10 <- stats::quantile(att10_B, c(.025,.975), na.rm=TRUE, names=FALSE)
  ci25 <- stats::quantile(att25_B, c(.025,.975), na.rm=TRUE, names=FALSE)
  se10 <- stats::sd(att10_B, na.rm=TRUE)
  se25 <- stats::sd(att25_B, na.rm=TRUE)
  n_boot_10 <- sum(is.finite(att10_B))
  n_boot_25 <- sum(is.finite(att25_B))
  
  # ============================= OUTPUT BUNDLES =============================
  att_over <- data.frame(
    pollutant=c("PM10","PM25"),
    point=c(point10, point25),
    lo95=c(ci10[1], ci25[1]),
    hi95=c(ci10[2], ci25[2]),
    se=c(se10, se25)
  )
  
  smd_table  <- if (length(smd_list)) rbindlist(smd_list, use.names=TRUE, fill=TRUE) else data.table()
  stab_table <- if (length(stab_list)) rbindlist(stab_list, use.names=TRUE, fill=TRUE) else data.table()
  if (nrow(stab_table)) setorder(stab_table, date)
  drop_log   <- if (length(drop_log_list)) rbindlist(drop_log_list, fill=TRUE) else data.table()
  
  # Attach meta tags
  if (nrow(smd_table))  smd_table[,  `:=`(q = q, unit_tag = unit_tag, window = day_w, region_var = region_var)]
  if (nrow(stab_table)) stab_table[, `:=`(q = q, unit_tag = unit_tag, window = day_w, region_var = region_var)]
  if (nrow(drop_log))   drop_log[,   `:=`(q = q, unit_tag = unit_tag, window = day_w, region_var = region_var)]
  if (!is.null(influence_table) && nrow(influence_table))
    influence_table[, `:=`(q = q, unit_tag = unit_tag, window = day_w, region_var = region_var)]
  if (!is.null(influence_quant_dt) && nrow(influence_quant_dt))
    influence_quant_dt[, `:=`(q = q, unit_tag = unit_tag, window = day_w, region_var = region_var)]
  
  # Save a comprehensive RData bundle
  save(TALL, CALL, att_over, att10_B, att25_B,
       smd_table, stab_table, drop_log,
       influence_table, influence_quant_dt,
       se10, se25, n_boot_10, n_boot_25, B, region_var,
       file=file.path(root_out, sprintf("cbps_UNIT_ATT_boot_%s_WITH_SMD_DIAG.RData", fname_tag)))
  
  # Optional: write diagnostics as CSVs for quick inspection
  if (nrow(smd_table))         fwrite(smd_table,         file.path(root_out, sprintf("SMD_%s.csv", fname_tag)))
  if (nrow(stab_table))        fwrite(stab_table,        file.path(root_out, sprintf("STAB_%s.csv", fname_tag)))
  if (nrow(drop_log))          fwrite(drop_log,          file.path(root_out, sprintf("DROP_%s.csv", fname_tag)))
  if (!is.null(influence_table) && nrow(influence_table))
    fwrite(influence_table,    file.path(root_out, sprintf("INFLUENCE_%s.csv", fname_tag)))
  if (!is.null(influence_quant_dt) && nrow(influence_quant_dt))
    fwrite(influence_quant_dt, file.path(root_out, sprintf("INFLU_Q_%s.csv", fname_tag)))
  
  list(
    q=q, unit_tag=unit_tag, window=day_w, region_key=region_var,
    PM10=point10, PM25=point25,
    lo10=ci10[1], hi10=ci10[2],
    lo25=ci25[1], hi25=ci25[2],
    se10=se10, se25=se25, n_boot_10=n_boot_10, n_boot_25=n_boot_25,
    ESS_mean       = if (nrow(stab_table)) mean(stab_table$ESS_C, na.rm=TRUE) else NA_real_,
    ESS_median     = if (nrow(stab_table)) stats::median(stab_table$ESS_C, na.rm=TRUE) else NA_real_,
    ESS_min        = if (nrow(stab_table)) safe_min(stab_table$ESS_C) else NA_real_,
    max_w_rel_max  = if (nrow(stab_table)) safe_max(stab_table$max_w_rel) else NA_real_,
    ESS_after_mean = if (nrow(stab_table)) mean(stab_table$ESS_C_after, na.rm=TRUE) else NA_real_,
    error=FALSE, message=NA_character_
  )
}

# Safe wrapper to prevent the whole job from failing on one combo
safe_run <- function(i){
  tryCatch(
    run_combo(i),
    error = function(e) list(
      q = combos$q[i], unit_tag = combos$unit_tag[i], window = combos$day_window[i], region_key=region_var,
      PM10=NA_real_, PM25=NA_real_, lo10=NA_real_, hi10=NA_real_, lo25=NA_real_, hi25=NA_real_,
      se10=NA_real_, se25=NA_real_, n_boot_10=NA_integer_, n_boot_25=NA_integer_,
      ESS_mean=NA_real_, ESS_median=NA_real_, ESS_min=NA_real_, max_w_rel_max=NA_real_,
      ESS_after_mean=NA_real_,
      error = TRUE, message = conditionMessage(e)
    )
  )
}

## =============================== EXECUTION ================================
do_all <- function() {
  if (USE_PARALLEL) {
    workers <- max(1, parallel::detectCores() - 1)
    plan(multisession, workers = workers)
    on.exit(plan(sequential), add = TRUE)
    withCallingHandlers(
      future_lapply(seq_len(nrow(combos)), safe_run, future.seed = TRUE),
      warning = function(w){
        # Suppress benign "Future was canceled" warnings in some parallel backends
        if (grepl("Future .* was canceled", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  } else {
    lapply(seq_len(nrow(combos)), safe_run)
  }
}

att_collect <- do_all()
att_table <- data.table::rbindlist(att_collect, fill=TRUE)

# Final summary across all combos; filename includes region_var
out_csv <- file.path(root_out, sprintf("cbps_UNIT_ATT_summary_%s.csv", region_var))
data.table::fwrite(att_table, out_csv)

print(att_table)
