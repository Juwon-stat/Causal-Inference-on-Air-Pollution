###############################################################################
# Title   : Block-wise PSM (Propensity Score Matching) with unit-level
#           Multiplier Bootstrap (ATT, Control − Treat) + Diagnostics
# Author  : Juwon Jung
# Version : 1.0
# Date    : 2025-08-19
#
# Purpose :
#   - For each treatment date (block), run nearest-neighbor PSM within a
#     symmetric ±k-day window and identical spatial units.
#   - Compute block-level paired contributions d_i = (Control − Treat) for PM10/PM2.5,
#     then summarize across blocks (simple mean over matched sets).
#   - Use unit-level multiplier bootstrap (Exp/Gamma with mean=1) to obtain
#     percentile 95% CIs and bootstrap standard errors (se10, se25).
#   - Provide covariate balance diagnostics (categorical & continuous SMD),
#     post-match weight stability (ESS, max_w_rel), and influence snapshots.
#
# Inputs  :
#   - Three RData files located in `root_in`:
#       seoul_q3.RData, seoul_q4.RData, seoul_q5.RData
#     Each file must contain `seoul_q3` (or `_q4`, `_q5`) with at least:
#       date (Date), year, month, day, hour, treated (0/1 for 2020 target days),
#       PM10, PM25,
#       temp, temp_quantile,
#       rain_binary, wind_sp_binary, wind_dir_8, time_of_day,
#       CO_quantile, O3_quantile,
#       area (5 zones) and dist (25 districts),
#       wday (weekday, precomputed).
#
# Outputs :
#   - For each (q, unit_tag, day_window, region_var) combination:
#       * RData: "psm_UNIT_ATT_boot_{q}_{unit_tag}_{region_var}_d{d}_WITH_SMD_DIAG.RData"
#                containing DALL, att10_B, att25_B, smd_table, stab_table,
#                influence_table, top_control_reuse, top_pairs_d10/d25,
#                drop_log, and run settings (K_MAX, CALIPER_SD, ...),
#                along with point estimates, CI, SE, B, n_boot_used.
#       * CSV summary across combos:
#           "psm_UNIT_ATT_summary_{region_var}.csv" (in `root_out/{region_var}`)
#
# Key conventions :
#   - ATT direction is Control − Treat (positive ⇒ controls larger).
#   - Exact matching on {wday, temp_match, time_of_day, region_var}; caliper on logit-PS.
#   - When MatchIt fails to produce subclasses, we reconstruct pairs via match.matrix.
###############################################################################

suppressPackageStartupMessages({
  library(MatchIt)
  library(data.table)
  library(future)
  library(future.apply)
})

## ============================== USER PATHS =================================
# Replace with your own directories.
root_in  <- "C:/path/to/Seoul_RData"                  # contains seoul_q{3,4,5}.RData
root_out <- "C:/path/to/Matching_Results/PSM"         # outputs grouped by region_var

## ============================= GLOBAL SETTINGS =============================
set.seed(20250810)   # reproducibility for deterministic parts

# Treatment period (inclusive). Adjust if needed.
start_date <- as.Date("2020-01-24")
end_date   <- as.Date("2020-02-09")
date_seq   <- as.Date(seq(start_date, end_date, by = "day"))

# Spatial unit switch: "area" (5 zones) or "dist" (25 districts).
region_var <- "area"

# Temperature matching mode:
#   "tq" → use temp_quantile; "rt" → use rounded temp as exact match factor.
unit_tags   <- c("tq", "rt")

# Symmetric time windows (±k days) for block-wise matching
day_windows <- c(3, 7, 14)

# Quantile-binning variants (input datasets)
q_list      <- c("q3", "q4", "q5")

# Parallelization and bootstrap multiplier distribution
USE_PARALLEL <- TRUE
MULT_ALPHA   <- 1       # 1: Exp(1); >1: Gamma(shape=α, rate=α) with mean 1

## =============================== PSM OPTIONS ===============================
K_MAX        <- 3       # max controls per treated
CALIPER_SD   <- 0.2     # caliper = 0.2 × sd(logit-PS)
REPLACE_CTL  <- TRUE    # allow replacement
TOPN         <- 25      # top lists for influence snapshots

## ============================ OUTPUT STRUCTURE =============================
root_out_var <- file.path(root_out, region_var)
dir.create(root_out_var, recursive = TRUE, showWarnings = FALSE)

## ============================== UTILITIES ==================================
to_num  <- function(x) suppressWarnings(as.numeric(as.character(x)))
safe_min <- function(x){ x <- x[is.finite(x)]; if (!length(x)) return(NA_real_); suppressWarnings(min(x, na.rm=TRUE)) }
safe_max <- function(x){ x <- x[is.finite(x)]; if (!length(x)) return(NA_real_); suppressWarnings(max(x, na.rm=TRUE)) }
logit   <- function(p) qlogis(p)
draw_mult <- function(n){ if (MULT_ALPHA==1) rexp(n,1) else rgamma(n,shape=MULT_ALPHA,rate=MULT_ALPHA) }

# Categorical SMD: before (unweighted) vs after (control weights).
compute_smd_cat <- function(pool_df, covars, w_col="w"){
  out_list <- vector("list", length(covars))
  Trows <- pool_df$treated == 1L
  Crows <- !Trows
  nT <- sum(Trows, na.rm=TRUE); nC <- sum(Crows, na.rm=TRUE)
  wC_sum <- if (w_col %in% names(pool_df)) sum(as.double(pool_df[[w_col]][Crows]), na.rm=TRUE) else NA_real_
  for (j in seq_along(covars)){
    v <- covars[j]; if (!v %in% names(pool_df)) next
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
    } else pC_a <- rep(NA_real_, length(levs))
    den_a <- sqrt((pT_b*(1-pT_b) + pC_a*(1-pC_a))/2)
    smd_a <- (pC_a - pT_b)/den_a
    smd_b[!is.finite(smd_b)] <- NA_real_; smd_a[!is.finite(smd_a)] <- NA_real_
    out_list[[j]] <- data.table(variable=v, level=levs,
                                pT_before=pT_b, pC_before=pC_b, pC_after=pC_a,
                                smd_before=smd_b, smd_after=smd_a,
                                abs_before=abs(smd_b), abs_after=abs(smd_a),
                                nT=nT, nC=nC, wC=wC_sum)
  }
  rbindlist(out_list, use.names=TRUE, fill=TRUE)
}

# Continuous SMD: here only for logit-PS (lp), before vs after (weighted).
compute_smd_cont <- function(df, covars, w_col="w"){
  out <- vector("list", length(covars))
  Trows <- df$treated==1L
  Crows <- df$treated==0L
  for (j in seq_along(covars)){
    v <- covars[j]; if (!v %in% names(df)) next
    xT <- as.numeric(df[[v]][Trows])
    xC <- as.numeric(df[[v]][Crows])
    mT <- mean(xT, na.rm=TRUE); vT <- stats::var(xT, na.rm=TRUE)
    if (!is.null(w_col) && w_col %in% names(df)) {
      wC <- as.numeric(df[[w_col]][Crows]); wC[!is.finite(wC)] <- NA_real_
      sw <- sum(wC, na.rm=TRUE)
      if (is.finite(sw) && sw>0){
        p  <- wC / sw
        mC <- sum(p * xC, na.rm=TRUE)
        vC <- sum(p * (xC - mC)^2, na.rm=TRUE)
      } else { mC <- mean(xC, na.rm=TRUE); vC <- stats::var(xC, na.rm=TRUE) }
    } else { mC <- mean(xC, na.rm=TRUE); vC <- stats::var(xC, na.rm=TRUE) }
    den <- sqrt((vT + vC)/2)
    smd <- (mC - mT)/den
    out[[j]] <- data.table(variable=v, smd=smd, abs=abs(smd),
                           meanT=mT, meanC=mC, sd_pool=den)
  }
  rbindlist(out, use.names=TRUE, fill=TRUE)
}

# Post-match control weights:
#   within each matched set, weight = 1/K (K = #controls in the set);
#   then aggregate to rid-level and scale so that sum_w(control) = n_treated(block).
build_after_weights <- function(matched_dt) {
  ct <- matched_dt[treated==0L, .(w_use = sum(1/ratio_in_set)), by=rid]
  nT <- matched_dt[treated==1L, .N]
  if (!nrow(ct) || !is.finite(nT) || nT<=0) return(data.table(rid=integer(), w_after=numeric()))
  s  <- nT / sum(ct$w_use)
  ct[, w_after := s * w_use][, .(rid, w_after)]
}

# Pair/set contributions: compute d = Control − Treat per matched set.
pair_contrib <- function(md) {
  md[, .(
    d10  = mean(Y10[treated==0L]) - {t <- Y10[treated==1L]; if (length(t)) t[1L] else NA_real_},
    d25  = mean(Y25[treated==0L]) - {t <- Y25[treated==1L]; if (length(t)) t[1L] else NA_real_},
    ridT = {t <- rid[treated==1L]; if (length(t)) t[1L] else NA_integer_}
  ), by=.(set_id)]
}

# Rebuild matched sets from MatchIt's match.matrix (robust when subclass is absent).
# Requires rownames(pool) = rid so match.matrix rows index treated rid explicitly.
build_md_from_match_matrix <- function(m, pool) {
  mm <- m$match.matrix
  if (is.null(mm)) return(NULL)
  pool_sub <- as.data.table(pool[, c("rid","treated","PM10","PM25")])
  rows <- list()
  t_ids <- rownames(mm)
  for (i in seq_along(t_ids)) {
    ridT_chr <- t_ids[i]
    ridT <- suppressWarnings(as.numeric(ridT_chr))
    t_row <- pool_sub[rid == ridT][1]
    if (nrow(t_row)==0) next
    t_row[, `:=`(set_id = ridT, Y10 = PM10, Y25 = PM25)]
    t_row <- t_row[, .(set_id, rid, treated, Y10, Y25)]
    
    ctr_ids_chr <- as.character(mm[i, ])
    ctr_ids_chr <- ctr_ids_chr[!is.na(ctr_ids_chr) & nzchar(ctr_ids_chr)]
    if (length(ctr_ids_chr)==0) next
    ctr_ids <- suppressWarnings(as.numeric(ctr_ids_chr))
    c_rows  <- pool_sub[rid %in% ctr_ids]
    if (!nrow(c_rows)) next
    c_rows[, `:=`(set_id = ridT, Y10 = PM10, Y25 = PM25)]
    c_rows <- c_rows[, .(set_id, rid, treated, Y10, Y25)]
    
    rows[[length(rows)+1L]] <- rbindlist(list(t_row, c_rows))
  }
  if (!length(rows)) return(NULL)
  md <- rbindlist(rows, use.names=TRUE, fill=TRUE)
  rat <- md[, .(ratio_in_set = sum(treated==0L)), by=.(set_id)]
  md <- merge(md, rat, by="set_id", all.x=TRUE)
  md[]
}

## ============================ EXPERIMENT GRID ==============================
combos <- data.table::CJ(q=q_list, unit_tag=unit_tags, day_window=day_windows)

## ============================== CORE ROUTINE ===============================
run_combo_psm <- function(idx){
  q        <- combos$q[idx]
  unit_tag <- combos$unit_tag[idx]
  day_w    <- combos$day_window[idx]
  
  # ---- Load & prepare data for this q ----
  load(file.path(root_in, sprintf("seoul_%s.RData", q)))
  seoul <- get(sprintf("seoul_%s", q))
  seoul <- as.data.frame(seoul, stringsAsFactors=FALSE)
  seoul$date <- as.Date(seoul$date)
  if (!("rid" %in% names(seoul))) seoul$rid <- seq_len(nrow(seoul))
  
  # Defensive numeric coercion
  num_cols0 <- intersect(c("PM10","PM25","temp","wind_sp_binary","hour","day","month","year"),
                         names(seoul))
  for (nm in num_cols0) seoul[[nm]] <- to_num(seoul[[nm]])
  
  # Temperature matching factor
  if (unit_tag=="tq") seoul$temp_match <- as.character(seoul$temp_quantile)
  else                seoul$temp_match <- as.character(round(seoul$temp))
  
  # Covariate set (assumes wday already exists in data)
  covars_all <- c("wday","temp_match","rain_binary","CO_quantile","O3_quantile",
                  "wind_sp_binary","wind_dir_8","time_of_day", region_var)
  covars <- intersect(covars_all, names(seoul))
  ps_form <- as.formula(paste("treated ~", paste(covars, collapse=" + ")))
  
  # Containers
  DALL_list <- list(); smd_list <- list(); stab_list <- list()
  infl_list <- list(); drop_log <- list(); reuse_raw_list <- list()
  
  # ============================ BLOCK LOOP ==================================
  for (td in date_seq){
    td <- as.Date(td); win <- as.Date(td) + (-day_w:day_w)
    
    trt <- seoul[seoul$date==td & seoul$year==2020, , drop=FALSE]
    if (!nrow(trt)) { drop_log[[length(drop_log)+1L]] <- data.table(date=td, reason="no_treated"); next }
    trt$treated <- 1L
    
    mdays <- format(win, "%m-%d")
    ctrl <- seoul[format(seoul$date,"%m-%d") %in% mdays & seoul$year!=2020, , drop=FALSE]
    ctrl <- ctrl[ ctrl[[region_var]] %in% unique(trt[[region_var]]), , drop=FALSE]
    if (!nrow(ctrl)) { drop_log[[length(drop_log)+1L]] <- data.table(date=td, reason="no_controls"); next }
    ctrl$treated <- 0L
    
    pool <- rbind(trt, ctrl)
    
    # Complete-case on covariates used in PS model
    cc_idx <- stats::complete.cases(pool[, covars])
    if (!all(cc_idx)) drop_log[[length(drop_log)+1L]] <- data.table(date=td, reason="covariate_missing", n_drop=sum(!cc_idx))
    pool <- pool[cc_idx, , drop=FALSE]
    if (!nrow(pool)) { drop_log[[length(drop_log)+1L]] <- data.table(date=td, reason="all_dropped_cc"); next }
    
    # Ensure factors for covariates; numeric outcomes; treated integer
    for (nm in covars) pool[[nm]] <- as.factor(pool[[nm]])
    pool$PM10 <- to_num(pool$PM10); pool$PM25 <- to_num(pool$PM25); pool$treated <- as.integer(pool$treated)
    
    # Critical: rownames(pool) = rid so match.matrix rows map to treated rid
    rownames(pool) <- as.character(pool$rid)
    
    # ---- PS modeling (logit link) ----
    ps_fit <- tryCatch(glm(ps_form, data=pool, family=binomial()), error=function(e) NULL)
    if (is.null(ps_fit)) { drop_log[[length(drop_log)+1L]] <- data.table(date=td, reason="glm_fail"); next }
    pool$ps <- pmin(pmax(fitted(ps_fit), 1e-6), 1-1e-6)
    pool$lp <- logit(pool$ps)
    sd_lp   <- stats::sd(pool$lp, na.rm=TRUE)
    cali    <- if (is.finite(sd_lp)) CALIPER_SD * sd_lp else 0.2
    
    # ---- PSM via MatchIt (distance = logit-PS) ----
    m <- tryCatch(
      matchit(
        treated ~ 1, data = pool, method = "nearest",
        distance = pool$lp, replace = REPLACE_CTL, ratio = K_MAX,
        exact = intersect(c("wday","temp_match","time_of_day", region_var), names(pool)),
        caliper = cali, estimand = "ATT"
      ),
      error=function(e) NULL
    )
    if (is.null(m) || is.null(m$match.matrix)) {
      drop_log[[length(drop_log)+1L]] <- data.table(date=td, reason="matchit_fail_or_nomatch")
      next
    }
    
    # ---- Rebuild matched sets; compute pair contributions (Control − Treat) ----
    md <- build_md_from_match_matrix(m, pool)
    if (is.null(md) || !nrow(md) || md[treated==1L, .N]==0L) {
      drop_log[[length(drop_log)+1L]] <- data.table(date=td, reason="build_md_fail")
      next
    }
    
    # ---- SMD (before/after): categorical and continuous (lp) ----------------
    smd_b_cat <- compute_smd_cat(pool, covars, w_col="w")[, `:=`(date=td, stage="before")]
    
    w_after <- build_after_weights(md)                        # (rid, w_after); sum to nT
    reuse_raw_list[[length(reuse_raw_list)+1L]] <- data.table(date=td, w_after)
    
    pool_after <- merge(pool[, c("rid","treated", covars, "lp")], w_after, by="rid", all.x=TRUE)
    pool_after$w <- ifelse(pool_after$treated==0L, pool_after$w_after, NA_real_)
    matched_tids <- unique(md[treated==1L, rid])
    
    setDT(pool_after)
    pool_after <- pool_after[(treated==1L & rid %in% matched_tids) | (treated==0L & !is.na(w))]
    
    smd_a_cat <- compute_smd_cat(as.data.frame(pool_after), covars, w_col="w")[, `:=`(date=td, stage="after")]
    smd_b_cont <- compute_smd_cont(pool[, c("treated","lp")], covars="lp")[, `:=`(date=td, stage="before")]
    smd_a_cont <- compute_smd_cont(as.data.frame(pool_after)[, c("treated","lp","w")], covars="lp", w_col="w")[, `:=`(date=td, stage="after")]
    
    smd_list[[length(smd_list)+1L]] <- rbind(
      cbind(type="categorical", smd_b_cat,  region_key = region_var),
      cbind(type="categorical", smd_a_cat,  region_key = region_var),
      cbind(type="continuous_before",  smd_b_cont, region_key = region_var),
      cbind(type="continuous_after",   smd_a_cont, region_key = region_var),
      fill=TRUE
    )
    
    # ---- Stability & influence snapshots (after weights) --------------------
    if (nrow(w_after)){
      sumw <- sum(w_after$w_after); sumw2 <- sum((w_after$w_after)^2)
      ess  <- if (sumw2>0) (sumw^2)/sumw2 else NA_real_
      stab <- data.table(
        date=td, ESS_C_after=ess,
        max_w_after = safe_max(w_after$w_after),
        max_w_rel_after = if (sumw>0) safe_max(w_after$w_after)/sumw else NA_real_,
        nT = md[treated==1L, .N], region_key = region_var
      )
      stab_list[[length(stab_list)+1L]] <- stab
      qv <- as.list(quantile(w_after$w_after, probs=c(.5,.9,.95,.99), na.rm=TRUE))
      infl_list[[length(infl_list)+1L]] <- data.table(date=td, metric="w_after_quantile",
                                                      p50=qv[[1]], p90=qv[[2]], p95=qv[[3]], p99=qv[[4]],
                                                      region_key = region_var)
    }
    
    # ---- Pair contributions (Control − Treat) -------------------------------
    dsub <- pair_contrib(md)[, `:=`(date = td, region_key = region_var)]
    DALL_list[[length(DALL_list)+1L]] <- dsub
  } # end block loop
  
  # If nothing matched across all blocks, return sentinel summary
  if (!length(DALL_list)) {
    return(list(q=q, unit_tag=unit_tag, window=day_w, region_key=region_var,
                error=TRUE, message="no valid blocks"))
  }
  
  # ======================== AGGREGATE ACROSS BLOCKS ==========================
  DALL <- rbindlist(DALL_list)
  point10 <- mean(DALL$d10, na.rm=TRUE)   # overall ATT (PM10), Control − Treat
  point25 <- mean(DALL$d25, na.rm=TRUE)   # overall ATT (PM25), Control − Treat
  
  # ======================= MULTIPLIER BOOTSTRAP (unit) ========================
  B <- 1000L
  att10_B <- att25_B <- numeric(B)
  for (b in seq_len(B)){
    m <- draw_mult(nrow(DALL))                # one weight per matched set
    wsum <- sum(m); if (!is.finite(wsum) || wsum<=0) { att10_B[b] <- NA; att25_B[b] <- NA; next }
    att10_B[b] <- sum(m * DALL$d10, na.rm=TRUE) / wsum
    att25_B[b] <- sum(m * DALL$d25, na.rm=TRUE) / wsum
  }
  ci10 <- stats::quantile(att10_B, c(.025,.975), na.rm=TRUE, names=FALSE)
  ci25 <- stats::quantile(att25_B, c(.025,.975), na.rm=TRUE, names=FALSE)
  se10 <- stats::sd(att10_B, na.rm=TRUE)
  se25 <- stats::sd(att25_B, na.rm=TRUE)
  n_boot_used <- sum(is.finite(att10_B) & is.finite(att25_B))
  
  # Diagnostics tables
  smd_table      <- if (length(smd_list)) rbindlist(smd_list, fill=TRUE) else data.table()
  stab_table     <- if (length(stab_list)) rbindlist(stab_list, fill=TRUE) else data.table()
  influence_table<- if (length(infl_list)) rbindlist(infl_list, fill=TRUE) else data.table()
  
  # Top influences
  top_control_reuse <- if (length(reuse_raw_list)) {
    reuse_all <- rbindlist(reuse_raw_list, fill=TRUE)
    reuse_all[, .(w_after_sum = sum(w_after)), by=rid][order(-w_after_sum)][1:min(TOPN, .N)]
  } else data.table()
  top_pairs_d10 <- DALL[order(-abs(d10))][1:min(TOPN, .N), .(set_id, ridT, d10)]
  top_pairs_d25 <- DALL[order(-abs(d25))][1:min(TOPN, .N), .(set_id, ridT, d25)]
  
  # ================================ SAVE ====================================
  rds_name <- sprintf("psm_UNIT_ATT_boot_%s_%s_%s_d%d_WITH_SMD_DIAG.RData",
                      q, unit_tag, region_var, day_w)
  rds_path <- file.path(root_out_var, rds_name)
  
  si <- utils::capture.output(sessionInfo())
  save(DALL, att10_B, att25_B, smd_table, stab_table, influence_table,
       top_control_reuse, top_pairs_d10, top_pairs_d25,
       drop_log, K_MAX, CALIPER_SD, REPLACE_CTL, MULT_ALPHA, si, region_var,
       point10, point25, ci10, ci25, se10, se25, B, n_boot_used,
       file = rds_path)
  
  list(
    q=q, unit_tag=unit_tag, window=day_w, region_key=region_var,
    PM10=point10, PM25=point25,
    lo10=ci10[1], hi10=ci10[2], lo25=ci25[1], hi25=ci25[2],
    se10=se10, se25=se25, n_boot=n_boot_used,
    ESS_after_mean = if (nrow(stab_table)) mean(stab_table$ESS_C_after, na.rm=TRUE) else NA_real_,
    error=FALSE, message=NA_character_
  )
}

# Safe wrapper to avoid aborting the whole run on one failure
safe_run_psm <- function(i){
  tryCatch(run_combo_psm(i),
           error=function(e) list(q=combos$q[i], unit_tag=combos$unit_tag[i], window=combos$day_window[i],
                                  region_key=region_var, error=TRUE, message=conditionMessage(e)))
}

## ================================ DRIVER ===================================
do_all_psm <- function(){
  if (USE_PARALLEL) {
    workers <- max(1, parallel::detectCores() - 2)
    plan(multisession, workers=workers); on.exit(plan(sequential), add=TRUE)
    withCallingHandlers(
      future_lapply(seq_len(nrow(combos)), safe_run_psm, future.seed=TRUE),
      warning=function(w){
        if (grepl("Future .* was canceled", conditionMessage(w))) invokeRestart("muffleWarning")
      }
    )
  } else {
    lapply(seq_len(nrow(combos)), safe_run_psm)
  }
}

psm_collect <- do_all_psm()
psm_table   <- data.table::rbindlist(psm_collect, fill=TRUE)

out_csv <- file.path(root_out_var, sprintf("psm_UNIT_ATT_summary_%s.csv", region_var))
data.table::fwrite(psm_table, out_csv)

print(psm_table)
