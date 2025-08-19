###############################################################################
# Title   : Block-wise Exact Matching (±k days) with Temperature tq/rt option
#           + Post-match Balance (SMD) + Outcome Differences (Control − Treat)
# Author  : (Your Name)
# Version : 1.1
# Date    : 2025-08-19
#
# Purpose :
#   - For each combination of {quantile version q ∈ {q3,q4,q5}, unit_tag ∈ {tq,rt},
#     day_window ∈ {3,7,14}}, perform block-wise exact matching:
#       * Treated: dates between start_date and end_date in year 2020
#       * Controls: same calendar month-day within ±k days, years ≠ 2020,
#         same spatial unit (region_var), and exact equality on core covariates.
#       * Temperature handling:
#           - "tq": match on temp_quantile
#           - "rt": match on round(temp)
#   - Compute Control − Treat differences for PM10/PM25 per matched control row.
#   - Compute SMDs (before vs after matching) on categorical covariates via
#     dummy coding; store all diagnostics and matched datasets.
#
# Notes   :
#   - This is a *deterministic exact-matching* routine (not propensity-based).
#   - Variables used for exact matching will often have SMD≈0 after matching
#     by construction. That is expected and indicates design-based balance.
#
# Inputs  :
#   - RData files with pre-processed lists (one per quantile scheme):
#       seoul_q3.RData  (object: seoul_q3)
#       seoul_q4.RData  (object: seoul_q4)
#       seoul_q5.RData  (object: seoul_q5)
#     Each object must include at least:
#       date (Date), year, month, day, hour,
#       PM10, PM25, temp, temp_quantile,
#       rain_binary, wind_sp_binary, wind_dir_8, time_of_day,
#       CO_quantile, O3_quantile, wday,
#       area (5 zones) and dist (25 districts).
#
# Outputs :
#   - RData per combo: "cem_{region_var}_{unit_tag}_{q}_d{day_window}.RData"
#       containing: matched_df (controls only with attached target), final_matched_df
#       (controls + matched treated), smd_before, smd_after,
#       and simple ATT-like summaries (att_mean10/25, att_sd10/25).
#   - Optional CSV exports for SMDs can be added if desired.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(lubridate)
  library(progress)
  library(ggplot2)
})

## ============================== USER PATHS =================================
# Replace with your own directories.
root_in   <- "C:/path/to/Seoul_RData"    # where seoul_q3/4/5.RData live
root_out  <- "C:/path/to/ExactMatch_Out" # where outputs will be saved
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

## ============================= GLOBAL SETTINGS =============================
# Load pre-processed datasets (one object per file)
load(file.path(root_in, "seoul_q3.RData"))  # provides object `seoul_q3`
load(file.path(root_in, "seoul_q4.RData"))  # provides object `seoul_q4`
load(file.path(root_in, "seoul_q5.RData"))  # provides object `seoul_q5`

seoul_list   <- list(q3 = seoul_q3, q4 = seoul_q4, q5 = seoul_q5)

unit_tags    <- c("tq", "rt")           # tq = temp_quantile; rt = rounded temperature
day_windows  <- c(3, 7, 14)             # symmetric windows (±k days)
region_var   <- "area"                  # "area" (5 zones) or "dist" (25 districts)
start_date   <- as.Date("2020-01-24")   # treatment start (lockdown)
end_date     <- as.Date("2020-02-09")   # treatment end

# Filename prefix (kept as 'cem_' for pipeline compatibility; change as needed)
file_prefix  <- "cem"

## ================================ HELPERS ==================================
# Build a dummy (0/1) design matrix for categorical balance diagnostics.
# Numeric binaries are auto-treated as factors to ensure proper dummy coding.
get_dummy_matrix <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (!length(vars)) return(data.frame())
  df_local <- df[, vars, drop = FALSE]
  df_local[] <- lapply(df_local, function(x) {
    # Convert to factor for consistent dummy creation (binary/ordinal ok)
    if (is.factor(x)) x else as.factor(x)
  })
  form <- as.formula(paste("~", paste(vars, collapse = " + ")))
  mm <- model.matrix(form, data = df_local)
  # Drop intercept if present
  if (colnames(mm)[1] == "(Intercept)") mm <- mm[, -1, drop = FALSE]
  as.data.frame(mm, stringsAsFactors = FALSE)
}

# Standardized Mean Difference for numeric columns x1, x2
diff_smd <- function(x1, x2) {
  m1 <- mean(x1, na.rm = TRUE); m2 <- mean(x2, na.rm = TRUE)
  s1 <- stats::sd(x1, na.rm = TRUE); s2 <- stats::sd(x2, na.rm = TRUE)
  sp <- sqrt((s1^2 + s2^2) / 2)
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (m1 - m2) / sp
}

# Safe mean/SD for vectors with NAs
smean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
ssd   <- function(x) if (all(is.na(x))) NA_real_ else stats::sd(x, na.rm = TRUE)

## =============================== MAIN LOOP =================================
for (q in names(seoul_list)) {
  seoul <- seoul_list[[q]]
  # Light coercions & checks
  seoul$date <- as.Date(seoul$date)
  if (!("wday" %in% names(seoul))) stop("`wday` must exist in input data.")
  if (!(region_var %in% names(seoul))) stop(sprintf("`%s` not found in data.", region_var))
  
  # Predefine covariates used in exact matching / SMD
  covariates_all <- c("temp_quantile", "rain_binary", "CO_quantile", "O3_quantile",
                      "wind_sp_binary", "wind_dir_8", "wday", "time_of_day")
  
  for (unit_tag in unit_tags) {
    for (day_window in day_windows) {
      
      # --------------------- Define treated/control pools --------------------
      target_data     <- seoul %>% filter(date >= start_date & date <= end_date)
      comparison_data <- seoul %>% filter(!(date >= start_date & date <= end_date))
      
      if (!nrow(target_data) || !nrow(comparison_data)) next
      
      # Progress bar
      pb <- progress::progress_bar$new(total = nrow(target_data),
                                       format = "[:bar] :percent eta: :eta  (q={q} {unit_tag} ±{k}d)")
      pb$tick(0, tokens = list(q = q, unit_tag = unit_tag, k = day_window))
      
      matching_results <- vector("list", nrow(target_data))  # per treated row controls
      
      # ====================== Block-wise exact matching ======================
      for (i in seq_len(nrow(target_data))) {
        target_row <- target_data[i, , drop = FALSE]
        
        # ±k days window by month-day (e.g., "01-25"), independent of year
        date_seq_md <- format(seq(target_row$date - days(day_window),
                                  target_row$date + days(day_window), by = "day"), "%m-%d")
        
        # Base exact constraints (excluding temperature at this step)
        base_match <- comparison_data %>%
          filter(format(date, "%m-%d") %in% date_seq_md,
                 .data[["time_of_day"]] == target_row$time_of_day,
                 .data[[region_var]]    == target_row[[region_var]],
                 wday            == target_row$wday,
                 CO_quantile     == target_row$CO_quantile,
                 O3_quantile     == target_row$O3_quantile,
                 wind_sp_binary  == target_row$wind_sp_binary,
                 wind_dir_8      == target_row$wind_dir_8,
                 rain_binary     == target_row$rain_binary)
        
        # Temperature matching mode
        matched_data <-
          if (unit_tag == "tq") {
            base_match %>% filter(temp_quantile == target_row$temp_quantile)
          } else {
            base_match %>% filter(round(temp) == round(target_row$temp))
          }
        
        # Keep if any controls found
        if (nrow(matched_data) > 0) {
          matched_data$matched_index <- seq_len(nrow(matched_data))  # per-target ranking (optional)
          matched_data$target_index  <- i
          matching_results[[i]] <- matched_data
        }
        
        pb$tick(tokens = list(q = q, unit_tag = unit_tag, k = day_window))
      } # end treated loop
      
      # Skip if no matches at all
      non_empty <- which(lengths(matching_results) > 0L)
      if (!length(non_empty)) next
      
      # Merge all matched controls
      matched_df <- data.table::rbindlist(matching_results[non_empty], use.names = TRUE, fill = TRUE)
      
      # Attach treated outcomes to each matched control row (for differences)
      target_data_labeled <- target_data %>% mutate(treated = 1L, target_index = row_number())
      setDT(matched_df)
      matched_df <- merge(
        matched_df,
        target_data_labeled[, .(target_index, PM10_target = PM10, PM25_target = PM25)],
        by = "target_index",
        all.x = TRUE
      )
      
      # Outcome differences (Control − Treat)
      matched_df[, `:=`(
        PM10_diff = PM10 - PM10_target,
        PM25_diff = PM25 - PM25_target
      )]
      
      # Reconstruct *treated* rows that actually found controls
      target_matched <- unique(matched_df[, .(target_index)])[
        target_data_labeled, on = .(target_index)
      ]
      
      # Final matched dataset = controls (with attached targets) + matched treated rows
      final_matched_df <- rbindlist(list(
        matched_df[, c(names(matched_df), "treated") := .(matched_df, 0L)][[]],
        target_matched
      ), use.names = TRUE, fill = TRUE)
      
      # ======================== Balance diagnostics ==========================
      # Only use covariates that exist and vary in the (sub)samples
      covariates <- intersect(covariates_all, names(final_matched_df))
      setDT(final_matched_df)
      
      # After matching
      if (!all(c("treated") %in% names(final_matched_df))) {
        stop("`treated` column is missing in final_matched_df.")
      }
      # Keep covariates with >1 level in both groups
      valid_cov_after <- covariates[
        vapply(covariates, function(v) {
          all(table(final_matched_df$treated, final_matched_df[[v]]) > 0) &&
            length(unique(final_matched_df[[v]])) > 1
        }, logical(1))
      ]
      
      X_treat_after   <- get_dummy_matrix(final_matched_df[treated == 1], valid_cov_after)
      X_control_after <- get_dummy_matrix(final_matched_df[treated == 0], valid_cov_after)
      smd_after <- if (ncol(X_treat_after) && ncol(X_control_after)) {
        sapply(intersect(names(X_treat_after), names(X_control_after)),
               function(v) diff_smd(X_control_after[[v]], X_treat_after[[v]]))
      } else numeric(0)
      
      # Before matching (full treated vs non-treated outside window)
      unmatched_control <- seoul %>%
        filter(!(date >= start_date & date <= end_date)) %>%
        mutate(treated = 0L)
      unmatched_treated <- target_data_labeled
      
      valid_cov_before <- intersect(covariates_all, names(unmatched_treated))
      valid_cov_before <- intersect(valid_cov_before, names(unmatched_control))
      valid_cov_before <- valid_cov_before[
        vapply(valid_cov_before, function(v) {
          (length(unique(unmatched_treated[[v]]))  > 1) &&
            (length(unique(unmatched_control[[v]])) > 1)
        }, logical(1))
      ]
      
      X_treat_before   <- get_dummy_matrix(unmatched_treated,  valid_cov_before)
      X_control_before <- get_dummy_matrix(unmatched_control, valid_cov_before)
      smd_before <- if (ncol(X_treat_before) && ncol(X_control_before)) {
        sapply(intersect(names(X_treat_before), names(X_control_before)),
               function(v) diff_smd(X_control_before[[v]], X_treat_before[[v]]))
      } else numeric(0)
      
      # =========================== Simple summaries ==========================
      # ATT-like summaries (mean of Control − Treat differences over matched rows)
      att_mean10 <- smean(matched_df$PM10_diff)
      att_mean25 <- smean(matched_df$PM25_diff)
      att_sd10   <- ssd(matched_df$PM10_diff)
      att_sd25   <- ssd(matched_df$PM25_diff)
      
      # ================================ SAVE =================================
      save_dir <- file.path(root_out, region_var)
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
      
      save_path <- file.path(
        save_dir,
        sprintf("%s_%s_%s_%s_d%d.RData", file_prefix, region_var, unit_tag, q, day_window)
      )
      
      # Package SMDs as named vectors for convenience
      smd_before_vec <- smd_before
      smd_after_vec  <- smd_after
      
      si <- utils::capture.output(sessionInfo())
      save(matched_df, final_matched_df,
           smd_before_vec, smd_after_vec,
           att_mean10, att_mean25, att_sd10, att_sd25,
           start_date, end_date, region_var, unit_tag, day_window, q, si,
           file = save_path)
      
      message(sprintf("[Saved] %s", save_path))
    } # day_window
  }   # unit_tag
}     # q
