###############################################################################
# Title   : CBPS Balance Diagnostics — 2×3 Love-Plot Panels + CSV Exports
# Author  : Juwon Jung
# Version : 1.1
# Date    : 2025-08-19
#
# Purpose :
#   - Load CBPS diagnostic bundles saved per combination of:
#       quantile scheme q ∈ {q3, q4, q5},
#       time window d ∈ {±3, ±7, ±14},
#       region_var ∈ {"dist", "area"},
#       temperature mode ∈ {"rt","tq"} (embedded in filenames).
#   - Summarize variable-level |SMD| across dates (median or nT-weighted mean)
#     into a single value per variable (max across levels).
#   - Plot 2×3 Love-plot panels (rows: region, columns: d ∈ {3,7,14}) for each q.
#   - Save condition-specific CSVs (wide & tidy) and optional raw SMD tables.
#
# Inputs  :
#   - RData files created by your CBPS pipeline, each containing `smd_table`
#     with columns at least:
#       variable, level, abs_before, abs_after, date, nT  (and possibly others)
#
# Outputs :
#   - Panels:   "CBPS_SMD_Panel_{q}_{TEMP_MODE}.png"
#   - CSV wide: base_dir/SMD_Values/SMD_values_{region}_{unit}_{q}_d{d}.csv
#   - CSV tidy: base_dir/SMD_Values_tidy/SMD_tidy_{region}_{unit}_{q}_d{d}.csv
#   - CSV raw (optional): base_dir/SMD_Values_raw/SMD_raw_{region}_{unit}_{q}_d{d}.csv
#
# Notes   :
#   - TEMP_MODE controls label text; files are auto-detected for both rt/tq in case
#     your RData filenames encode the mode. Set TEMP_MODE to the target mode for
#     labeling/filenames of CSV/PNG outputs.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(ggpubr)
})

## ============================== USER PATHS =================================
# Change only this root directory to where your CBPS results live.
base_dir  <- "C:/path/to/Matching_Results/CBPS"

## ============================= GLOBAL SETTINGS =============================
wins      <- c(3, 7, 14)                # ±k-day windows
qs        <- c("q3", "q4", "q5")        # quantile schemes
DATE_AGG  <- "weighted_mean"            # "median" or "weighted_mean"
TEMP_MODE <- "rt"                       # "rt" or "tq" (affects labels/outputs only)

## ======================== FILE DISCOVERY / LOADING =========================
# Candidate file patterns: supports both with/without explicit region suffix, and rt/tq.
cand_paths <- function(q, region, win) {
  dirp <- file.path(base_dir, region)
  file.path(dirp, c(
    sprintf("cbps_UNIT_ATT_boot_%s_rt_d%d_%s_WITH_SMD_DIAG.RData", q, win, region),
    sprintf("cbps_UNIT_ATT_boot_%s_tq_d%d_%s_WITH_SMD_DIAG.RData", q, win, region),
    sprintf("cbps_UNIT_ATT_boot_%s_rt_d%d_WITH_SMD_DIAG.RData", q, win),
    sprintf("cbps_UNIT_ATT_boot_%s_tq_d%d_WITH_SMD_DIAG.RData", q, win)
  ))
}

# Load the first existing file from candidates and return smd_table as data.table
load_smd_from_first_exist <- function(paths) {
  for (p in paths) if (file.exists(p)) {
    e <- new.env()
    load(p, envir = e)
    if (exists("smd_table", envir = e)) {
      smd <- as.data.table(e$smd_table)
      attr(smd, "src_file") <- p
      return(smd)
    }
  }
  NULL
}

## =============================== SUMMARIZATION ==============================
# Summarize SMD across dates:
#  - For each variable-level: {abs_before, abs_after} aggregated by date (median or nT-weighted mean).
#  - Then collapse to variable-level by taking max across levels.
summarize_smd <- function(smd, date_agg = c("median","weighted_mean")) {
  date_agg <- match.arg(date_agg)
  # Defensive columns
  if (!("date" %in% names(smd))) smd[, date := NA]
  if (!("nT"   %in% names(smd))) smd[, nT := 1]
  
  if (date_agg == "median") {
    lvl <- smd[, .(
      abs_b = median(abs_before, na.rm = TRUE),
      abs_a = median(abs_after,  na.rm = TRUE),
      n_dates = uniqueN(date)
    ), by = .(variable, level)]
  } else {
    lvl <- smd[, .(
      abs_b = {
        w <- nT; if (!any(is.finite(w))) w[] <- 1
        weighted.mean(abs_before, w, na.rm = TRUE)
      },
      abs_a = {
        w <- nT; if (!any(is.finite(w))) w[] <- 1
        weighted.mean(abs_after,  w, na.rm = TRUE)
      },
      n_dates = uniqueN(date)
    ), by = .(variable, level)]
  }
  
  # Clean non-finite values and ensure zeros for missing levels
  lvl[, c("abs_b","abs_a") := lapply(.SD, function(z) ifelse(is.finite(z), z, NA_real_)),
      .SDcols = c("abs_b","abs_a")]
  lvl[, abs_b := fifelse(is.na(abs_b), 0, abs_b)]
  lvl[, abs_a := fifelse(is.na(abs_a), 0, abs_a)]
  
  # Collapse to variable-level (max across levels)
  lvl[, .(
    before   = max(abs_b, na.rm = TRUE),
    after    = max(abs_a, na.rm = TRUE),
    n_levels = .N,
    n_dates  = max(n_dates, na.rm = TRUE)
  ), by = variable]
}

## ============================== LABEL HELPERS ==============================
label_vars <- function(v, temp_mode = c("rt","tq"), q_tag = c("q3","q4","q5")){
  temp_mode <- match.arg(temp_mode); q_tag <- match.arg(q_tag)
  qN <- as.integer(sub("^q", "", q_tag))
  v <- as.character(v)
  v[v == "wind_sp_binary"] <- "Wind speed (binary)"
  v[v == "rain_binary"]    <- "Precipitation occurrence"
  v[v == "time_of_day"]    <- "Time of day"
  v[v == "wind_dir_8"]     <- "Wind directions"
  v[v == "wday"]           <- "Weekday"
  v[v == "dist"]           <- "Administrative districts"
  v[v == "area"]           <- "Residential zones"
  v[v == "temp_match"]     <- if (temp_mode == "rt") "Temperature (rounded)"
  else sprintf("Temperature (quantiled = %d)", qN)
  v[v == "CO_quantile"]    <- sprintf("CO (quantiled = %d)",  qN)
  v[v == "O3_quantile"]    <- sprintf("O3 (quantiled = %d)",  qN)
  v
}

## ============================== LOVE PLOT ==================================
build_love_plot <- function(var_df, title_text, xlim_max = NULL, temp_mode = TEMP_MODE, q_tag){
  df_long <- melt(var_df, id.vars = "variable",
                  measure.vars = c("before","after"),
                  variable.name = "stage", value.name = "abs_smd")
  df_long[, stage := factor(stage, levels = c("before","after"),
                            labels = c("Before","After"))]
  df_long[, var_show := label_vars(variable, temp_mode, q_tag)]
  df_long[, var_show := fct_reorder(var_show, abs_smd, .fun = function(z) max(z, na.rm = TRUE))]
  
  col_before <- "#1b7837"; col_after <- "#d73027"
  xmax <- if (is.null(xlim_max)) max(df_long$abs_smd, na.rm = TRUE) else xlim_max
  if (!is.finite(xmax)) xmax <- 0.2
  xmax_break <- ceiling(xmax * 10) / 10
  breaks_seq <- seq(0, xmax_break, by = 0.1)
  
  ggplot(df_long, aes(x = abs_smd, y = var_show)) +
    geom_point(aes(color = stage), shape = 16, size = 2.6, alpha = 0.95) +
    geom_vline(xintercept = 0.1, linetype = "dashed", linewidth = 0.6, color = "#555555") +
    scale_color_manual(values = c("Before" = col_before, "After" = col_after)) +
    scale_x_continuous(limits = c(0, xmax_break), breaks = breaks_seq, minor_breaks = NULL) +
    labs(title = title_text, x = "|SMD|", y = NULL, color = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank(),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
}

## =========================== CSV EXPORT OPTIONS ============================
SAVE_CSV_VALUES  <- TRUE   # wide format (per variable: before/after)
SAVE_CSV_TIDY    <- TRUE   # long format (stage = before/after)
SAVE_RAW_SMD_CSV <- FALSE  # optional: raw variable-level SMD by date/level

csv_dir_values <- file.path(base_dir, "SMD_Values")
csv_dir_tidy   <- file.path(base_dir, "SMD_Values_tidy")
csv_dir_raw    <- file.path(base_dir, "SMD_Values_raw")
if (SAVE_CSV_VALUES) dir.create(csv_dir_values, showWarnings = FALSE, recursive = TRUE)
if (SAVE_CSV_TIDY)   dir.create(csv_dir_tidy,   showWarnings = FALSE, recursive = TRUE)
if (SAVE_RAW_SMD_CSV)dir.create(csv_dir_raw,    showWarnings = FALSE, recursive = TRUE)

# Accumulators for master CSVs
.all_values <- list()
.all_tidy   <- list()

## =============================== MAIN PANELS ===============================
for (q in qs) {
  # 6 panels: (dist, d=3/7/14) + (area, d=3/7/14)
  plots <- vector("list", 6L)
  names(plots) <- c("dist_3","dist_7","dist_14","area_3","area_7","area_14")
  
  xmax_global <- 0
  k <- 1L
  
  for (reg in c("dist","area")) {
    reg_label <- if (reg == "dist") "Administrative districts" else "Residential zones"
    
    for (win in wins) {
      smd <- load_smd_from_first_exist(cand_paths(q, reg, win))
      if (is.null(smd)) {
        # Placeholder for missing inputs
        plots[[k]] <- ggplot() + theme_void() +
          labs(title = sprintf("%s — ±%d time window (missing)", reg_label, win))
      } else {
        var_df <- summarize_smd(smd, DATE_AGG)
        xmax_global <- max(xmax_global, var_df[, max(c(before, after), na.rm = TRUE)], na.rm = TRUE)
        
        ## ---------------------------- CSV EXPORTS ---------------------------
        base_tag <- sprintf("%s_%s_%s_d%d", reg, TEMP_MODE, q, win)
        
        # Wide (values)
        if (SAVE_CSV_VALUES) {
          values_out <- copy(var_df)[
            , .(variable, before, after, n_levels, n_dates)
          ][
            , `:=`(variable_label = label_vars(variable, TEMP_MODE, q),
                   region = reg, unit = TEMP_MODE, q = q, d_win = win)
          ]
          fwrite(values_out, file.path(csv_dir_values, sprintf("SMD_values_%s.csv", base_tag)))
          .all_values[[length(.all_values) + 1L]] <- values_out
        }
        
        # Tidy (long)
        if (SAVE_CSV_TIDY) {
          tidy <- melt(var_df, id.vars = "variable",
                       measure.vars = c("before","after"),
                       variable.name = "stage", value.name = "SMD")
          tidy[, `:=`(variable_label = label_vars(variable, TEMP_MODE, q),
                      region = reg, unit = TEMP_MODE, q = q, d_win = win)]
          fwrite(tidy, file.path(csv_dir_tidy, sprintf("SMD_tidy_%s.csv", base_tag)))
          .all_tidy[[length(.all_tidy) + 1L]] <- tidy
        }
        
        # Raw (optional)
        if (SAVE_RAW_SMD_CSV) {
          raw_keep <- intersect(c("variable","level","abs_before","abs_after","date","nT"), names(smd))
          raw_df <- copy(smd)[, ..raw_keep][
            , `:=`(region = reg, unit = TEMP_MODE, q = q, d_win = win)
          ]
          fwrite(raw_df, file.path(csv_dir_raw, sprintf("SMD_raw_%s.csv", base_tag)))
        }
        ## -------------------------------------------------------------------
        
        # Panel plot
        plots[[k]] <- build_love_plot(
          var_df,
          title_text = sprintf("%s — ±%d time window", reg_label, win),
          xlim_max   = NULL,
          temp_mode  = TEMP_MODE,
          q_tag      = q
        )
        message(sprintf("[OK] %s", attr(smd, "src_file")))
      }
      k <- k + 1L
    }
  }
  
  # Harmonize x-axis across all six subplots for this q
  if (!is.finite(xmax_global) || xmax_global <= 0) xmax_global <- 0.2
  xmax_break <- ceiling(xmax_global * 10) / 10
  breaks_seq <- seq(0, xmax_break, by = 0.1)
  apply_common_scale <- function(p) {
    if (!inherits(p, "gg")) return(p)
    p + scale_x_continuous(limits = c(0, xmax_break),
                           breaks = breaks_seq,
                           minor_breaks = NULL)
  }
  plots <- lapply(plots, apply_common_scale)
  
  # Arrange 2×3 grid (row1: dist; row2: area; cols: 3,7,14)
  panel <- ggarrange(
    plots[[1]], plots[[2]], plots[[3]],
    plots[[4]], plots[[5]], plots[[6]],
    ncol = 3, nrow = 2,
    labels = c("A","B","C","D","E","F"),
    font.label = list(size = 14, face = "bold"),
    common.legend = TRUE, legend = "right",
    align = "hv"
  )
  
  out_png <- file.path(base_dir, sprintf("CBPS_SMD_Panel_%s_%s.png", q, TEMP_MODE))
  ggsave(filename = out_png, plot = panel,
         width = 18, height = 10, units = "in", dpi = 1200, bg = "white")
  message("Saved: ", out_png)
}

## ============================ MASTER CSV EXPORTS ===========================
if (SAVE_CSV_VALUES && length(.all_values)) {
  master_values <- rbindlist(.all_values, fill = TRUE)
  fwrite(master_values, file.path(csv_dir_values,
                                  sprintf("CBPS_SMD_values_master_%s_%s.csv", TEMP_MODE, DATE_AGG)))
}
if (SAVE_CSV_TIDY && length(.all_tidy)) {
  master_tidy <- rbindlist(.all_tidy, fill = TRUE)
  fwrite(master_tidy, file.path(csv_dir_tidy,
                                sprintf("CBPS_SMD_tidy_master_%s_%s.csv", TEMP_MODE, DATE_AGG)))
}

message("CSV(values): ", csv_dir_values)
message("CSV(tidy):   ", csv_dir_tidy)
if (SAVE_RAW_SMD_CSV) message("CSV(raw): ", csv_dir_raw)

