###############################################################################
# Title   : PSM Balance Diagnostics — 2×3 Love-Plot Panels + CSV Exports
# Author  : Juwon Jung
# Version : 1.0
# Date    : 2025-08-19
#
# Purpose :
#   - Load PSM diagnostic bundles saved per combination of:
#       q ∈ {q3,q4,q5}, region ∈ {"dist","area"}, unit ∈ {"rt","tq"},
#       d ∈ {±3, ±7, ±14} (encoded as d3, d7, d14 in filenames).
#   - Summarize |SMD| across dates (median or nT-weighted mean) first at
#     level-level, then collapse to variable-level (max across levels).
#   - Produce 2×3 Love-plot panels per q: rows = {dist, area}, cols = {3,7,14}.
#   - Export per-condition CSVs (wide & tidy) and optional raw SMD tables;
#     also export master CSVs aggregated over all conditions.
#
# Input files (created by your PSM pipeline):
#   base_dir/<region>/psm_UNIT_ATT_boot_{q}_{unit}_{region}_d{win}_WITH_SMD_DIAG.RData
#   Each RData must contain `smd_table` with at least:
#     variable, level, abs_before, abs_after, date, nT, type (categorical/...)
#
# Outputs:
#   Figures:
#     base_dir/PSM_SMD_Panel_{q}_{TEMP_MODE}.png
#   CSV (per-condition):
#     base_dir/SMD_Values/SMD_values_{region}_{unit}_{q}_d{d}.csv    (wide)
#     base_dir/SMD_Values_tidy/SMD_tidy_{region}_{unit}_{q}_d{d}.csv (tidy)
#     base_dir/SMD_Values_raw/SMD_raw_{region}_{unit}_{q}_d{d}.csv   (optional)
#   CSV (master over conditions):
#     base_dir/SMD_Values/PSM_SMD_values_master_{TEMP_MODE}_{DATE_AGG}.csv
#     base_dir/SMD_Values_tidy/PSM_SMD_tidy_master_{TEMP_MODE}_{DATE_AGG}.csv
#
# Notes:
#   - TEMP_MODE affects labeling and output filenames; loading uses that unit.
#   - Only `type == "categorical"` rows are used for panel plots.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(ggpubr)
})

## ============================== USER PATHS =================================
# Change this to your PSM results root directory.
base_dir  <- "C:/path/to/Matching_Results/PSM"

## ============================= GLOBAL SETTINGS =============================
wins      <- c(3, 7, 14)                 # ±k-day windows
qs        <- c("q3", "q4", "q5")         # quantile schemes
DATE_AGG  <- "weighted_mean"             # "median" or "weighted_mean" (weight = nT)
TEMP_MODE <- "rt"                        # "rt" or "tq" (unit in file names & labels)

## =========================== FILE PATH / LOADING ===========================
# PSM file naming pattern (region + unit are encoded in the filename)
psm_path <- function(q, region, win, unit = TEMP_MODE) {
  # psm_UNIT_ATT_boot_{q}_{unit}_{region}_d{win}_WITH_SMD_DIAG.RData
  file.path(base_dir, region,
            sprintf("psm_UNIT_ATT_boot_%s_%s_%s_d%d_WITH_SMD_DIAG.RData",
                    q, unit, region, win))
}

# Load `smd_table` if file exists; otherwise return NULL
load_smd_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  e <- new.env()
  load(path, envir = e)
  if (!exists("smd_table", envir = e)) return(NULL)
  smd <- as.data.table(e$smd_table)
  attr(smd, "src_file") <- path
  smd
}

## =============================== SUMMARIZATION ==============================
# Aggregate SMD across dates → levels → variables.
# - date_agg: "median" or "weighted_mean" (w = nT).
# - Returns a list: var (variable-level summary) and raw_lvl (level-level summary).
summarize_smd <- function(smd, date_agg = c("median","weighted_mean")) {
  date_agg <- match.arg(date_agg)
  
  # Use categorical covariates only
  if ("type" %in% names(smd)) smd <- smd[type == "categorical"]
  
  # Defensive defaults
  if (!("date" %in% names(smd))) smd[, date := NA]
  if (!("nT"   %in% names(smd))) smd[, nT := 1]
  
  # Level-level aggregation across dates
  if (date_agg == "median") {
    lvl <- smd[, .(
      abs_b   = median(abs_before, na.rm = TRUE),
      abs_a   = median(abs_after,  na.rm = TRUE),
      n_dates = uniqueN(date)
    ), by = .(variable, level)]
  } else {
    lvl <- smd[, .(
      abs_b   = { w <- nT; if (!any(is.finite(w))) w[] <- 1; weighted.mean(abs_before, w, na.rm = TRUE) },
      abs_a   = { w <- nT; if (!any(is.finite(w))) w[] <- 1; weighted.mean(abs_after,  w, na.rm = TRUE) },
      n_dates = uniqueN(date)
    ), by = .(variable, level)]
  }
  
  # Clean non-finite; set missing to 0 to avoid dropping variables
  lvl[, c("abs_b","abs_a") := lapply(.SD, function(z) ifelse(is.finite(z), z, NA_real_)),
      .SDcols = c("abs_b","abs_a")]
  lvl[, abs_b := fifelse(is.na(abs_b), 0, abs_b)]
  lvl[, abs_a := fifelse(is.na(abs_a), 0, abs_a)]
  
  # Collapse to variable level (max across levels)
  var <- lvl[, .(
    before   = max(abs_b, na.rm = TRUE),
    after    = max(abs_a, na.rm = TRUE),
    n_levels = .N,
    n_dates  = max(n_dates, na.rm = TRUE)
  ), by = variable]
  
  list(var = var, raw_lvl = lvl)
}

## ============================== LABEL HELPERS ==============================
label_vars <- function(v, temp_mode = c("rt","tq"), q_tag = c("q3","q4","q5")){
  temp_mode <- match.arg(temp_mode); q_tag <- match.arg(q_tag)
  N <- as.integer(sub("^q","", q_tag))
  v <- as.character(v)
  
  repl <- function(x, from, to){ x[x == from] <- to; x }
  
  v <- repl(v, "wind_sp_binary", "Wind speed (binary)")
  v <- repl(v, "rain_binary",    "Precipitation occurrence")
  v <- repl(v, "time_of_day",    "Time of day")
  v <- repl(v, "wind_dir_8",     "Wind directions")
  v <- repl(v, "wday",           "Weekday")
  v <- repl(v, "dist",           "Administrative districts")
  v <- repl(v, "area",           "Residential zones")
  v <- repl(v, "CO_quantile",    sprintf("CO (quantiled = %d)", N))
  v <- repl(v, "O3_quantile",    sprintf("O3 (quantiled = %d)", N))
  v <- repl(v, "temp_match",
            if (temp_mode == "rt") "Temperature (rounded)"
            else sprintf("Temperature (quantiled = %d)", N))
  v
}

## ============================== LOVE PLOT ==================================
# Single Love-plot builder (points only; common 0.1 ticks; shared legend).
build_love_plot <- function(var_df, title_text, xlim_max, q_tag) {
  df_long <- melt(var_df, id.vars = "variable",
                  measure.vars = c("before","after"),
                  variable.name = "Stage", value.name = "SMD")
  df_long[, Stage := factor(Stage, levels = c("before","after"),
                            labels = c("Before","After"))]
  df_long[, var_show := label_vars(variable, TEMP_MODE, q_tag)]
  df_long[, var_show := fct_reorder(var_show, SMD, .fun = function(z) max(z, na.rm = TRUE))]
  
  xmax_break <- ceiling(xlim_max * 10) / 10
  breaks_seq <- seq(0, xmax_break, by = 0.1)
  
  ggplot(df_long, aes(x = SMD, y = var_show)) +
    geom_point(aes(color = Stage), shape = 16, size = 2.6, alpha = 0.95) +
    geom_vline(xintercept = 0.1, linetype = "dashed", linewidth = 0.6, color = "#555555") +
    scale_color_manual(values = c("Before" = "#1b7837", "After" = "#d73027")) +
    scale_x_continuous(limits = c(0, xmax_break), breaks = breaks_seq, minor_breaks = NULL) +
    labs(title = title_text, x = "|SMD|", y = NULL, color = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
}

## =========================== CSV EXPORT OPTIONS ============================
SAVE_CSV_VALUES  <- TRUE     # wide format (per variable: before/after)
SAVE_CSV_TIDY    <- TRUE     # long format (stage = before/after)
SAVE_RAW_SMD_CSV <- FALSE    # optional: raw (variable-level by date/level)

csv_dir_values <- file.path(base_dir, "SMD_Values")
csv_dir_tidy   <- file.path(base_dir, "SMD_Values_tidy")
csv_dir_raw    <- file.path(base_dir, "SMD_Values_raw")
if (SAVE_CSV_VALUES)  dir.create(csv_dir_values, showWarnings = FALSE, recursive = TRUE)
if (SAVE_CSV_TIDY)    dir.create(csv_dir_tidy,   showWarnings = FALSE, recursive = TRUE)
if (SAVE_RAW_SMD_CSV) dir.create(csv_dir_raw,    showWarnings = FALSE, recursive = TRUE)

# Accumulators for master CSVs
.all_values <- list()
.all_tidy   <- list()

## =============================== MAIN PANELS ===============================
for (q in qs) {
  cells <- list()
  
  for (reg in c("dist","area")) {
    reg_title <- if (reg == "dist") "Administrative districts" else "Residential zones"
    
    for (win in wins) {
      smd <- load_smd_if_exists(psm_path(q, reg, win, TEMP_MODE))
      key <- paste(reg, win)
      
      if (is.null(smd)) {
        cells[[key]] <- list(kind = "empty",
                             title = sprintf("%s — ±%d time window (missing)", reg_title, win))
      } else {
        smd_summ <- summarize_smd(smd, DATE_AGG)
        var_df   <- smd_summ$var
        raw_lvl  <- smd_summ$raw_lvl
        
        ## ----------------------------- CSVs --------------------------------
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
          keep <- intersect(c("variable","level","abs_before","abs_after","date","nT","type"), names(smd))
          raw_df <- copy(smd)[, ..keep][
            , `:=`(region = reg, unit = TEMP_MODE, q = q, d_win = win)
          ]
          fwrite(raw_df, file.path(csv_dir_raw, sprintf("SMD_raw_%s.csv", base_tag)))
        }
        ## -------------------------------------------------------------------
        
        cells[[key]] <- list(kind = "plot",
                             var_df = var_df,
                             title  = sprintf("%s — ±%d time window", reg_title, win))
        message(sprintf("[OK] %s", attr(smd, "src_file")))
      }
    }
  }
  
  # Common x-axis upper bound across all six cells for this q
  xmax <- max(vapply(
    cells,
    function(x) if (!identical(x$kind, "plot")) NA_real_
    else max(c(x$var_df$before, x$var_df$after), na.rm = TRUE),
    numeric(1)
  ), na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) xmax <- 0.2
  
  make_cell <- function(reg, win) {
    x <- cells[[paste(reg, win)]]
    if (identical(x$kind, "empty")) {
      ggplot() + theme_void() + labs(title = x$title)
    } else {
      build_love_plot(x$var_df, x$title, xlim_max = xmax, q_tag = q)
    }
  }
  
  panel <- ggarrange(
    make_cell("dist", 3), make_cell("dist", 7), make_cell("dist", 14),
    make_cell("area", 3), make_cell("area", 7), make_cell("area", 14),
    ncol = 3, nrow = 2,
    labels = c("A","B","C","D","E","F"),
    font.label = list(size = 14, face = "bold"),
    common.legend = TRUE, legend = "right",
    align = "hv"
  )
  
  out_png <- file.path(base_dir, sprintf("PSM_SMD_Panel_%s_%s.png", q, TEMP_MODE))
  ggsave(out_png, panel, width = 18, height = 10, units = "in", dpi = 1200, bg = "white")
  message("Saved: ", out_png)
}

## ============================ MASTER CSV EXPORTS ===========================
if (SAVE_CSV_VALUES && length(.all_values)) {
  master_values <- rbindlist(.all_values, fill = TRUE)
  fwrite(master_values,
         file.path(csv_dir_values, sprintf("PSM_SMD_values_master_%s_%s.csv", TEMP_MODE, DATE_AGG)))
}
if (SAVE_CSV_TIDY && length(.all_tidy)) {
  master_tidy <- rbindlist(.all_tidy, fill = TRUE)
  fwrite(master_tidy,
         file.path(csv_dir_tidy, sprintf("PSM_SMD_tidy_master_%s_%s.csv", TEMP_MODE, DATE_AGG)))
}

message("CSV(values): ", csv_dir_values)
message("CSV(tidy):   ", csv_dir_tidy)
if (SAVE_RAW_SMD_CSV) message("CSV(raw): ", csv_dir_raw)

