###############################################################################
# Title   : CEM Balance Diagnostics — Recompute SMDs + 2×3 Panels + CSV Exports
# Author  : (Your Name)
# Version : 1.0
# Date    : 2025-08-19
#
# Purpose :
#   - Iterate over CEM output files (pattern: cem_{region}_{unit}_{q}_d{win}.RData).
#   - Recompute dummy-level SMDs consistently (Before vs After), aggregate to
#     variable-level via max |SMD| across dummies, and save:
#       * condition-wise CSVs (wide, tidy, optional raw),
#       * master CSVs (merged across all conditions),
#       * optional single love plots,
#       * 2×3 panel plots per q (rows={dist, area}, cols={3,7,14}).
#   - Backup original RData before overwriting with updated SMD vectors.
#
# Inputs  :
#   - CEM RData files under `cem_dir` matching ^cem_.*\\.RData$, each containing:
#       final_matched_df (required), matched_df (optional), and possibly previous
#       smd_before / smd_after (to be recomputed in this script).
#   - Source data RData files under `data_dir`:
#       seoul_q3.RData (object `seoul_q3`), seoul_q4.RData (`seoul_q4`), seoul_q5.RData (`seoul_q5`)
#
# Outputs :
#   - Updated RData (same filenames) with recomputed `smd_before`, `smd_after`.
#   - PNG panels: plot_dir/CEM_SMD_Panel_{q}_{TEMP_MODE}.png
#   - CSVs (per condition):
#       csv_dir_values/SMD_values_{region}_{unit}_{q}_d{win}.csv          (wide)
#       csv_dir_tidy/SMD_tidy_{region}_{unit}_{q}_d{win}.csv              (tidy)
#       csv_dir_raw/SMD_raw_{region}_{unit}_{q}_d{win}.csv                (optional)
#   - Master CSVs:
#       csv_dir_values/SMD_values_master_{TEMP_MODE}.csv
#       csv_dir_tidy/SMD_tidy_master_{TEMP_MODE}.csv
#
# Conventions :
#   - “Before” = full sample split by treatment period; “After” = CEM-matched sample.
#   - Variable-level summary = max |SMD| across the variable’s dummy columns.
#   - TEMP_MODE controls axis labels for temperature (“rounded” vs “quantiled = N”).
###############################################################################

suppressPackageStartupMessages({
  library(data.table); library(lubridate); library(ggplot2)
  library(dplyr); library(tidyr); library(ggpubr); library(tools); library(forcats)
})

## ============================== USER PATHS =================================
# Change these three paths for your environment.
cem_dir    <- "C:/path/to/Matching_Results/CEM"
data_dir   <- "C:/path/to/Seoul_RData"           # where seoul_q{3,4,5}.RData live
plot_dir   <- file.path(cem_dir, "SMD_Plots_fixed")
backup_dir <- file.path(cem_dir, "backup_before_smdfix")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE)

## ============================== CSV OPTIONS ================================
SAVE_SINGLE_PNG  <- FALSE  # save per-file single love plots
SAVE_CSV_VALUES  <- TRUE   # save variable-level |SMD| summary (wide)
SAVE_CSV_TIDY    <- TRUE   # save long format (stage = before/after)
SAVE_RAW_SMD_CSV <- FALSE  # save raw dummy-level SMDs (optional)

csv_dir_values <- file.path(cem_dir, "SMD_Values")
csv_dir_tidy   <- file.path(cem_dir, "SMD_Values_tidy")
csv_dir_raw    <- file.path(cem_dir, "SMD_Values_raw")
dir.create(csv_dir_values, showWarnings = FALSE, recursive = TRUE)
dir.create(csv_dir_tidy,   showWarnings = FALSE, recursive = TRUE)
if (SAVE_RAW_SMD_CSV) dir.create(csv_dir_raw, showWarnings = FALSE, recursive = TRUE)

## ============================== GLOBAL SETTINGS ============================
qs        <- c("q3","q4","q5")
wins      <- c(3,7,14)
TEMP_MODE <- "rt"     # "rt" → Temperature (rounded), "tq" → Temperature (quantiled = N)

# Treatment period used to form the Before groups
start_date <- as.Date("2020-01-24")
end_date   <- as.Date("2020-02-09")

# Base covariates (region variable will be appended per-file)
covariates_base <- c("temp_quantile","rain_binary","CO_quantile","O3_quantile",
                     "wind_sp_binary","wind_dir_8","wday","time_of_day")

## ================================ HELPERS ==================================
# Load seoul_q* object from RData (expects object name "seoul_{q}")
load_q <- function(q) {
  f <- file.path(data_dir, sprintf("seoul_%s.RData", q))
  e <- new.env(); load(f, envir = e)
  get(sprintf("seoul_%s", q), envir = e)
}

# SMD (treated − control); if pooled SD is 0 or non-finite, return 0
smd_signed <- function(x1, x0) {
  m1 <- mean(x1, na.rm = TRUE); m0 <- mean(x0, na.rm = TRUE)
  s1 <- sd(x1, na.rm = TRUE);   s0 <- sd(x0, na.rm = TRUE)
  den <- sqrt((s1^2 + s0^2) / 2)
  if (!is.finite(den) || den == 0) return(0)
  (m1 - m0) / den
}

# Pretty variable labels for figures (temperature label depends on TEMP_MODE and q)
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
  v <- repl(v, "dist",           "Districts")
  v <- repl(v, "area",           "Residential zones")
  v <- repl(v, "CO_quantile",    sprintf("CO (quantiled = %d)", N))
  v <- repl(v, "O3_quantile",    sprintf("O3 (quantiled = %d)", N))
  v <- repl(v, "temp_quantile",
            if (temp_mode == "rt") "Temperature (rounded)"
            else sprintf("Temperature (quantiled = %d)", N))
  v
}

pretty_region <- function(region) {
  if (region == "area") "Residential Zones"
  else if (region == "dist") "Administrative Districts"
  else toupper(region)
}

# Collapse dummy-level SMDs to variable-level by max |SMD|
to_var_level <- function(smd_vec, covariates_vars) {
  out <- numeric(length(covariates_vars)); names(out) <- covariates_vars
  nm  <- names(smd_vec)
  for (v in covariates_vars) {
    cols <- grep(paste0("^", v), nm, value = TRUE)
    if (!length(cols)) next
    out[v] <- max(abs(smd_vec[cols]), na.rm = TRUE)
    if (!is.finite(out[v])) out[v] <- NA_real_
  }
  out
}

# CBPS-style Love plot (points only; 0.1 grid; shared legend)
build_love_plot <- function(var_df, title_text, xlim_max, q_tag) {
  df_long <- melt(var_df, id.vars="variable",
                  measure.vars=c("before","after"),
                  variable.name="Stage", value.name="SMD")
  df_long[, Stage := factor(Stage, levels=c("before","after"),
                            labels=c("Before","After"))]
  df_long[, var_show := label_vars(variable, TEMP_MODE, q_tag)]
  df_long[, var_show := fct_reorder(var_show, SMD, .fun=function(z) max(z, na.rm=TRUE))]
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

## ===================== STORAGE FOR PANEL ASSEMBLY ==========================
# key format: sprintf("%s_%s_%s_d%d", region, unit, q, d)
var_store <- new.env(parent = emptyenv())

## =================== MAIN LOOP: RECOMPUTE & SAVE SMDs ======================
files <- list.files(cem_dir, pattern = "^cem_.*\\.RData$", full.names = TRUE)

for (f in files) {
  tag <- file_path_sans_ext(basename(f))  # e.g., cem_area_rt_q5_d3
  message("Processing: ", tag)
  
  toks <- strsplit(tag, "_", fixed = TRUE)[[1]]
  if (length(toks) < 5) next
  region <- toks[2]    # "area" / "dist" (region variable name for this file)
  unit   <- toks[3]    # "rt" / "tq"
  q      <- toks[4]    # "q3" / "q4" / "q5"
  d_win  <- as.integer(sub("^d", "", toks[5]))
  
  # Covariates used for this file = base + region
  covariates_this <- unique(c(covariates_base, region))
  
  # Load RData & check contents
  env <- new.env(); load(f, envir = env)
  if (!exists("final_matched_df", env)) { warning("final_matched_df not found: ", tag); next }
  final_matched_df <- as.data.table(env$final_matched_df)
  
  # BEFORE: from source `seoul_q*` with period split
  seoul <- as.data.table(load_q(q))
  covars_in_seoul <- intersect(covariates_this, names(seoul))
  if (!length(covars_in_seoul)) { warning("No covariates in seoul data: ", tag); next }
  
  unmatched_treated <- seoul[date >= start_date & date <= end_date, ..covars_in_seoul]
  unmatched_treated[, treated := 1L]
  unmatched_control <- seoul[!(date >= start_date & date <= end_date), ..covars_in_seoul]
  unmatched_control[, treated := 0L]
  before_df <- rbindlist(list(unmatched_treated, unmatched_control), use.names = TRUE)
  before_df[, stage := "Before"]
  
  # AFTER: from CEM result (use intersection of available covariates)
  covars_in_after <- intersect(covariates_this, names(final_matched_df))
  if (!length(covars_in_after)) {
    warning("No overlapping covariates in final_matched_df: ", tag)
    next
  }
  after_df <- copy(final_matched_df)[, .SD, .SDcols = c(covars_in_after, "treated")]
  after_df[, stage := "After"]
  
  # Common covariates for both stages
  covariates_use <- Reduce(intersect, list(covariates_this, names(before_df), names(after_df)))
  if (!length(covariates_use)) { warning("No usable covariates: ", tag); next }
  
  # Build a unified design matrix of dummy variables
  both <- rbindlist(list(before_df, after_df), use.names = TRUE, fill = TRUE)
  both <- both[complete.cases(both[, ..covariates_use]) & !is.na(treated) & !is.na(stage)]
  for (v in covariates_use) both[[v]] <- as.factor(both[[v]])
  mm_formula_use <- as.formula(paste("~ 0 +", paste(sprintf("%s", covariates_use), collapse = " + ")))
  mm <- model.matrix(mm_formula_use, data = both)
  
  # Index masks
  idxB1 <- both$stage == "Before" & both$treated == 1L
  idxB0 <- both$stage == "Before" & both$treated == 0L
  idxA1 <- both$stage == "After"  & both$treated == 1L
  idxA0 <- both$stage == "After"  & both$treated == 0L
  
  # Recompute dummy-level SMDs
  smd_before <- apply(mm, 2, function(col) smd_signed(col[idxB1], col[idxB0]))
  smd_after  <- apply(mm, 2, function(col) smd_signed(col[idxA1], col[idxA0]))
  
  # Backup then overwrite original RData with updated SMDs (keep original matched_df)
  file.copy(f, file.path(backup_dir, basename(f)), overwrite = FALSE)
  matched_df <- if (exists("matched_df", env)) env$matched_df else NULL
  save(matched_df, final_matched_df, smd_before, smd_after, file = f)
  
  # Variable-level summary (max |SMD| over dummy columns); include region variable
  vb <- to_var_level(smd_before, covariates_use)
  va <- to_var_level(smd_after,  covariates_use)
  var_df <- data.table(variable = names(vb), before = as.numeric(vb), after = as.numeric(va))
  var_df <- var_df[is.finite(before) | is.finite(after)]
  
  # Store for panel assembly
  assign(sprintf("%s_%s_%s_d%d", region, unit, q, d_win), var_df, envir = var_store)
  
  # ---------- Condition-wise CSV exports ----------
  base_tag <- sub("^cem_", "", tag)  # e.g., area_rt_q5_d3
  if (SAVE_CSV_VALUES) {
    out_values <- copy(var_df)[, `:=`(region = region, unit = unit, q = q, d_win = d_win)]
    fwrite(out_values, file.path(csv_dir_values, sprintf("SMD_values_%s.csv", base_tag)))
  }
  if (SAVE_CSV_TIDY) {
    tidy <- melt(var_df, id.vars = "variable",
                 measure.vars = c("before","after"),
                 variable.name = "stage", value.name = "SMD")
    tidy[, `:=`(region = region, unit = unit, q = q, d_win = d_win)]
    fwrite(tidy, file.path(csv_dir_tidy, sprintf("SMD_tidy_%s.csv", base_tag)))
  }
  if (SAVE_RAW_SMD_CSV) {
    raw_df <- data.table(name = colnames(mm),
                         before = as.numeric(smd_before[colnames(mm)]),
                         after  = as.numeric(smd_after[colnames(mm)]))
    raw_df[, `:=`(region = region, unit = unit, q = q, d_win = d_win)]
    fwrite(raw_df, file.path(csv_dir_raw, sprintf("SMD_raw_%s.csv", base_tag)))
  }
  
  # Optional: single-plot export
  if (SAVE_SINGLE_PNG) {
    xmax_single <- max(c(var_df$before, var_df$after), na.rm = TRUE)
    if (!is.finite(xmax_single)) xmax_single <- 0.2
    p_single <- build_love_plot(
      var_df,
      title_text = sprintf("%s — ±%d Time Window", pretty_region(region), d_win),
      xlim_max = xmax_single, q_tag = q
    )
    ggsave(filename = file.path(plot_dir, paste0("SMD_unified_", tag, ".png")),
           plot = p_single, width = 8, height = 10, dpi = 600, bg = "white")
  }
}

## =========================== MASTER CSV EXPORTS ============================
keys <- ls(envir = var_store)
if (length(keys)) {
  all_values <- rbindlist(lapply(keys, function(k) {
    v <- get(k, envir = var_store)
    tt <- strsplit(k, "_", fixed = TRUE)[[1]]  # c(region, unit, q, dX)
    data.table(
      variable = v$variable, before = v$before, after = v$after,
      region = tt[1], unit = tt[2], q = tt[3], d_win = as.integer(sub("^d","", tt[4]))
    )
  }), fill = TRUE)
  
  if (SAVE_CSV_VALUES) {
    fwrite(all_values, file.path(csv_dir_values, sprintf("SMD_values_master_%s.csv", TEMP_MODE)))
  }
  if (SAVE_CSV_TIDY) {
    all_tidy <- melt(all_values,
                     id.vars = c("variable","region","unit","q","d_win"),
                     measure.vars = c("before","after"),
                     variable.name = "stage", value.name = "SMD")
    fwrite(all_tidy, file.path(csv_dir_tidy, sprintf("SMD_tidy_master_%s.csv", TEMP_MODE)))
  }
}

## ===================== 2×3 PANELS PER q (dist/area × 3/7/14) =================
for (q in qs) {
  get_cell <- function(region, d) {
    key <- sprintf("%s_%s_%s_d%d", region, TEMP_MODE, q, d)
    if (exists(key, envir = var_store, inherits = FALSE)) get(key, envir = var_store) else NULL
  }
  
  cells <- list(
    list(reg="dist", win=3,  var_df=get_cell("dist",3)),
    list(reg="dist", win=7,  var_df=get_cell("dist",7)),
    list(reg="dist", win=14, var_df=get_cell("dist",14)),
    list(reg="area", win=3,  var_df=get_cell("area",3)),
    list(reg="area", win=7,  var_df=get_cell("area",7)),
    list(reg="area", win=14, var_df=get_cell("area",14))
  )
  
  # Shared x-axis bound across the six plots for a given q
  xmax <- max(vapply(
    cells,
    function(x) {
      if (is.null(x$var_df)) return(NA_real_)
      vals <- c(x$var_df$before, x$var_df$after)
      vals <- vals[is.finite(vals)]
      if (!length(vals)) NA_real_ else max(vals)
    },
    numeric(1)
  ), na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) xmax <- 0.2
  
  mk <- function(cell){
    reg_title <- if (cell$reg=="dist") "Administrative Districts" else "Residential Zones"
    if (is.null(cell$var_df)) {
      ggplot()+theme_void()+labs(title=sprintf("%s — ±%d Time Window (missing)", reg_title, cell$win))
    } else {
      build_love_plot(cell$var_df,
                      title_text = sprintf("%s — ±%d Time Window", reg_title, cell$win),
                      xlim_max = xmax, q_tag = q)
    }
  }
  plots <- lapply(cells, mk)
  
  panel <- ggarrange(
    plots[[1]], plots[[2]], plots[[3]],
    plots[[4]], plots[[5]], plots[[6]],
    ncol = 3, nrow = 2,
    labels = c("A","B","C","D","E","F"),
    font.label = list(size = 14, face = "bold"),
    common.legend = TRUE, legend = "right",
    align = "hv"
  )
  
  out_png <- file.path(plot_dir, sprintf("CEM_SMD_Panel_%s_%s.png", q, TEMP_MODE))
  ggsave(out_png, panel, width = 18, height = 10, units = "in", dpi = 1200, bg = "white")
  message("Saved panel: ", out_png)
}

message("Single-plot directory: ", plot_dir)
message("Backup directory: ", backup_dir)
message("CSV(values) directory: ", csv_dir_values)
message("CSV(tidy) directory:   ", csv_dir_tidy)
if (SAVE_RAW_SMD_CSV) message("CSV(raw) directory: ", csv_dir_raw)
