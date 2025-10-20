###############################################################################
# Title   : Combined SMD Panels (3×3) — CEM vs PSM vs CBPS by Window
# Author  : Juwon Jung
# Version : 1.0
# Date    : 2025-08-19
#
# Purpose :
#   - For each quantile scheme q ∈ {q3, q4, q5}, assemble a 3×3 panel of
#     |SMD| “love plots”, comparing CEM/PSM/CBPS (rows) across windows
#     {±3, ±7, ±14} (columns), for both region types {dist, area} and units
#     {rt, tq}.
#   - Inputs are the *master CSVs* previously exported by each pipeline:
#       CEM : SMD_values_master_{unit}.csv
#       PSM : PSM_SMD_values_master_{unit}_{DATE_AGG}.csv
#       CBPS: CBPS_SMD_values_master_{unit}_{DATE_AGG}.csv
#
# Inputs (directory layout):
#   root/
#     CEM/   SMD_values_master_{unit}.csv
#     PSM/   PSM_SMD_values_master_{unit}_{DATE_AGG}.csv
#     CBPS/  CBPS_SMD_values_master_{unit}_{DATE_AGG}.csv
#
# Expected CSV schema:
#   columns: variable, before, after, region, unit, q, d_win
#
# Outputs:
#   root/Combined_Panels_3x3_by_q/SMD_3x3_{region}_{unit}_{q}.png
#
# Notes:
#   - “before/after” are *variable-level* |SMD| summaries (max over dummy levels),
#     as produced by the earlier scripts.
#   - Temperature labels follow the chosen unit: rt = “rounded”, tq = “quantiled = N”.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(ggpubr)
  library(ggtext)   
  library(grid)     
})

## ============================== USER PATHS =================================
root      <- "Your Path"
cem_dir   <- file.path(root, "CEM")
psm_dir   <- file.path(root, "PSM")
cbps_dir  <- file.path(root, "CBPS")
out_dir   <- file.path(root, "Combined_Panels_3x3_by_q")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ========================== MASTER CSV FILENAMES ===========================
cem_master  <- function(unit)
  file.path(cem_dir,  sprintf("SMD_values_master_%s.csv", unit))

psm_master  <- function(unit, date_agg = c("weighted_mean","median")) {
  date_agg <- match.arg(date_agg)
  file.path(psm_dir,  sprintf("PSM_SMD_values_master_%s_%s.csv",  unit, date_agg))
}
cbps_master <- function(unit, date_agg = c("weighted_mean","median")) {
  date_agg <- match.arg(date_agg)
  file.path(cbps_dir, sprintf("CBPS_SMD_values_master_%s_%s.csv", unit, date_agg))
}

## ================================ PARAMETERS ===============================
wins      <- c(3, 7, 14)
methods   <- c("CEM","PSM","CBPS")
units     <- c("rt","tq")
regions   <- c("dist","area")
qs        <- c("q3","q4","q5")
DATE_AGG  <- "weighted_mean"

## ================================ STYLING ==================================
TEXT_Y_SIZE <- 14
DOT_SIZE    <- 3.6
LEGEND_TEXT_SIZE   <- 20
LEGEND_DOT_SIZE    <- 4.5
LEGEND_KEY_SIZE_PT <- 24
LABEL_APPEND_WINDOW <- FALSE

METHOD_TITLE_PT <- 18  
make_panel_title <- function(method, w){
  sprintf(
    "<span style='font-size:%dpt; font-weight:700'>%s</span> — ±%d-day Time Window",
    METHOD_TITLE_PT, method, w
  )
}

## =============================== LABEL HELPERS =============================
label_vars_agg <- function(v, temp_mode = c("rt","tq"), q_tag = NULL, d_win = NULL){
  temp_mode <- match.arg(temp_mode)
  q_num <- if (!is.null(q_tag)) as.integer(sub("^q","", q_tag)) else NA_integer_
  qsuf  <- if (is.finite(q_num)) sprintf(" (quantiled = %d)", q_num) else ""
  wsuf  <- if (!is.null(d_win) && LABEL_APPEND_WINDOW) sprintf("; ±%d", d_win) else ""
  suf_qw <- paste0(qsuf, wsuf)
  
  v <- as.character(v)
  v[v %in% c("temp_match","temp_quantile","temp")] <-
    if (temp_mode == "rt") "Temperature (rounded)" else paste0("Temperature", suf_qw)
  v[v %in% c("CO_quantile","CO")] <- paste0("CO",  suf_qw)
  v[v %in% c("O3_quantile","O3")] <- paste0("O3",  suf_qw)
  v[v == "wind_sp_binary"] <- "Wind speed (binary)"
  v[v == "rain_binary"]    <- "Precipitation occurrence"
  v[v == "time_of_day"]    <- "Time of day"
  v[v == "wind_dir_8"]     <- "Wind directions"
  v[v == "wday"]           <- "Weekday"
  v[v == "dist"]           <- "Districts"
  v[v == "area"]           <- "Residential zones"
  v
}

pretty_region <- function(r) if (r == "dist") "Administrative Districts" else "Residential Zones"

## ============================== MASTER LOADER ==============================
load_master <- function(method, unit){
  f <- switch(method,
              "CEM"  = cem_master(unit),
              "PSM"  = psm_master(unit, DATE_AGG),
              "CBPS" = cbps_master(unit, DATE_AGG))
  if (!file.exists(f)) stop("Master CSV not found: ", f)
  dt <- fread(f)
  need <- c("variable","before","after","region","unit","q","d_win")
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop(sprintf("[%s] missing columns: %s", method, paste(miss, collapse=", ")))
  dt[, method := method]
  dt
}

## ============================= CELL DATA BUILDER ===========================
make_cell_var_df <- function(master_dt, region_val, unit_val, win_val, q_tag){
  sub <- master_dt[region == region_val & unit == unit_val & d_win == win_val & q == q_tag]
  if (nrow(sub) == 0) return(NULL)
  sub[, .(before = max(before, na.rm = TRUE),
          after  = max(after,  na.rm = TRUE)), by = .(variable)]
}

## ================================ LOVE PLOT ================================
build_love_plot <- function(var_df, title_text, temp_mode, q_tag, xlim_max, d_win = NULL){
  if (is.null(var_df) || !nrow(var_df)) {
    return(
      ggplot() + theme_void() + labs(title = paste0(title_text, " (missing)")) +
        theme(plot.title = ggtext::element_markdown())
    )
  }
  df_long <- melt(copy(var_df), id.vars = "variable",
                  measure.vars = c("before","after"),
                  variable.name = "Stage", value.name = "SMD")
  df_long[, Stage := factor(Stage, levels = c("before","after"),
                            labels = c("Before","After"))]
  df_long[, var_show := label_vars_agg(variable, temp_mode = temp_mode, q_tag = q_tag, d_win = d_win)]
  df_long[, var_show := forcats::fct_reorder(var_show, SMD, .fun = function(z) max(z, na.rm = TRUE))]
  
  xmax_break <- ceiling(xlim_max * 10) / 10
  breaks_seq <- seq(0, xmax_break, by = 0.1)
  
  ggplot(df_long, aes(x = SMD, y = var_show)) +
    geom_point(aes(shape = Stage), size = DOT_SIZE, color = "black", alpha = 0.95) +
    scale_shape_manual(
      values = c("Before" = 1, "After" = 16),
      guide = guide_legend(override.aes = list(size = LEGEND_DOT_SIZE, color = "black"))
    ) +
    geom_vline(xintercept = 0.1, linetype = "dashed", linewidth = 0.6, color = "#555555") +
    scale_x_continuous(limits = c(0, xmax_break), breaks = breaks_seq, minor_breaks = NULL) +
    labs(title = title_text, x = "|SMD|", y = NULL, shape = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position   = "right",
      legend.text       = element_text(size = LEGEND_TEXT_SIZE),
      legend.title      = element_text(size = LEGEND_TEXT_SIZE),
      legend.key.height = grid::unit(LEGEND_KEY_SIZE_PT, "pt"),
      legend.key.width  = grid::unit(LEGEND_KEY_SIZE_PT, "pt"),
      panel.grid.minor  = element_blank(),
      panel.grid.major.y= element_blank(),
      axis.text.y       = element_text(size = TEXT_Y_SIZE, face = "bold"),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.title        = ggtext::element_markdown()
    )
}

## ================================ MAIN LOOP ================================
for (reg in regions) {
  for (unit in units) {
    cem  <- load_master("CEM",  unit)
    psm  <- load_master("PSM",  unit)
    cbps <- load_master("CBPS", unit)
    master_all <- rbindlist(list(cem, psm, cbps), use.names = TRUE, fill = TRUE)
    
    for (q_tag in qs) {
      cell_var <- list(); xmax <- 0
      
      for (m in methods) {
        for (w in wins) {
          key <- paste(m, w, sep = "_")
          var_df <- make_cell_var_df(master_all[method == m],
                                     region_val = reg, unit_val = unit, win_val = w, q_tag = q_tag)
          cell_var[[key]] <- var_df
          if (!is.null(var_df) && nrow(var_df)) {
            xmax <- max(xmax, var_df[, max(c(before, after), na.rm = TRUE)], na.rm = TRUE)
          }
        }
      }
      if (!is.finite(xmax) || xmax <= 0) xmax <- 0.2
      
      mk <- function(m, w){
        ttl <- make_panel_title(m, w)  
        build_love_plot(cell_var[[paste(m, w, sep = "_")]],
                        ttl, temp_mode = unit, q_tag = q_tag, xlim_max = xmax, d_win = w)
      }
      
      p11 <- mk("CEM",3);  p12 <- mk("CEM",7);  p13 <- mk("CEM",14)
      p21 <- mk("PSM",3);  p22 <- mk("PSM",7);  p23 <- mk("PSM",14)
      p31 <- mk("CBPS",3); p32 <- mk("CBPS",7); p33 <- mk("CBPS",14)
      
      panel <- ggarrange(
        p11, p12, p13,
        p21, p22, p23,
        p31, p32, p33,
        ncol = 3, nrow = 3,
        font.label   = list(size = 14, face = "bold"),
        common.legend = TRUE, legend = "right",
        align = "hv"
      )
      
      out_png <- file.path(out_dir, sprintf("SMD_3x3_%s_%s_%s.png", reg, unit, q_tag))
      ggsave(out_png, panel, width = 20, height = 14, units = "in", dpi = 600, bg = "white")
      message("Saved: ", out_png,
              "  (", if (reg == "dist") "Administrative Districts" else "Residential Zones",
              ", unit=", unit, ", q=", q_tag, ")")
    }
  }
}

