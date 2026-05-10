###############################################################################
# Updated SMD 3x3 Panel Generator
# rows    = CEM, PSM, CBPS
# columns = ±3, ±7, ±14 day windows
#
# Input folders:
#   CEM_results  : all* SMD csv
#   PSM_results  : all* SMD csv
#   CBPS_results : all* SMD csv
#
# Output:
#   SMD_3x3_dist_rt_q3.png / .pdf
#   SMD_3x3_dist_rt_q4.png / .pdf
#   SMD_3x3_dist_rt_q5.png / .pdf
#   SMD_3x3_area_rt_q3.png / .pdf
#   SMD_3x3_area_rt_q4.png / .pdf
#   SMD_3x3_area_rt_q5.png / .pdf
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(ggpubr)
  library(ggtext)
  library(grid)
})

## =============================================================================
## 0. USER PATHS
## =============================================================================

psm_dir  <- "dir/PSM_results"
cem_dir  <- "dir/CEM_results"
cbps_dir <- "dir/CBPS_results"

out_dir <- "dir/out"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## 1. PARAMETERS
## =============================================================================

wins    <- c(3, 7, 14)
methods <- c("CEM", "PSM", "CBPS")

regions_to_plot <- c("dist", "area")
q_to_plot <- c("q3", "q4", "q5")

ONLY_RT <- TRUE

MIN_XMAX <- 0.4

X_BREAK_BY <- 0.1

## =============================================================================
## 2. PLOT STYLE
## =============================================================================

TEXT_Y_SIZE        <- 13
TEXT_X_SIZE        <- 10
AXIS_TITLE_SIZE    <- 11
DOT_SIZE           <- 3.4
LEGEND_TEXT_SIZE   <- 18
LEGEND_DOT_SIZE    <- 4.5
LEGEND_KEY_SIZE_PT <- 22
METHOD_TITLE_PT    <- 17

LABEL_APPEND_WINDOW <- FALSE

make_panel_title <- function(method, w) {
  sprintf(
    "<span style='font-size:%dpt; font-weight:700'>%s</span> — ±%d-day Time Window",
    METHOD_TITLE_PT, method, w
  )
}

## =============================================================================
## 3. FILE FINDER
## =============================================================================

find_smd_file <- function(method_dir, method) {
  
  if (!dir.exists(method_dir)) {
    stop("Directory does not exist: ", method_dir)
  }
  
  csv_files <- list.files(
    method_dir,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(csv_files) == 0) {
    stop("No csv files found in: ", method_dir)
  }
  
  bn <- basename(csv_files)
  
  cand_plot <- csv_files[
    grepl("all", bn, ignore.case = TRUE) &
      grepl("smd", bn, ignore.case = TRUE) &
      grepl("plot", bn, ignore.case = TRUE) &
      grepl("summary", bn, ignore.case = TRUE)
  ]
  
  cand_var <- csv_files[
    grepl("all", bn, ignore.case = TRUE) &
      grepl("smd", bn, ignore.case = TRUE) &
      grepl("variable", bn, ignore.case = TRUE) &
      grepl("summary", bn, ignore.case = TRUE)
  ]
  
  cand_master <- csv_files[
    grepl("SMD_values_master|PSM_SMD_values_master|CBPS_SMD_values_master",
          bn, ignore.case = TRUE)
  ]
  
  if (length(cand_plot) > 0) {
    message(sprintf("[%s] using all_smd_plot_summary: %s",
                    method, paste(basename(cand_plot), collapse = ", ")))
    return(cand_plot)
  }
  
  if (length(cand_var) > 0) {
    message(sprintf("[%s] using all_smd_variable_summary: %s",
                    method, paste(basename(cand_var), collapse = ", ")))
    return(cand_var)
  }
  
  if (length(cand_master) > 0) {
    wm <- cand_master[grepl("weighted_mean", basename(cand_master), ignore.case = TRUE)]
    if (length(wm) > 0) cand_master <- wm
    
    message(sprintf("[%s] using old master SMD file: %s",
                    method, paste(basename(cand_master), collapse = ", ")))
    return(cand_master)
  }
  
  stop(
    sprintf(
      "[%s] Cannot find SMD summary file in %s. Expected all_smd_plot_summary, all_smd_variable_summary, or SMD_values_master.",
      method, method_dir
    )
  )
}

## =============================================================================
## 4. NORMALIZE SMD DATA
##    target schema:
##    method, variable, before, after, region, unit, q, d_win
## =============================================================================

standardize_names <- function(dt) {
  setDT(dt)
  
  old_names <- names(dt)
  new_names <- old_names
  
  new_names <- gsub("\\s+", "_", new_names)
  new_names <- gsub("\\.", "_", new_names)
  new_names <- gsub("-", "_", new_names)
  
  setnames(dt, old_names, new_names)
  dt
}

normalize_one_smd_file <- function(file, method) {
  
  dt <- fread(file)
  dt <- standardize_names(dt)
  
  if (all(c("variable", "before", "after", "region", "unit", "q", "d_win") %in% names(dt))) {
    
    out <- copy(dt)
    
    if (!("method" %in% names(out))) out[, method := method]
    
    out[, method := toupper(as.character(method))]
    out[, variable := as.character(variable)]
    out[, region := as.character(region)]
    out[, unit := as.character(unit)]
    out[, q := as.character(q)]
    out[!grepl("^q", q), q := paste0("q", q)]
    out[, d_win := as.integer(d_win)]
    out[, before := abs(as.numeric(before))]
    out[, after  := abs(as.numeric(after))]
    
    return(out[, .(method, variable, before, after, region, unit, q, d_win)])
  }
  
  before_candidates <- c("before_smd", "smd_before", "pre_smd", "unmatched_smd")
  after_candidates  <- c("after_smd", "smd_after", "post_smd", "matched_smd")
  
  before_col <- intersect(before_candidates, names(dt))
  after_col  <- intersect(after_candidates, names(dt))
  
  if (length(before_col) > 0 && length(after_col) > 0) {
    
    if (!("method" %in% names(dt))) dt[, method := method]
    
    if (!("variable" %in% names(dt))) {
      if ("covariate" %in% names(dt)) setnames(dt, "covariate", "variable")
    }
    
    if (!("region" %in% names(dt))) {
      if ("region_var" %in% names(dt)) setnames(dt, "region_var", "region")
    }
    
    if (!("d_win" %in% names(dt))) {
      if ("day_window" %in% names(dt)) setnames(dt, "day_window", "d_win")
    }
    
    if (!("unit" %in% names(dt))) {
      if ("temp_var" %in% names(dt)) {
        dt[, unit := fifelse(
          tolower(as.character(temp_var)) %in% c("temp_quantile", "tq", "temp_q"),
          "tq",
          "rt"
        )]
      } else {
        dt[, unit := "rt"]
      }
    }
    
    if (!("q" %in% names(dt))) stop("[", basename(file), "] missing q column.")
    if (!("variable" %in% names(dt))) stop("[", basename(file), "] missing variable/covariate column.")
    if (!("region" %in% names(dt))) stop("[", basename(file), "] missing region/region_var column.")
    if (!("d_win" %in% names(dt))) stop("[", basename(file), "] missing d_win/day_window column.")
    
    setnames(dt, before_col[1], "before")
    setnames(dt, after_col[1],  "after")
    
    out <- dt[
      ,
      .(
        before = max(abs(as.numeric(before)), na.rm = TRUE),
        after  = max(abs(as.numeric(after)),  na.rm = TRUE)
      ),
      by = .(method, variable, region, unit, q, d_win)
    ]
    
    out[is.infinite(before), before := NA_real_]
    out[is.infinite(after),  after  := NA_real_]
    
    out[, method := toupper(as.character(method))]
    out[, q := as.character(q)]
    out[!grepl("^q", q), q := paste0("q", q)]
    out[, d_win := as.integer(d_win)]
    
    return(out[, .(method, variable, before, after, region, unit, q, d_win)])
  }
  
  if (!("method" %in% names(dt))) dt[, method := method]
  
  if (!("variable" %in% names(dt))) {
    if ("covariate" %in% names(dt)) {
      setnames(dt, "covariate", "variable")
    } else {
      stop("[", basename(file), "] missing variable/covariate column.")
    }
  }
  
  if (!("stage" %in% names(dt))) {
    if ("sample" %in% names(dt)) {
      setnames(dt, "sample", "stage")
    } else if ("status" %in% names(dt)) {
      setnames(dt, "status", "stage")
    } else {
      stop("[", basename(file), "] missing stage/sample/status column.")
    }
  }
  
  if (!("SMD" %in% names(dt))) {
    if ("max_abs_smd" %in% names(dt)) {
      setnames(dt, "max_abs_smd", "SMD")
    } else if ("smd" %in% names(dt)) {
      setnames(dt, "smd", "SMD")
    } else if ("abs_smd" %in% names(dt)) {
      setnames(dt, "abs_smd", "SMD")
    } else {
      stop("[", basename(file), "] missing SMD/max_abs_smd/smd/abs_smd column.")
    }
  }
  
  if (!("region" %in% names(dt))) {
    if ("region_var" %in% names(dt)) {
      setnames(dt, "region_var", "region")
    } else {
      stop("[", basename(file), "] missing region/region_var column.")
    }
  }
  
  if (!("d_win" %in% names(dt))) {
    if ("day_window" %in% names(dt)) {
      setnames(dt, "day_window", "d_win")
    } else {
      stop("[", basename(file), "] missing d_win/day_window column.")
    }
  }
  
  if (!("q" %in% names(dt))) {
    stop("[", basename(file), "] missing q column.")
  }
  
  if (!("unit" %in% names(dt))) {
    if ("temp_var" %in% names(dt)) {
      dt[, unit := fifelse(
        tolower(as.character(temp_var)) %in% c("temp_quantile", "tq", "temp_q"),
        "tq",
        "rt"
      )]
    } else {
      dt[, unit := "rt"]
    }
  }
  
  dt[, method := toupper(as.character(method))]
  dt[, variable := as.character(variable)]
  dt[, region := as.character(region)]
  dt[, unit := as.character(unit)]
  dt[, q := as.character(q)]
  dt[!grepl("^q", q), q := paste0("q", q)]
  dt[, d_win := as.integer(d_win)]
  
  dt[, stage := tolower(as.character(stage))]
  dt[, stage := gsub("\\s+", "_", stage)]
  
  dt[stage %in% c("after", "matched", "after_matching", "post", "weighted"), stage := "after"]
  dt[stage %in% c("before", "unmatched", "before_matching", "pre", "unweighted"), stage := "before"]
  
  dt[, SMD := abs(as.numeric(SMD))]
  
  long <- dt[
    stage %in% c("before", "after"),
    .(
      SMD = max(SMD, na.rm = TRUE)
    ),
    by = .(method, variable, region, unit, q, d_win, stage)
  ]
  
  long[is.infinite(SMD), SMD := NA_real_]
  
  wide <- dcast(
    long,
    method + variable + region + unit + q + d_win ~ stage,
    value.var = "SMD",
    fun.aggregate = function(x) {
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      max(x)
    }
  )
  
  if (!("before" %in% names(wide))) wide[, before := NA_real_]
  if (!("after"  %in% names(wide))) wide[, after  := NA_real_]
  
  wide[, .(method, variable, before, after, region, unit, q, d_win)]
}

load_method_smd <- function(method, method_dir) {
  
  files <- find_smd_file(method_dir, method)
  
  out <- rbindlist(
    lapply(files, normalize_one_smd_file, method = method),
    use.names = TRUE,
    fill = TRUE
  )
  
  out[, method := factor(toupper(as.character(method)), levels = methods)]
  out[]
}

## =============================================================================
## 5. VARIABLE LABELS
## =============================================================================

label_vars_agg <- function(v, temp_mode = c("rt", "tq"), q_tag = NULL, d_win = NULL) {
  
  temp_mode <- match.arg(temp_mode)
  
  q_num <- if (!is.null(q_tag)) as.integer(sub("^q", "", q_tag)) else NA_integer_
  qsuf <- if (is.finite(q_num)) sprintf(" (quantiled = %d)", q_num) else ""
  wsuf <- if (!is.null(d_win) && LABEL_APPEND_WINDOW) sprintf("; ±%d", d_win) else ""
  suf_qw <- paste0(qsuf, wsuf)
  
  v <- as.character(v)
  
  v[v %in% c("temp_match", "temp_quantile", "temp")] <-
    if (temp_mode == "rt") "Temperature (rounded)" else paste0("Temperature", suf_qw)
  
  v[v %in% c("temp_round", "round_temp", "temp_rounded")] <- "Temperature (rounded)"
  
  v[v %in% c("CO_quantile", "CO", "co_quantile", "co")] <- paste0("CO", suf_qw)
  v[v %in% c("O3_quantile", "O3", "o3_quantile", "o3")] <- paste0("O3", suf_qw)
  
  v[v %in% c("wind_sp_binary", "wind_sp_ind", "wind_speed_binary")] <- "Wind speed (binary)"
  v[v %in% c("rain_binary", "rain_indicator", "precipitation")] <- "Precipitation occurrence"
  v[v %in% c("time_of_day", "hour_group")] <- "Time of day"
  v[v %in% c("wind_dir_8", "wind_direction", "wind_dir")] <- "Wind directions"
  v[v %in% c("wday", "weekday")] <- "Weekday"
  v[v %in% c("dist", "district")] <- "Districts"
  v[v %in% c("area", "region5", "residential_zone")] <- "Residential zones"
  
  v
}

pretty_region <- function(r) {
  if (r == "dist") "Administrative Districts" else "Residential Zones"
}

## =============================================================================
## 6. CELL DATA BUILDER
## =============================================================================

make_cell_var_df <- function(master_dt, method_val, region_val, unit_val, win_val, q_tag) {
  
  sub <- master_dt[
    as.character(method) == method_val &
      region == region_val &
      unit == unit_val &
      d_win == win_val &
      q == q_tag
  ]
  
  if (nrow(sub) == 0) return(NULL)
  
  out <- sub[
    ,
    .(
      before = max(before, na.rm = TRUE),
      after  = max(after,  na.rm = TRUE)
    ),
    by = .(variable)
  ]
  
  out[is.infinite(before), before := NA_real_]
  out[is.infinite(after),  after  := NA_real_]
  
  out[]
}

## =============================================================================
## 7. LOVE PLOT BUILDER
## =============================================================================

build_love_plot <- function(var_df, title_text, temp_mode, q_tag, xlim_max, d_win = NULL) {
  
  if (is.null(var_df) || !nrow(var_df)) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste0(title_text, " (missing)")) +
        theme(plot.title = ggtext::element_markdown())
    )
  }
  
  df_long <- melt(
    copy(var_df),
    id.vars = "variable",
    measure.vars = c("before", "after"),
    variable.name = "Stage",
    value.name = "SMD"
  )
  
  df_long <- df_long[is.finite(SMD)]
  
  df_long[, Stage := factor(
    Stage,
    levels = c("before", "after"),
    labels = c("Before", "After")
  )]
  
  df_long[, var_show := label_vars_agg(
    variable,
    temp_mode = temp_mode,
    q_tag = q_tag,
    d_win = d_win
  )]
  
  df_long[, var_show := forcats::fct_reorder(
    var_show,
    SMD,
    .fun = function(z) max(z, na.rm = TRUE)
  )]
  
  xmax_break <- ceiling(max(MIN_XMAX, xlim_max, na.rm = TRUE) * 10) / 10
  major_breaks <- seq(0, xmax_break, by = X_BREAK_BY)
  
  ggplot(df_long, aes(x = SMD, y = var_show)) +
    geom_point(
      aes(shape = Stage),
      size = DOT_SIZE,
      color = "black",
      alpha = 0.95
    ) +
    scale_shape_manual(
      values = c("Before" = 1, "After" = 16),
      guide = guide_legend(
        override.aes = list(size = LEGEND_DOT_SIZE, color = "black")
      )
    ) +
    geom_vline(
      xintercept = 0.1,
      linetype = "dashed",
      linewidth = 0.6,
      color = "#555555"
    ) +
    scale_x_continuous(
      limits = c(0, xmax_break),
      breaks = major_breaks,
      labels = function(x) sprintf("%.1f", x),
      minor_breaks = NULL,
      expand = expansion(mult = c(0.01, 0.03)),
      guide = guide_axis(check.overlap = TRUE)
    ) +
    labs(
      title = title_text,
      x = "|SMD|",
      y = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = LEGEND_TEXT_SIZE),
      legend.title = element_text(size = LEGEND_TEXT_SIZE),
      legend.key.height = grid::unit(LEGEND_KEY_SIZE_PT, "pt"),
      legend.key.width  = grid::unit(LEGEND_KEY_SIZE_PT, "pt"),
      
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "grey85"),
      
      axis.text.x = element_text(size = TEXT_X_SIZE),
      axis.text.y = element_text(size = TEXT_Y_SIZE, face = "bold"),
      axis.title.x = element_text(size = AXIS_TITLE_SIZE),
      
      plot.title = ggtext::element_markdown(hjust = 0.5),
      plot.margin = margin(5, 8, 5, 8),
      
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

## =============================================================================
## 8. LOAD ALL SMD OUTPUTS
## =============================================================================

cem_smd  <- load_method_smd("CEM",  cem_dir)
psm_smd  <- load_method_smd("PSM",  psm_dir)
cbps_smd <- load_method_smd("CBPS", cbps_dir)

smd_all <- rbindlist(
  list(cem_smd, psm_smd, cbps_smd),
  use.names = TRUE,
  fill = TRUE
)

smd_all[, method := factor(toupper(as.character(method)), levels = methods)]
smd_all[, region := as.character(region)]
smd_all[, unit := as.character(unit)]
smd_all[, q := as.character(q)]
smd_all[!grepl("^q", q), q := paste0("q", q)]
smd_all[, d_win := as.integer(d_win)]
smd_all[, before := abs(as.numeric(before))]
smd_all[, after  := abs(as.numeric(after))]

if (ONLY_RT) {
  smd_all <- smd_all[unit == "rt"]
}

diag_summary <- smd_all[
  ,
  .(
    n_variables = uniqueN(variable),
    n_rows = .N,
    max_before = max(before, na.rm = TRUE),
    max_after  = max(after,  na.rm = TRUE)
  ),
  by = .(method, region, unit, q, d_win)
][order(region, unit, q, d_win, method)]

diag_summary[is.infinite(max_before), max_before := NA_real_]
diag_summary[is.infinite(max_after),  max_after  := NA_real_]

fwrite(
  diag_summary,
  file.path(out_dir, "SMD_3x3_input_diagnostic_summary.csv")
)

print(diag_summary)

## =============================================================================
## 9. GENERATE 3x3 PANELS
## =============================================================================

units_to_plot <- sort(unique(smd_all$unit))
if (ONLY_RT) units_to_plot <- "rt"

for (reg in regions_to_plot) {
  for (unit_val in units_to_plot) {
    for (q_tag in q_to_plot) {
      
      cell_var <- list()
      xmax <- 0
      
      for (m in methods) {
        for (w in wins) {
          
          key <- paste(m, w, sep = "_")
          
          var_df <- make_cell_var_df(
            master_dt  = smd_all,
            method_val = m,
            region_val = reg,
            unit_val   = unit_val,
            win_val    = w,
            q_tag      = q_tag
          )
          
          cell_var[[key]] <- var_df
          
          if (!is.null(var_df) && nrow(var_df)) {
            local_max <- var_df[, max(c(before, after), na.rm = TRUE)]
            if (is.finite(local_max)) xmax <- max(xmax, local_max)
          }
        }
      }
      
      if (!is.finite(xmax) || xmax <= 0) xmax <- MIN_XMAX
      xmax <- max(xmax, MIN_XMAX)
      
      mk <- function(m, w) {
        ttl <- make_panel_title(m, w)
        
        build_love_plot(
          var_df     = cell_var[[paste(m, w, sep = "_")]],
          title_text = ttl,
          temp_mode  = unit_val,
          q_tag      = q_tag,
          xlim_max   = xmax,
          d_win      = w
        )
      }
      
      p11 <- mk("CEM",  3)
      p12 <- mk("CEM",  7)
      p13 <- mk("CEM", 14)
      
      p21 <- mk("PSM",  3)
      p22 <- mk("PSM",  7)
      p23 <- mk("PSM", 14)
      
      p31 <- mk("CBPS",  3)
      p32 <- mk("CBPS",  7)
      p33 <- mk("CBPS", 14)
      
      panel <- ggarrange(
        p11, p12, p13,
        p21, p22, p23,
        p31, p32, p33,
        ncol = 3,
        nrow = 3,
        common.legend = TRUE,
        legend = "right",
        align = "hv"
      )
      
      out_base <- sprintf("SMD_3x3_%s_%s_%s", reg, unit_val, q_tag)
      
      out_png <- file.path(out_dir, paste0(out_base, ".png"))
      out_pdf <- file.path(out_dir, paste0(out_base, ".pdf"))
      
      ggsave(
        filename = out_png,
        plot = panel,
        width = 22,
        height = 14,
        units = "in",
        dpi = 600,
        bg = "white"
      )
      
      ggsave(
        filename = out_pdf,
        plot = panel,
        width = 22,
        height = 14,
        units = "in",
        device = cairo_pdf,
        bg = "white"
      )
      
      message(
        "Saved: ", out_png,
        " | region=", reg,
        ", unit=", unit_val,
        ", q=", q_tag
      )
    }
  }
}

message("\nAll updated SMD 3x3 panels completed.")
message("Output directory: ", out_dir)
