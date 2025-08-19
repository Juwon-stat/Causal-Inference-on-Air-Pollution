# Causal-Inference-on-Air-Pollution
Causal Inference on Long-Range Air Pollution: A Matching-Based Study of Seoul’s Particulate Matter During China’s COVID-19 Lockdown (for Scientific Reports)

- **seoul_q3.RData**: Data used for matching based on quantile grouping (Q3, ±3 days)
- **seoul_q4.RData**: Data for Q4 matching window
- **seoul_q5.RData**: Data for Q5 matching window

## Each seoul_q*.RData file contains data.frame named:
- `seoul_q*`: Air pollution and weather data (2017–2023)


## 📂 Repository Structure
This repository contains R code and materials for causal inference analysis of air pollution in Seoul, focusing on matching methods and balance diagnostics.

Causal-Inference-on-Air-Pollution/
├── R/
│ ├── matching_code/ # Matching algorithms (PSM, CBPS, CEM)
│ ├── smd_code/ # Standardized Mean Difference (SMD) summaries and panel plots
│ └── README.md
└── README.md

## 🔎 Contents

- **Matching codes**  
  - Propensity Score Matching (PSM)  
  - Covariate Balancing Propensity Score (CBPS)  
  - Coarsened Exact Matching (CEM)  

- **Balance diagnostics**  
  - Standardized Mean Difference (SMD) calculation  
  - Panel visualization (Love plots, 2×3 and 3×3 panels)  
  - CSV exports of summarized balance tables  

## ⚙️ Requirements
- R (≥ 4.0)
- Required R packages:
  - `data.table`, `ggplot2`, `ggpubr`, `forcats`, `dplyr`, `tidyr`, `lubridate`

## 🚀 Usage
1. Run matching codes (`R/matching_code/`) to generate matched datasets.  
2. Run SMD codes (`R/smd_code/`) to compute balance diagnostics and plots.  
3. Combined 3×3 panel plots can be created to compare methods (CEM / PSM / CBPS).
