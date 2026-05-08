# Block-wise Intertemporal Matching for Causal Inference on Meteorological Data
## : Quantifying the Impact of Upwind Lockdowns on Downwind Air Quality**

This repository provides the full reproducible R code and analysis materials for a methodological causal inference study on transboundary air pollution. The study develops and applies a **block-wise intertemporal matching framework** to estimate the causal impact of upstream lockdown-related emission reductions on downstream particulate matter concentrations in Seoul.

The main objective of this repository is to document the complete computational workflow used in the paper, including block-wise matching, covariate balance assessment, bootstrap-based uncertainty evaluation, and visualization of standardized mean differences.

---

## Overview

This study treats the early 2020 COVID-19 lockdown period in China as a quasi-experimental intervention and evaluates its downstream impact on air quality in Seoul. Because air pollution is strongly affected by meteorological conditions, seasonal patterns, temporal dependence, and local emissions, the analysis uses a block-wise intertemporal matching design rather than a simple before-after or regression-based comparison.

The repository contains full implementation code for the three matching strategies considered in the paper:

- Coarsened Exact Matching (CEM)
- Propensity Score Matching (PSM)
- Covariate Balancing Propensity Score (CBPS)

Each method-specific code file contains the complete workflow from matching construction to covariate balance diagnostics, including SMD calculation and bootstrap-based evaluation.

---

## Data Files

The repository includes the following preprocessed Seoul air pollution and meteorological datasets:

- **seoul_q3.RData**  
  Data used for matching based on 3-quantile grouping.

- **seoul_q4.RData**  
  Data used for matching based on quartile grouping.

- **seoul_q5.RData**  
  Data used for matching based on 5-quantile grouping.

Each `.RData` file contains a data frame named according to its quantile level:

- `seoul_q3`
- `seoul_q4`
- `seoul_q5`

These datasets contain hourly air pollution and meteorological information for Seoul from January to February across 2017–2023.

---

## Repository Contents

### 1. Method-specific full analysis code

The `R/matching_code/` directory contains the full implementation code for each matching method:

- CEM full code
- PSM full code
- CBPS full code

Each method-specific code includes:

- construction of treatment and control observations
- block-wise intertemporal matching
- covariate balance assessment
- standardized mean difference calculation
- bootstrap-based SMD calculation
- export of matched datasets and diagnostic summaries

### 2. SMD and balance diagnostic code

The `R/smd_code/` directory contains code for producing balance diagnostic outputs from the matched datasets.

This includes:

- SMD calculation before and after matching
- bootstrap-based balance summaries
- Love plots
- method-specific SMD plots
- combined panel plots comparing CEM, PSM, and CBPS

### 3. Visualization outputs

The code can be used to generate graphical diagnostics for comparing covariate balance across methods, quantile groupings, and matching specifications.

The main visualization outputs include:

- method-specific SMD plots
- before-and-after matching balance plots
- 3×3 panel plots comparing CEM, PSM, and CBPS

---

## Methodological Focus

This repository is intended to support the methodological contribution of the paper. The central idea is not merely to apply existing matching methods, but to adapt them to a time-indexed environmental causal inference setting where treatment and control observations are constructed across comparable temporal blocks.

The proposed framework emphasizes:

- treatment-day-centered construction of comparison sets
- block-wise matching under temporal and meteorological comparability
- balance diagnostics tailored to reused control observations
- bootstrap-based assessment of balance stability
- transparent comparison of CEM, PSM, and CBPS within the same empirical design

This structure allows the analysis to address key challenges in air pollution causal inference, including seasonality, meteorological confounding, temporal dependence, and limited treated periods.

---

## Requirements

The analysis was conducted in R.

Required R packages include:

- `data.table`
- `dplyr`
- `tidyr`
- `lubridate`
- `ggplot2`
- `ggpubr`
- `forcats`
- `CBPS`
- `MatchIt`

Additional packages may be required depending on the specific matching or plotting script.

---

## Usage

1. Load one of the preprocessed datasets:

   - `seoul_q3.RData`
   - `seoul_q4.RData`
   - `seoul_q5.RData`

2. Run the method-specific full code in `R/matching_code/`.

3. Use the SMD diagnostic scripts in `R/smd_code/` to compute and visualize covariate balance.

4. Generate method-specific or combined SMD plots to compare the performance of CEM, PSM, and CBPS.

---

## Citation

This repository accompanies the methodological paper:

**Causal Inference on Transboundary Air Pollution: Quantifying the Impact of Upwind Lockdowns on Downwind Air Quality**

Please cite the paper when using this repository or adapting the code.
