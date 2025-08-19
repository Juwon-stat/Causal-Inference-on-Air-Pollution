# SMD Codes

This folder contains R scripts for balance diagnostics using Standardized Mean Difference (SMD).

## Files
- **CBPS_SMD.R** : Summarizes SMD results and generates 2×3 panel plots for CBPS.
- **PSM_SMD.R**  : Summarizes SMD results and generates 2×3 panel plots for PSM.
- **CEM_SMD.R**  : Summarizes SMD results and generates 2×3 panel plots for CEM.
- **Combined_3x3_Panel.R** : Combines results from CEM, PSM, and CBPS into 3×3 comparison panels.

## Features
- Exports SMD values as CSV (wide format, tidy format, optional raw level).  
- Generates **Love plots** (before/after balance).  
- Produces **2×3 panels** (region × time window).  
- Produces **3×3 panels** (methods × time window) for cross-method comparison.

## Workflow
1. After running the matching scripts, run the corresponding SMD script.  
2. CSV files and panel PNG plots are saved into the appropriate directories:
   - `SMD_Values/`
   - `SMD_Values_tidy/`
   - `SMD_Plots/`
