# Matching Codes

This folder contains R scripts for applying different matching algorithms to air pollution data.

## Files
- **CBPS.R** : Covariate Balancing Propensity Score (CBPS) implementation.
- **PSM.R**  : Propensity Score Matching (PSM) implementation.
- **CEM.R**  : Coarsened Exact Matching (CEM) implementation.

## Description
Each script loads raw data, applies the specified matching method, and outputs matched datasets.  
The matched datasets are later used for balance diagnostics (see `../smd_code/`).

## Notes
- Input data should be prepared beforehand (see `data_dir` path settings inside the scripts).  
- Matching results are stored as `.RData` files in the corresponding method folder under `Matching_Results/`.
