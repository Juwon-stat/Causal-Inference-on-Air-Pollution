How to run

Requirements

R (≥ 4.1 recommended)

Optional: data.table (script will fall back to base R if not installed)

Expected bundle contents

Each RData must contain:

TALL: treated-side unit table with columns rid, Y10, Y25, N

CALL: control-side unit table with columns rid, wY10, wY25, w

(Optional) att_over: any object you previously saved (preserved as-is)

Semantics:

rid = unit/cluster identifier used for unit-level resampling or reweighting

N = treated-side denominator contribution (e.g., unit weight or count)

w = control-side denominator contribution (e.g., weighted count)

Y10, Y25 = treated outcomes for PM10 / PM2.5

wY10, wY25 = control outcomes already multiplied by their control weights

File layout

Set in_dir to the folder containing your CBPS bundles whose names match:

^cbps_UNIT_ATT_boot_.*_WITH_SMD_DIAG\\.RData$

Outputs will be written to file.path(in_dir, "posthoc_se"):

One enriched *_POSTSE.RData per input

One CSV summary across all inputs:

POSTHOC_ATT_SE_summary_dist.csv (change filename if you like)

Run
