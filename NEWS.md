# PTPM 0.1.5

* Six model functions share one protocol: `run_ptpm()` (genome-wide),
  `run_ptpm_egfr()` (EGFR-axis genes), `run_ptpm_pathways()` (non-outcome
  Hallmark pathways), `run_ptpm_weighted_score()` (PC1-weighted outcome),
  `run_ptpm_hybrid()` (pathways + EGFR-axis genes) and
  `run_ptpm_hybrid_interactions()` (pathway interactions + EGFR-axis genes).
* Nested, patient-grouped cross-validation (outer and inner folds).
* Coupling adjustment re-estimated within outer and inner training folds.
* R² defined as the mean over repeats of 1 - SSE/SST on out-of-fold predictions.
* Patient-level bootstrap confidence intervals and corrected resampled t
  intervals; paired contrasts of the 2 x 2 design.
* Sensitivity analyses: training-reference GSVA (validated against
  `GSVA::gsva()`), unnormalised ssGSEA, pathway-specific adjustment, KRAS-free
  phenotype, KRAS-gene and pathway-overlap exclusion, coupling null.
* Shared GSVA helpers moved to a single internal file (`R/gsva_reference.R`).
* Regimes renamed: standard, coupling-adjusted, fold-internal, fully controlled.
