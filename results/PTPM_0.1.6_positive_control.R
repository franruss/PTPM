## =====================================================================
## PTPM 0.1.6: positive control for feature-selection leakage
## =====================================================================
## Genes are ranked by absolute correlation with the outcome instead of
## by variance (fs_method = "correlation"). Ranked on the full dataset,
## this supervised selection uses the outcome of held-out patients, so
## the 2 x 2 design should show a large feature-selection component.
## Reproduces the last row of Table 1 (user guide, section 6.1).
## =====================================================================

library(PTPM)
example_data <- load_ptpm_example()

pc <- run_ptpm(
  expression   = example_data$expression,
  metadata     = example_data$metadata,
  top_k_values = c(3000, 300),
  fs_method    = "correlation",
  k = 3, repeats = 20, n_boot = 2000, seed = 123,
  n_cores = parallel::detectCores() - 1
)
saveRDS(pc, "ptpm_positive_control.rds")   # used by PTPM_0.1.6_scale_of_S_and_Oldham.R

pc$confidence_intervals
pc$contrast_intervals   # Delta_FS_no_coupling_adj: leakage from supervised selection
