
###Positive control for leakage

library(PTPM)
example_data <- load_ptpm_example()
common_args <- list(
  expression = example_data$expression,
  metadata = example_data$metadata,
  k = 3,
  repeats = 20,
  n_boot = 2000,
  seed = 123,
  n_cores = parallel::detectCores() - 1
)

models <- list(
  genome_wide        = list(fun = run_ptpm, extra = list(top_k_values = c(3000, 300))),
  egfr_axis          = list(fun = run_ptpm_egfr),
  pathways           = list(fun = run_ptpm_pathways),
  weighted_score     = list(fun = run_ptpm_weighted_score),
  hybrid             = list(fun = run_ptpm_hybrid),
  hybrid_interaction = list(fun = run_ptpm_hybrid_interactions)
)
pc <- run_ptpm(expression = example_data$expression, metadata = example_data$metadata,
               top_k_values = c(3000, 300), fs_method = "correlation",
               k = 3, repeats = 20, n_boot = 2000, seed = 123,
               n_cores = parallel::detectCores() - 1)
pc$contrast_intervals 