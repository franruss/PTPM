## Full analysis (manuscript settings)

#The manuscript settings are `k = 3`, `repeats = 20`, `n_boot = 2000` and
#`seed = 123`, with each function's default K grid. For `run_ptpm()` the K values
#are `top_k_values = c(3000, 300)`.

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

results <- list()
for (nm in names(models)) {
  cat("\n#####", nm, "#####\n")
  t0 <- Sys.time()
  results[[nm]] <- do.call(models[[nm]]$fun, c(common_args, models[[nm]]$extra))
  print(Sys.time() - t0)
  saveRDS(results[[nm]], paste0("ptpm_", nm, ".rds"))
}

results

