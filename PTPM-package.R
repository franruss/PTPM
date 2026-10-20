#' PTPM: Paired Transcriptomic Perturbation Modeling
#'
#' Protocol for estimating how well baseline transcriptomic profiles predict
#' a paired perturbation phenotype without optimism from feature-selection
#' leakage or mathematical coupling. Six functions apply the same protocol to
#' different predictor or outcome models: [run_ptpm()], [run_ptpm_egfr()],
#' [run_ptpm_pathways()], [run_ptpm_weighted_score()], [run_ptpm_hybrid()] and
#' [run_ptpm_hybrid_interactions()]. [run_ptpm_benchmark()] re-evaluates the
#' genome-wide model with alternative learners on the same folds.
#'
#' The functions accept any dataset with a baseline and a treated expression
#' profile per paired unit through the `design` argument. [load_ptpm_example()]
#' loads the XENTURION example data (GSE204805), and [load_ptpm_geo()] /
#' [load_ptpm_gse109756()] prepare paired GEO microarray series.
#'
#' @keywords internal
"_PACKAGE"

# Loop and iterator variables used inside foreach() expressions
utils::globalVariables(c("i", "t"))
