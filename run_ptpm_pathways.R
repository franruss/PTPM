#' Paired Transcriptomic Perturbation Modeling: pathway-level predictor model
#'
#' Same protocol as run_ptpm(), with a different predictor model: the
#' elastic net uses as predictors the baseline (NT) GSVA scores of the
#' Hallmark pathways that do NOT define the phenotype (by default the 50
#' Hallmark sets minus E2F_TARGETS, KRAS_SIGNALING_UP and MTORC1_SIGNALING,
#' i.e. 47 pathways), instead of baseline gene expression.
#'
#' Excluding the three outcome pathways removes the direct algebraic link
#' between predictors and outcome, but not all coupling: (i) predictor and
#' outcome pathways share genes, and (ii) predictor pathways are correlated
#' with the outcome pathways. The coupling adjustment, the coupling
#' diagnostics and a robustness analysis excluding predictor pathways with
#' high gene overlap with the outcome sets address this residual coupling.
#'
#' Sensitivity on K: the top-K pathways by variance are retained for a grid
#' of K values (default 10, 20, 30, 40 and all predictor pathways). For K
#' equal to all pathways no feature selection takes place, so the standard
#' and fold-internal regimes (and the coupling-adjusted and fully controlled
#' regimes) coincide by construction.
#'
#' Predictor scores are computed in the same GSVA run as the outcome scores
#' (GSVA scores of a gene set do not depend on the other sets analysed). In
#' the training-reference GSVA sensitivity analysis and in the ssGSEA
#' sensitivity analysis, predictors are recomputed with the same scoring as
#' the outcome.
#'
#' Runs the PTPM workflow on a transcriptomic expression matrix and sample
#' metadata. Paired NT/CET samples are identified, pathway scores are obtained
#' with GSVA, the suppression phenotype is defined as NT minus CET, and
#' baseline expression is evaluated with repeated, patient-grouped, nested
#' cross-validation and elastic-net regression.
#'
#' EVALUATION DESIGN
#' -----------------
#' Four regimes form a 2 x 2 factorial design. All regimes share the same
#' patient-grouped outer folds, the same patient-grouped inner folds and the
#' same elastic-net engine, so they differ only in the two factors below.
#'
#'   Regime (key)                          Feature selection (FS)        Coupling adjustment
#'   Standard           (standard)          pre-CV, on the full dataset   no
#'   Coupling-adjusted  (coupling_adjusted) pre-CV, on the full dataset   yes
#'   Fold-internal      (fold_internal)     inside outer training folds   no
#'   Fully controlled   (fully_controlled)  inside outer training folds   yes
#'
#' The term "leakage-free" is deliberately not used: GSVA scores are computed
#' on the whole cohort (see below), so no regime is strictly free of
#' cross-sample information.
#'
#' Coupling adjustment (fully nested): S is regressed on the baseline pathway
#' score b in the outer training fold; the fitted coefficients are applied
#' unchanged to the held-out fold, and the model predicts the residual. During
#' hyperparameter tuning the regression is re-estimated within each inner
#' training fold and applied to the corresponding inner validation fold, so no
#' residual used for tuning depends on the observations it is evaluated on.
#'
#' DEFINITION OF "FULLY CONTROLLED"
#' --------------------------------
#' An estimate is called fully controlled when every supervised or
#' data-dependent step that could carry information about a held-out patient
#' is estimated on training patients only:
#'   (i)   outer train/test splits at patient level;
#'   (ii)  inner (tuning) splits at patient level, within the outer training set;
#'   (iii) gene ranking (FS) on the outer training set;
#'   (iv)  elastic-net hyperparameters (alpha, lambda) by inner CV only;
#'   (v)   coupling regression fitted on the outer training set only.
#' "Fold-internal" satisfies (i)-(iv) without coupling adjustment.
#' Documented exception: GSVA pathway scores, which define S and b, are
#' computed once on the whole cohort (unsupervised, transductive). Their
#' influence is quantified by the training-reference GSVA sensitivity analysis.
#'
#' GSVA SENSITIVITY (TRAINING-REFERENCE GSVA)
#' ------------------------------------------
#' GSVA has one cross-sample step (per-gene CDF, Gaussian kernel or direct
#' ECDF) and per-sample steps (within-sample ranking, random walk). In each
#' outer fold the cross-sample step is estimated on the NT + CET profiles of
#' the training models only and applied to all samples, so training and test
#' samples are scored with the same training-derived transformation and
#' held-out patients do not contribute to it. With reference = all samples
#' the implementation must reproduce GSVA::gsva() (checked at run time; the
#' matching kernel type and tie handling are then used in every fold).
#' Running GSVA::gsva() separately on training and test samples would instead
#' place them on different scales and is not a valid sensitivity analysis.
#'
#' PERFORMANCE METRIC
#' ------------------
#' For each repeat, R2_r = 1 - SSE_r/SST_r on the pooled out-of-fold
#' predictions of that repeat (its k folds; each sample predicted exactly once;
#' SST_r around the mean of the observed values of that repeat). The reported
#' R2 is the arithmetic MEAN of the 'repeats' values R2_r; it is not a single
#' SSE/SST pooled over all repeats, nor a mean of per-fold R2. The
#' squared correlation (cor2) is also reported, and the legacy pooled cor2
#' over all repeats is kept for comparison with earlier results. Differences
#' between regimes are computed per repeat (paired: same folds) and summarised
#' by mean, SD, 2.5-97.5% range and the share of repeats with difference > 0.
#' Repeats reuse the same data and are not independent: these summaries
#' describe partition variability, not sampling uncertainty.
#'
#' UNCERTAINTY (CONFIDENCE INTERVALS)
#' ----------------------------------
#' Two complementary intervals are reported for the R2 of each regime and for
#' each paired contrast between regimes:
#'   (a) Patient-level cluster bootstrap of the out-of-fold predictions
#'       (n_boot resamples of patients with replacement, all models of a
#'       patient kept together). In each resample R2 is recomputed per repeat
#'       and averaged over repeats; the same resample is used for all regimes,
#'       so contrasts stay paired. Percentile intervals. Fitted models are
#'       held fixed, so training variability is not included (intervals may be
#'       somewhat narrow).
#'   (b) Corrected resampled t interval for repeated k-fold CV
#'       (Nadeau & Bengio 2003; Bouckaert & Frank 2004): per-fold R2 values,
#'       variance inflated by (1/J + n_test/n_train) to account for
#'       overlapping training sets; t distribution with J - 1 df.
#'       Note: per-fold R2 uses the test-fold mean in SST and is computed on
#'       small test sets, so its mean (nb_mean) is lower and noisier than the
#'       pooled per-repeat R2. Use (a) as the primary interval and (b) as a
#'       conservative check.
#'
#' ROBUSTNESS (robustness = TRUE)
#' ------------------------------
#' All with fold-internal FS and the same folds as the main analysis:
#'   (1) pathway-specific coupling adjustment: S residualised on the baseline
#'       score of each pathway (instead of their mean), within training folds;
#'   (2) KRAS-free phenotype: S and b rebuilt without the KRAS gene set, so
#'       KRAS is no longer part of the outcome (with and without adjustment);
#'   (3) KRAS genes excluded from the baseline predictors;
#'   (4) partial correlation of the baseline KRAS score with the full and the
#'       KRAS-free phenotype, adjusted for the phenotype's own baseline score
#'       (mean, and, as a secondary analysis, each non-KRAS baseline pathway
#'       separately; in-sample, descriptive, patient-level bootstrap CI).
#'
#' SECONDARY SENSITIVITY: ssGSEA (ssgsea_sensitivity = TRUE)
#' ---------------------------------------------------------
#' S and b are rebuilt from ssGSEA scores without normalisation, a truly
#' single-sample method (no cross-sample information by construction). The
#' four regimes and the pathway-specific adjustment are re-run on the same
#' folds and compared with the GSVA-based R2 (paired, patient bootstrap).
#' The coupling null of cor(b, S) and the KRAS partial correlations are also
#' recomputed with ssGSEA scores.
#' Differences vs the main regimes are paired per repeat, with patient-level
#' bootstrap intervals.
#'
#' PERMUTATION NULL OF THE NT/CET PAIRING (optional, n_perm > 0)
#' --------------------------------------------------------------
#' CET scores are shuffled across models (S_perm = b - `c[perm]`), which keeps
#' the arithmetic coupling of S with b and destroys model-specific pairing.
#' The fold-internal and fully controlled regimes are re-run on S_perm. The
#' null R2 of the fold-internal regime is the predictability expected from
#' coupling alone; under the null the fully controlled R2 should be near 0.
#' A cheap null of cor(b, S) with n_perm_coupling permutations is always run.
#'
#' INPUT FOR OTHER DATASETS
#' ------------------------
#' Any dataset with a baseline and a treated expression profile for each
#' paired unit can be analysed by passing `design`. Expression must be on a
#' log-like scale (e.g. log2 microarray intensities or variance-stabilised
#' counts) with gene symbols as identifiers. Replicate samples of the same
#' pair and condition are averaged; pairs without both profiles are excluded.
#' The suppression phenotype is S = baseline score - treated score, so for a
#' treatment that activates the outcome gene sets S is negative. The default
#' outcome gene sets are those of the anti-EGFR analysis; supply `gene_sets`
#' (and `kras_set`, or set `robustness = FALSE`) for other treatments.
#' Internally, baseline and treated profiles are labelled "NT" and "CET".
#'
#' @param expression Expression matrix or data frame. By default, the first
#'   column contains gene identifiers and the remaining columns are samples.
#' @param metadata Data frame containing a `title` column with the sample
#'   identifiers/titles used to derive treatment and pair information. Used
#'   only for the XENTURION example data (`design = NULL`).
#' @param design Data frame describing the samples of any paired dataset (see
#'   INPUT FOR OTHER DATASETS). Columns: `sample` (expression column name),
#'   `pair` (paired unit: model, tumour or subject), `condition` (values
#'   `baseline_label` or `treated_label`; other values are ignored) and,
#'   optionally, `patient` (grouping of the cross-validation; defaults to
#'   `pair`). If `NULL`, pairs are parsed from `metadata$title` as in the
#'   XENTURION example.
#' @param baseline_label,treated_label Values of `design$condition` that
#'   identify the baseline and the treated profile of each pair.
#' @param balanced_folds Logical; build the outer folds by dealing the
#'   shuffled patients to the k folds in turn, so that every fold has the same
#'   number of patients (+/- 1). `NULL` (default) means `TRUE` with `design`
#'   and `FALSE` for the XENTURION example, where `caret::groupKFold()` is kept
#'   to reproduce the published results (it samples fold labels with
#'   replacement, so fold sizes vary and a fold can be empty in small cohorts).
#' @param gene_prefix Prefix that marks the human genes in the expression
#'   identifiers (removed before scoring). `NULL` (default) means `"H_"` for
#'   the XENTURION example and no prefix with `design`; use `""` for none.
#' @param gene_sets Optional named list of the OUTCOME gene sets. If `NULL`,
#'   the three Hallmark sets of the original PTPM analysis are retrieved with
#'   `msigdbr`.
#' @param predictor_gene_sets Optional named list of PREDICTOR gene sets. If
#'   `NULL`, all Hallmark gene sets except the outcome sets are used.
#' @param overlap_threshold Robustness analysis: predictor pathways whose
#'   fraction of genes shared with the outcome sets exceeds this value are
#'   excluded (default 0.2).
#' @param top_k_values Integer vector with the numbers of highest-variance
#'   predictor pathways evaluated. If `NULL` (default), c(10, 20, 30, 40,
#'   <number of predictor pathways>) is used; larger values are truncated.
#' @param k Number of outer grouped folds.
#' @param repeats Number of repeated outer cross-validation runs.
#' @param n_inner Number of inner grouped folds used for tuning (>= 3).
#' @param alpha_grid Elastic-net mixing values explored in the inner CV.
#' @param gsva_sensitivity Logical; run the training-reference GSVA
#'   sensitivity analysis (see GSVA SENSITIVITY).
#' @param gsva_sensitivity_repeats Number of repeats (from the first) used in
#'   the GSVA sensitivity analysis. Each fold requires one GSVA scoring.
#' @param robustness Logical; run the robustness analyses (pathway-specific
#'   coupling adjustment and KRAS circularity controls; see ROBUSTNESS).
#' @param robustness_repeats Number of repeats (from the first) used in the
#'   robustness analyses.
#' @param kras_set Name of the KRAS gene set in `gene_sets`.
#' @param ssgsea_sensitivity Logical; run the secondary sensitivity analysis
#'   with ssGSEA (normalize = FALSE) instead of GSVA.
#' @param ssgsea_repeats Number of repeats (from the first) used in the ssGSEA
#'   sensitivity analysis.
#' @param ssgsea_alpha Weight exponent of the ssGSEA random walk.
#' @param gsva_ref_tol Tolerance for the check that the reference-based GSVA,
#'   with reference = all samples, reproduces GSVA::gsva().
#' @param n_perm Number of permutations of the NT/CET pairing for the full
#'   cross-validated null (0 = skip; computationally heavy).
#' @param perm_repeats Number of repeats (from the first) used per permutation.
#' @param n_perm_coupling Number of permutations for the cor(b, S) null.
#' @param n_boot Number of patient-level bootstrap resamples (0 = skip).
#' @param conf_level Confidence level of the intervals.
#' @param seed Integer random seed.
#' @param n_cores Number of parallel workers.
#' @param expr_matrix Deprecated alias for `expression`.
#'
#' @return A list; see section 26 for its elements. Main tables:
#'   `results_genome_wide`, `results_summary`, `confidence_intervals`,
#'   `contrast_intervals`, `paired_differences`,
#'   `gsva_sensitivity`, `permutation_null`, `coupling_null`, `leakage_audit`.
#' @export
run_ptpm_pathways <- function(
    expression = NULL,
    metadata = NULL,
    design = NULL,
    baseline_label = "baseline",
    treated_label = "treated",
    gene_prefix = NULL,
    balanced_folds = NULL,
    gene_sets = NULL,
    top_k_values = NULL,
    predictor_gene_sets = NULL,
    overlap_threshold = 0.2,
    k = 3,
    repeats = 20,
    n_inner = 5,
    alpha_grid = seq(0, 1, length.out = 11),
    gsva_sensitivity = TRUE,
    gsva_sensitivity_repeats = repeats,
    gsva_ref_tol = 1e-6,
    robustness = TRUE,
    robustness_repeats = repeats,
    kras_set = "HALLMARK_KRAS_SIGNALING_UP",
    ssgsea_sensitivity = TRUE,
    ssgsea_repeats = repeats,
    ssgsea_alpha = 0.25,
    n_perm = 0,
    perm_repeats = 1,
    n_perm_coupling = 1000,
    n_boot = 2000,
    conf_level = 0.95,
    seed = 123,
    n_cores = max(1, parallel::detectCores() - 1),
    expr_matrix = NULL
) {
  # 0. ARGUMENT CHECKS ---------------------------------------------------
  if (!is.null(expr_matrix)) {
    if (!is.null(expression)) {
      stop("Supply either 'expression' or the deprecated 'expr_matrix', not both.")
    }
    expression <- expr_matrix
  }
  if (is.null(expression)) stop("'expression' is required.")
  if (!(is.matrix(expression) || is.data.frame(expression))) {
    stop("'expression' must be a matrix or data frame.")
  }
  if (is.null(design) && is.null(metadata)) {
    stop("Supply 'design' (any paired dataset) or 'metadata' (XENTURION example data).")
  }

  is_int_scalar <- function(v, min_v) {
    length(v) == 1L && is.finite(v) && v >= min_v && v == as.integer(v)
  }
  if (!is.null(top_k_values) &&
      (length(top_k_values) == 0 || any(!is.finite(top_k_values)) ||
       any(top_k_values < 1) || any(top_k_values != as.integer(top_k_values)))) {
    stop("'top_k_values' must be NULL or contain positive integers.")
  }
  if (!is_int_scalar(k, 2)) stop("'k' must be a single integer >= 2.")
  if (!is_int_scalar(repeats, 1)) stop("'repeats' must be a single integer >= 1.")
  if (!is_int_scalar(n_inner, 3)) stop("'n_inner' must be a single integer >= 3.")
  if (!is_int_scalar(n_cores, 1)) stop("'n_cores' must be a single integer >= 1.")
  if (!is_int_scalar(n_perm, 0)) stop("'n_perm' must be a single integer >= 0.")
  if (!is_int_scalar(n_perm_coupling, 0)) stop("'n_perm_coupling' must be a single integer >= 0.")
  if (!is_int_scalar(n_boot, 0)) stop("'n_boot' must be a single integer >= 0.")
  if (length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be a single number in (0, 1).")
  }
  if (!is_int_scalar(perm_repeats, 1) || perm_repeats > repeats) {
    stop("'perm_repeats' must be an integer between 1 and 'repeats'.")
  }
  if (!is_int_scalar(gsva_sensitivity_repeats, 1) || gsva_sensitivity_repeats > repeats) {
    stop("'gsva_sensitivity_repeats' must be an integer between 1 and 'repeats'.")
  }
  if (length(seed) != 1L || !is.finite(seed)) stop("'seed' must be a single finite number.")
  if (!is_int_scalar(ssgsea_repeats, 1) || ssgsea_repeats > repeats) {
    stop("'ssgsea_repeats' must be an integer between 1 and 'repeats'.")
  }
  if (!is_int_scalar(robustness_repeats, 1) || robustness_repeats > repeats) {
    stop("'robustness_repeats' must be an integer between 1 and 'repeats'.")
  }
  if (!is.numeric(alpha_grid) || any(alpha_grid < 0 | alpha_grid > 1)) {
    stop("'alpha_grid' must contain values in [0, 1].")
  }

  for (pkg in c("glmnet", "caret", "msigdbr", "GSVA", "foreach", "doParallel")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
  }
  `%dopar%` <- foreach::`%dopar%`

  cat("\n============================================================\n")
  cat("PTPM analysis\n")
  cat("============================================================\n")
  cat("predictor model   : baseline GSVA scores of non-outcome pathways\n")
  cat("top_k_values      :", if (is.null(top_k_values)) "auto (10, 20, 30, 40, all)" else
    paste(top_k_values, collapse = ", "), "\n")
  cat("outer folds       :", k, " x repeats:", repeats, "=", k * repeats, "splits\n")
  cat("inner folds       :", n_inner, "(grouped by patient)\n")
  cat("regimes           : 4 (FS position x coupling correction)\n")
  cat("GSVA sensitivity  :", if (gsva_sensitivity) paste(gsva_sensitivity_repeats, "repeats") else "off", "\n")
  cat("pairing null      :", if (n_perm > 0) paste(n_perm, "permutations x", perm_repeats, "repeat(s)") else "off", "\n")
  cat("seed / cores      :", seed, "/", n_cores, "\n")
  cat("============================================================\n\n")

  # 1. EXPRESSION MATRIX -------------------------------------------------
  expression <- as.data.frame(expression, check.names = FALSE,
                              stringsAsFactors = FALSE)
  if (ncol(expression) < 2L) {
    stop("'expression' must contain a gene-ID column and at least one sample column.")
  }
  rownames(expression) <- as.character(expression[[1L]])
  expression <- expression[, -1L, drop = FALSE]
  expression <- as.matrix(expression)
  suppressWarnings(mode(expression) <- "numeric")
  if (anyNA(expression)) {
    stop("The expression matrix contains non-numeric or missing values after conversion.")
  }
  if (anyDuplicated(rownames(expression))) {
    stop("Gene identifiers in the first expression column must be unique.")
  }
  expr_ids <- colnames(expression)
  if (is.null(expr_ids) || anyNA(expr_ids) || any(expr_ids == "")) {
    stop("Expression sample columns must have non-empty names.")
  }

  # 2-7. PAIRED PROFILES (AVERAGED PER PAIR AND CONDITION) --------------
  # Design mode (any dataset) or legacy XENTURION parsing of metadata$title;
  # see prepare_pairs.R. Columns of expr_avg: "<pair>_NT" (baseline) and
  # "<pair>_CET" (treated).
  prep <- .ptpm_prepare_pairs(expression, metadata = metadata, design = design,
                              baseline_label = baseline_label,
                              treated_label = treated_label,
                              gene_prefix = gene_prefix, k = k)
  expr_avg <- prep$expr_avg
  pair_patient <- prep$pair_patient
  gene_prefix <- prep$gene_prefix
  if (is.null(balanced_folds)) balanced_folds <- prep$mode == "design"
  if (!is.logical(balanced_folds) || length(balanced_folds) != 1L || is.na(balanced_folds)) {
    stop("'balanced_folds' must be TRUE, FALSE or NULL.")
  }
  cat("Input mode:", prep$mode,
      if (nzchar(gene_prefix)) paste0("| gene prefix '", gene_prefix, "'") else "| no gene prefix",
      "\n")

  # 8. GSVA (GLOBAL, MAIN ANALYSIS) --------------------------------------
  if (is.null(gene_sets)) {
    # msigdbr >= 10.0.0 renamed 'category' to 'collection'
    hallmark <- if ("collection" %in% names(formals(msigdbr::msigdbr))) {
      msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    } else {
      msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    }
    egfr_sets <- hallmark[hallmark$gs_name %in% c(
      "HALLMARK_KRAS_SIGNALING_UP",
      "HALLMARK_MTORC1_SIGNALING",
      "HALLMARK_E2F_TARGETS"
    ), , drop = FALSE]
    gene_sets <- split(egfr_sets$gene_symbol, egfr_sets$gs_name)
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("'gene_sets' must be a named list of gene vectors.")
  }

  human_genes <- .ptpm_has_prefix(rownames(expr_avg), gene_prefix)
  expr_human <- expr_avg[human_genes, , drop = FALSE]
  if (nrow(expr_human) == 0L) {
    stop("No expression rows beginning with '", gene_prefix, "' were found.")
  }
  rownames(expr_human) <- .ptpm_strip_prefix(rownames(expr_human), gene_prefix)

  # Predictor gene sets: all Hallmark sets except the outcome sets
  if (is.null(predictor_gene_sets)) {
    hall_all <- if ("collection" %in% names(formals(msigdbr::msigdbr))) {
      msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    } else {
      msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    }
    predictor_gene_sets <- split(hall_all$gene_symbol, hall_all$gs_name)
  }
  if (!is.list(predictor_gene_sets) || is.null(names(predictor_gene_sets))) {
    stop("'predictor_gene_sets' must be a named list of gene vectors.")
  }
  predictor_gene_sets <- predictor_gene_sets[setdiff(names(predictor_gene_sets),
                                                     names(gene_sets))]
  all_sets <- c(gene_sets, predictor_gene_sets)

  run_gsva <- function(expr_mat, sets = all_sets) {
    param <- GSVA::gsvaParam(expr_mat, sets)
    suppressMessages(GSVA::gsva(param))
  }
  gsva_all <- run_gsva(expr_human)
  outcome_rows <- intersect(names(gene_sets), rownames(gsva_all))
  predictor_rows <- intersect(names(predictor_gene_sets), rownames(gsva_all))
  if (length(outcome_rows) == 0L) stop("No outcome gene set could be scored by GSVA.")
  if (length(predictor_rows) < 2L) stop("Fewer than 2 predictor gene sets could be scored by GSVA.")
  gsva_scores <- gsva_all[outcome_rows, , drop = FALSE]

  # 9. NT / CET ----------------------------------------------------------
  nt_cols <- grepl("_NT$", colnames(gsva_scores))
  cet_cols <- grepl("_CET$", colnames(gsva_scores))
  if (!any(nt_cols) || !any(cet_cols)) stop("Both NT and CET groups are required for GSVA.")
  gsva_nt <- gsva_scores[, nt_cols, drop = FALSE]
  gsva_cet <- gsva_scores[, cet_cols, drop = FALSE]

  # 10. SUPPRESSION PHENOTYPE ------------------------------------------
  # S = mean pathway score in NT - mean pathway score in CET = b - c
  pair_ids <- intersect(sub("_NT$", "", colnames(gsva_nt)),
                        sub("_CET$", "", colnames(gsva_cet)))
  if (length(pair_ids) < k) stop("Too few paired GSVA profiles for grouped cross-validation.")
  score_pairs <- function(scores, pids, sets = rownames(scores)) {
    b_val <- colMeans(scores[sets, paste0(pids, "_NT"), drop = FALSE])
    c_val <- colMeans(scores[sets, paste0(pids, "_CET"), drop = FALSE])
    names(b_val) <- names(c_val) <- pids
    list(b = b_val, c = c_val, S = b_val - c_val)
  }

  # 11. BASELINE EXPRESSION ---------------------------------------------
  baseline_cols <- paste0(pair_ids, "_NT")
  expr_baseline <- expr_avg[, baseline_cols, drop = FALSE]
  colnames(expr_baseline) <- pair_ids
  rownames(expr_baseline) <- .ptpm_strip_prefix(rownames(expr_baseline), gene_prefix)
  x_baseline_matrix <- t(expr_baseline)
  if (anyDuplicated(colnames(x_baseline_matrix))) {
    stop("Gene identifiers are not unique after removing the gene prefix.")
  }
  n_obs <- nrow(x_baseline_matrix)

  # 11b. PATHWAY-LEVEL PREDICTORS: baseline GSVA scores of non-outcome sets --
  x_baseline_matrix <- t(gsva_all[predictor_rows, paste0(pair_ids, "_NT"), drop = FALSE])
  rownames(x_baseline_matrix) <- pair_ids
  n_pred <- ncol(x_baseline_matrix)
  if (is.null(top_k_values)) top_k_values <- c(10, 20, 30, 40, n_pred)
  top_k_values <- sort(unique(pmin(as.integer(top_k_values), n_pred)))
  outcome_genes <- unique(unlist(gene_sets, use.names = FALSE))
  pathway_overlap <- data.frame(
    pathway = predictor_rows,
    n_genes = vapply(predictor_gene_sets[predictor_rows], function(g) length(unique(g)), numeric(1)),
    frac_shared_with_outcome = vapply(predictor_gene_sets[predictor_rows], function(g) {
      g <- unique(g); length(intersect(g, outcome_genes)) / length(g)
    }, numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE)
  cat("Predictor pathways:", n_pred, "(outcome pathways excluded:",
      paste(outcome_rows, collapse = ", "), ") | K grid:",
      paste(top_k_values, collapse = ", "), "\n")
  if (n_pred %in% top_k_values) {
    cat("  note: for K =", n_pred, "(all pathways) no feature selection is performed;",
        "standard = fold-internal and coupling-adjusted = fully controlled by construction.\n")
  }

  # 12. PHENOTYPE, BASELINE SCORE AND CET SCORE (GLOBAL GSVA) -----------
  glob <- score_pairs(gsva_scores, pair_ids)
  y_main <- glob$S[rownames(x_baseline_matrix)]
  b_main <- glob$b[rownames(x_baseline_matrix)]
  c_main <- glob$c[rownames(x_baseline_matrix)]
  if (anyNA(y_main) || anyNA(b_main)) stop("Could not align scores to baseline samples.")

  # 13. PATIENT GROUPS ---------------------------------------------------
  groups <- unname(pair_patient[rownames(x_baseline_matrix)])
  if (anyNA(groups)) stop("Could not assign a patient to every paired model.")
  n_patients <- length(unique(groups))
  cat("Paired models:", n_obs, " | patients (groups):", n_patients, "\n")
  if (n_patients < k) stop("Fewer patients than outer folds.")

  # 14. OUTER CV (GROUPED BY PATIENT) -----------------------------------
  folds <- list()
  for (r in seq_len(repeats)) {
    set.seed(seed + r - 1)
    fold_ids <- if (balanced_folds) .ptpm_balanced_group_folds(groups, k) else
      caret::groupKFold(groups, k = k)
    names(fold_ids) <- paste0("Fold", seq_along(fold_ids), "_Rep", r)
    folds <- c(folds, fold_ids)
  }
  fold_rep <- as.integer(sub(".*_Rep", "", names(folds)))
  cat("Generated", length(folds), "outer train/test splits",
      if (balanced_folds) "(balanced patient folds).\n" else "(caret::groupKFold).\n")
  if (length(folds) < k * repeats) {
    warning("Only ", length(folds), " of ", k * repeats, " outer splits were generated: ",
            "some folds were empty. Consider balanced_folds = TRUE.", call. = FALSE)
  }

  # 15. INNER CV (GROUPED BY PATIENT, WITHIN EACH OUTER TRAINING SET) ----
  # All models of a patient share the same inner fold id, so hyperparameter
  # tuning is also evaluated on unseen patients.
  inner_foldid_list <- vector("list", length(folds))
  for (i in seq_along(folds)) {
    g_train <- groups[folds[[i]]]
    u <- unique(g_train)
    n_f <- min(n_inner, length(u))
    if (n_f < 3) stop("Outer training fold ", i, " has fewer than 3 patients for inner CV.")
    set.seed(seed + 10000 + i)
    patient_fold <- sample(rep(seq_len(n_f), length.out = length(u)))
    names(patient_fold) <- u
    inner_foldid_list[[i]] <- unname(patient_fold[g_train])
  }

  # 16. GLOBAL (PRE-CV) FEATURE RANKING ---------------------------------
  # Used only by the two regimes with FS before CV (full dataset, all patients).
  global_vars <- apply(x_baseline_matrix, 2, stats::var)
  global_gene_order <- names(sort(global_vars, decreasing = TRUE))

  # 17. CORE FOLD FIT ----------------------------------------------------
  # Environment set to globalenv so that workers do not receive a copy of
  # the whole run_ptpm frame (expression matrix) when the function is exported.
  fit_fold <- function(x_train, x_test, y_train, y_test, b_train, b_test,
                       fs_inside, global_gene_order, coupling_correct,
                       top_k_values, alpha_grid, inner_foldid) {
    if (fs_inside) {
      v <- apply(x_train, 2, stats::var)
      gene_order <- names(sort(v, decreasing = TRUE))
    } else {
      gene_order <- global_gene_order
    }
    # OLS coefficients of y on [1, B]; B is a vector (mean baseline score)
    # or a matrix (one column per baseline pathway score)
    ols <- function(B1, y) {
      cf <- stats::lm.fit(B1, y)$coefficients
      cf[is.na(cf)] <- 0
      cf
    }
    if (coupling_correct) {
      # Outer level: beta estimated on the outer training fold only and
      # applied unchanged to the held-out fold
      Btr <- cbind(1, as.matrix(b_train))
      Bte <- cbind(1, as.matrix(b_test))
      cf <- ols(Btr, y_train)
      r_train <- as.numeric(y_train - Btr %*% cf)
      r_test <- as.numeric(y_test - Bte %*% cf)
    }
    pred <- matrix(NA_real_, nrow(x_test), length(top_k_values))
    for (j in seq_along(top_k_values)) {
      top <- gene_order[seq_len(min(top_k_values[j], length(gene_order)))]
      xtr <- x_train[, top, drop = FALSE]
      xte <- x_test[, top, drop = FALSE]
      if (!coupling_correct) {
        fits <- lapply(alpha_grid, function(a) {
          glmnet::cv.glmnet(xtr, y_train, alpha = a, foldid = inner_foldid)
        })
        best <- fits[[which.min(vapply(fits, function(f) min(f$cvm), numeric(1)))]]
        pred[, j] <- as.numeric(stats::predict(best, newx = xte, s = "lambda.min"))
      } else {
        # Inner level: for each inner fold, beta is RE-ESTIMATED on the inner
        # training part only, and both inner training and inner validation
        # outcomes are residualised with it. The lambda path is taken from
        # the fit on the outer-training residuals (as in cv.glmnet), and the
        # CV error is the mean squared error over all inner-validation models.
        best <- list(cvm = Inf)
        for (a in alpha_grid) {
          full <- glmnet::glmnet(xtr, r_train, alpha = a)
          lam <- full$lambda
          err <- matrix(NA_real_, length(y_train), length(lam))
          for (f in unique(inner_foldid)) {
            it <- inner_foldid != f
            iv <- !it
            cf_f <- ols(Btr[it, , drop = FALSE], y_train[it])
            r_it <- as.numeric(y_train[it] - Btr[it, , drop = FALSE] %*% cf_f)
            r_iv <- as.numeric(y_train[iv] - Btr[iv, , drop = FALSE] %*% cf_f)
            fit_f <- glmnet::glmnet(xtr[it, , drop = FALSE], r_it, alpha = a, lambda = lam)
            pr <- stats::predict(fit_f, newx = xtr[iv, , drop = FALSE])
            err[iv, seq_len(ncol(pr))] <- (r_iv - pr)^2
          }
          cvm <- colMeans(err)
          cvm[is.na(cvm)] <- Inf
          i_min <- which.min(cvm)
          if (cvm[i_min] < best$cvm) {
            best <- list(cvm = cvm[i_min], fit = full, lambda = lam[i_min])
          }
        }
        pred[, j] <- as.numeric(stats::predict(best$fit, newx = xte, s = best$lambda))
      }
    }
    list(pred = pred, obs = unname(if (coupling_correct) r_test else y_test))
  }
  environment(fit_fold) <- globalenv()

  # 18. METRICS ----------------------------------------------------------
  r2_sse <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
  r2_cor <- function(obs, pred) {
    if (stats::sd(pred) == 0) return(0)
    stats::cor(obs, pred)^2
  }
  # Per-repeat metrics from a list of fold results (each with $fold, $pred, $obs)
  per_repeat_metrics <- function(fold_results) {
    f_idx <- vapply(fold_results, `[[`, numeric(1), "fold")
    reps <- fold_rep[f_idx]
    out <- list()
    for (r in sort(unique(reps))) {
      sub <- fold_results[reps == r]
      obs <- unlist(lapply(sub, `[[`, "obs"))
      pred <- do.call(rbind, lapply(sub, `[[`, "pred"))
      for (j in seq_along(top_k_values)) {
        out[[length(out) + 1L]] <- data.frame(
          repeat_id = r, K = top_k_values[j],
          R2 = r2_sse(obs, pred[, j]), cor2 = r2_cor(obs, pred[, j]))
      }
    }
    do.call(rbind, out)
  }
  pooled_cor2 <- function(fold_results) {
    obs <- unlist(lapply(fold_results, `[[`, "obs"))
    pred <- do.call(rbind, lapply(fold_results, `[[`, "pred"))
    stats::setNames(vapply(seq_along(top_k_values),
                           function(j) r2_cor(obs, pred[, j]), numeric(1)),
                    paste0("K_", top_k_values))
  }
  summarise_vec <- function(v) {
    c(mean = mean(v), sd = stats::sd(v),
      q025 = unname(stats::quantile(v, 0.025)), q975 = unname(stats::quantile(v, 0.975)),
      prop_positive = mean(v > 0))
  }

  # 19. PARALLEL BACKEND -------------------------------------------------
  cat("Starting parallel backend with", n_cores, "cores...\n")
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Evaluate one regime on a set of outer folds.
  # yb: optional list (per fold, named by fold index) with fold-specific
  # y and b for ALL samples (GSVA sensitivity: training-reference GSVA, so
  # training and test samples share the same training-derived transformation).
  eval_regime <- function(fs_inside, coupling_correct, fold_idx,
                          y_vec = y_main, b_vec = b_main, yb = NULL,
                          x_mat = x_baseline_matrix, gorder = global_gene_order) {
    # local copies: foreach auto-exports only objects of the calling frame
    x_full <- x_mat
    folds <- folds
    fit_fold <- fit_fold
    global_gene_order <- gorder
    top_k_values <- top_k_values
    alpha_grid <- alpha_grid
    inner_foldid_list <- inner_foldid_list
    foreach::foreach(i = fold_idx, .packages = "glmnet") %dopar% {
      tr <- folds[[i]]
      te <- setdiff(seq_len(nrow(x_full)), tr)
      if (is.null(yb)) {
        yv <- y_vec; bv <- b_vec
      } else {
        yv <- yb[[as.character(i)]]$y; bv <- yb[[as.character(i)]]$b
      }
      b_tr <- if (is.matrix(bv)) bv[tr, , drop = FALSE] else bv[tr]
      b_te <- if (is.matrix(bv)) bv[te, , drop = FALSE] else bv[te]
      # fold-specific predictors (training-reference GSVA sensitivity)
      xf <- if (!is.null(yb) && !is.null(yb[[as.character(i)]]$x)) yb[[as.character(i)]]$x else x_full
      res <- fit_fold(xf[tr, , drop = FALSE], xf[te, , drop = FALSE],
                      yv[tr], yv[te], b_tr, b_te,
                      fs_inside, global_gene_order, coupling_correct,
                      top_k_values, alpha_grid, inner_foldid_list[[i]])
      list(fold = i, pred = res$pred, obs = res$obs)
    }
  }

  # 20. FOUR REGIMES (2 x 2 FACTORIAL, GLOBAL GSVA) ----------------------
  regimes <- data.frame(
    key   = c("standard", "coupling_adjusted", "fold_internal", "fully_controlled"),
    label = c("Standard R\u00B2",
              "Coupling-adjusted R\u00B2",
              "Fold-internal R\u00B2",
              "Fully controlled R\u00B2"),
    fs_inside        = c(FALSE, FALSE, TRUE, TRUE),
    coupling_correct = c(FALSE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  all_folds <- seq_along(folds)
  fold_results <- list()
  per_repeat <- list()
  for (q in seq_len(nrow(regimes))) {
    cat("\n-- Regime:", regimes$key[q],
        "| FS", if (regimes$fs_inside[q]) "inside folds" else "pre-CV",
        "| coupling correction", if (regimes$coupling_correct[q]) "yes" else "no", "\n")
    fr <- eval_regime(regimes$fs_inside[q], regimes$coupling_correct[q], all_folds)
    fold_results[[regimes$key[q]]] <- fr
    pr <- per_repeat_metrics(fr)
    pr$regime <- regimes$key[q]
    per_repeat[[regimes$key[q]]] <- pr
    agg <- stats::aggregate(R2 ~ K, pr, mean)
    print(stats::setNames(round(agg$R2, 4), paste0("K_", agg$K)))
  }
  per_repeat_all <- do.call(rbind, per_repeat)
  rownames(per_repeat_all) <- NULL

  # 21. SUMMARY TABLES ---------------------------------------------------
  mean_r2 <- function(key, metric = "R2") {
    pr <- per_repeat[[key]]
    vapply(top_k_values, function(K) mean(pr[[metric]][pr$K == K]), numeric(1))
  }
  results_genome_wide <- data.frame(
    Model = paste0("pathways (K=", top_k_values,
                   ifelse(top_k_values == n_pred, ", all", ""), ")"),
    a = mean_r2("standard"), b = mean_r2("coupling_adjusted"),
    c = mean_r2("fold_internal"), d = mean_r2("fully_controlled"),
    row.names = paste0("K_", top_k_values), stringsAsFactors = FALSE
  )
  colnames(results_genome_wide) <- c("Model", regimes$label)

  results_summary <- do.call(rbind, lapply(regimes$key, function(key) {
    pr <- per_repeat[[key]]
    do.call(rbind, lapply(top_k_values, function(K) {
      s_r2 <- summarise_vec(pr$R2[pr$K == K])
      s_c2 <- summarise_vec(pr$cor2[pr$K == K])
      data.frame(regime = key, K = K,
                 R2_mean = s_r2["mean"], R2_sd = s_r2["sd"],
                 R2_q025 = s_r2["q025"], R2_q975 = s_r2["q975"],
                 cor2_mean = s_c2["mean"], cor2_sd = s_c2["sd"],
                 cor2_pooled_legacy = pooled_cor2(fold_results[[key]])[paste0("K_", K)],
                 row.names = NULL)
    }))
  }))

  # 22. PAIRED DIFFERENCES BETWEEN REGIMES (PER REPEAT) -----------------
  wide <- function(metric) {
    w <- stats::reshape(per_repeat_all[, c("repeat_id", "K", "regime", metric)],
                        idvar = c("repeat_id", "K"), timevar = "regime",
                        direction = "wide")
    names(w) <- sub(paste0("^", metric, "\\."), "", names(w))
    w
  }
  contrasts <- list(
    Delta_FS_no_coupling_adj    = function(w) w$standard - w$fold_internal,
    Delta_FS_with_coupling_adj  = function(w) w$coupling_adjusted - w$fully_controlled,
    Delta_coupling_FS_preCV     = function(w) w$standard - w$coupling_adjusted,
    Delta_coupling_FS_inside    = function(w) w$fold_internal - w$fully_controlled,
    Interaction                 = function(w) (w$standard - w$fold_internal) - (w$coupling_adjusted - w$fully_controlled),
    Delta_FS_shapley            = function(w) ((w$standard - w$fold_internal) + (w$coupling_adjusted - w$fully_controlled)) / 2,
    Delta_coupling_shapley      = function(w) ((w$standard - w$coupling_adjusted) + (w$fold_internal - w$fully_controlled)) / 2,
    Total_gap                   = function(w) w$standard - w$fully_controlled
  )
  paired_differences <- do.call(rbind, lapply(c("R2", "cor2"), function(metric) {
    w <- wide(metric)
    do.call(rbind, lapply(names(contrasts), function(nm) {
      d <- contrasts[[nm]](w)
      do.call(rbind, lapply(top_k_values, function(K) {
        s <- summarise_vec(d[w$K == K])
        data.frame(metric = metric, contrast = nm, K = K, t(s), row.names = NULL)
      }))
    }))
  }))
  paired_per_repeat <- wide("R2")

  # 22b. CONFIDENCE INTERVALS ---------------------------------------------
  alpha_ci <- (1 - conf_level) / 2
  nK <- length(top_k_values)
  keys <- regimes$key

  # (a) Patient-level cluster bootstrap of out-of-fold predictions
  oof_arrays <- function(fr, reps = seq_len(repeats)) {
    O <- matrix(NA_real_, n_obs, length(reps))
    P <- array(NA_real_, c(n_obs, length(reps), nK))
    for (z in fr) {
      i <- z$fold
      col <- match(fold_rep[i], reps)
      if (is.na(col)) next
      te <- setdiff(seq_len(n_obs), folds[[i]])
      O[te, col] <- z$obs
      P[te, col, ] <- z$pred
    }
    if (anyNA(O) || anyNA(P)) stop("Incomplete out-of-fold predictions.")
    list(O = O, P = P)
  }
  # Bootstrap distribution (n_boot x nK) of the mean-over-repeats R2
  boot_mean_r2 <- function(fr, reps = seq_len(repeats)) {
    arr <- oof_arrays(fr, reps)
    O <- arr$O
    mo <- crossprod(W, O) / sw              # n_boot x repeats (weighted means)
    sst <- crossprod(W, O^2) - sw * mo^2
    out <- vapply(seq_len(nK), function(j) {
      Pj <- matrix(arr$P[, , j], nrow = n_obs)
      sse <- crossprod(W, (O - Pj)^2)
      rowMeans(1 - sse / sst)               # mean over repeats, per resample
    }, numeric(n_boot))
    matrix(out, nrow = n_boot)
  }
  boot_R2 <- NULL
  if (n_boot > 0) {
    cat("\n-- Patient-level bootstrap:", n_boot, "resamples...\n")
    patients <- unique(groups)
    pat_idx <- match(groups, patients)
    set.seed(seed + 70000)
    W <- vapply(seq_len(n_boot), function(b) {
      cnt <- tabulate(sample(length(patients), replace = TRUE), nbins = length(patients))
      as.numeric(cnt[pat_idx])
    }, numeric(n_obs))                       # n_obs x n_boot matrix of weights
    sw <- colSums(W)
    boot_R2 <- lapply(stats::setNames(keys, keys), function(key) {
      boot_mean_r2(fold_results[[key]])
    })
  }

  # (b) Corrected resampled t (Nadeau & Bengio 2003) on per-fold R2
  per_fold_r2 <- lapply(stats::setNames(keys, keys), function(key) {
    t(vapply(fold_results[[key]], function(z) {
      vapply(seq_len(nK), function(j) r2_sse(z$obs, z$pred[, j]), numeric(1))
    }, numeric(nK)))                        # folds x nK (rows in fold order)
  })
  if (nK == 1L) per_fold_r2 <- lapply(per_fold_r2, function(m) matrix(m, ncol = 1L))
  ratio_te_tr <- mean(vapply(folds, function(tr) (n_obs - length(tr)) / length(tr), numeric(1)))
  nb_interval <- function(d) {
    J <- length(d); m <- mean(d)
    se <- sqrt((1 / J + ratio_te_tr) * stats::var(d))
    tq <- stats::qt(1 - alpha_ci, J - 1)
    c(nb_mean = m, nb_lower = m - tq * se, nb_upper = m + tq * se,
      nb_p = if (se > 0) 2 * stats::pt(-abs(m / se), J - 1) else NA_real_)
  }
  boot_q <- function(v) {
    c(boot_lower = unname(stats::quantile(v, alpha_ci)),
      boot_upper = unname(stats::quantile(v, 1 - alpha_ci)),
      boot_se = stats::sd(v))
  }

  confidence_intervals <- do.call(rbind, lapply(keys, function(key) {
    do.call(rbind, lapply(seq_len(nK), function(j) {
      est <- mean_r2(key)[j]
      bq <- if (!is.null(boot_R2)) boot_q(boot_R2[[key]][, j]) else
        c(boot_lower = NA, boot_upper = NA, boot_se = NA)
      nb <- nb_interval(per_fold_r2[[key]][, j])
      data.frame(regime = key, K = top_k_values[j], R2 = est, t(bq),
                 t(nb[c("nb_mean", "nb_lower", "nb_upper")]), row.names = NULL)
    }))
  }))

  contrast_intervals <- do.call(rbind, lapply(names(contrasts), function(nm) {
    do.call(rbind, lapply(seq_len(nK), function(j) {
      est_list <- lapply(stats::setNames(keys, keys), function(key) mean_r2(key)[j])
      est <- contrasts[[nm]](est_list)
      if (!is.null(boot_R2)) {
        bl <- lapply(stats::setNames(keys, keys), function(key) boot_R2[[key]][, j])
        bd <- contrasts[[nm]](bl)
        bq <- c(boot_q(bd), boot_prop_positive = mean(bd > 0))
      } else {
        bq <- c(boot_lower = NA, boot_upper = NA, boot_se = NA, boot_prop_positive = NA)
      }
      fl <- lapply(stats::setNames(keys, keys), function(key) per_fold_r2[[key]][, j])
      nb <- nb_interval(contrasts[[nm]](fl))
      data.frame(contrast = nm, K = top_k_values[j], estimate = est,
                 t(bq), t(nb), row.names = NULL)
    }))
  }))

  # 23. GSVA SENSITIVITY: TRAINING-REFERENCE GSVA -------------------------
  # Main analysis: GSVA on all samples (cohort-level, unsupervised).
  # Sensitivity: for each outer fold, the only cross-sample step of GSVA
  # (per-gene CDF) is estimated on the NT + CET profiles of the TRAINING
  # models and then applied to ALL samples; the per-sample steps (ranking,
  # random walk) are unchanged. Training and test samples are therefore
  # scored with the same, training-derived transformation, and held-out
  # patients do not contribute to it.
  # Validity check: with reference = all samples, the implementation must
  # reproduce GSVA::gsva() on these data (max |difference| < gsva_ref_tol);
  # the kernel type (Gaussian vs direct ECDF, which GSVA >= 2.0 selects
  # automatically from the number of samples) and tie handling are taken
  # from the variant that reproduces the global scores. If none does, the
  # sensitivity analysis is skipped with a warning.
  gsva_sens <- NULL
  per_repeat_gsva <- NULL
  gsva_ref_validation <- NULL
  gsva_score_agreement <- NULL
  gsva_ref_settings <- NULL
  if (gsva_sensitivity) {
    cat("\n-- GSVA sensitivity: training-reference GSVA\n")
    cat("   validating against GSVA::gsva() (reference = all samples)...\n")
    all_cols <- seq_len(ncol(expr_human))
    val_rows <- list()
    for (method in c("gaussian", "none")) {
      base <- .ptpm_gsva_ref(expr_human, all_sets, all_cols, method, "last", cl = cl)
      for (ties in c("last", "first")) {
        sc <- if (ties == "last") base$scores else
          .ptpm_gsva_ref(expr_human, all_sets, all_cols, method, ties, cdf = base$cdf)$scores
        d <- if (setequal(rownames(sc), rownames(gsva_all)) &&
                 setequal(colnames(sc), colnames(gsva_all))) {
          max(abs(sc[rownames(gsva_all), colnames(gsva_all)] - gsva_all))
        } else Inf
        val_rows[[length(val_rows) + 1L]] <- data.frame(kcdf = method, ties = ties,
                                                        max_abs_diff = d)
      }
    }
    gsva_ref_validation <- do.call(rbind, val_rows)
    ok <- which(gsva_ref_validation$max_abs_diff < gsva_ref_tol)
    print(gsva_ref_validation, digits = 3, row.names = FALSE)

    if (length(ok) == 0L) {
      warning("Reference-based GSVA does not reproduce GSVA::gsva() on these data ",
              "(installed GSVA ", as.character(utils::packageVersion("GSVA")),
              "); GSVA sensitivity analysis skipped.")
    } else {
      gsva_ref_settings <- gsva_ref_validation[ok[1], ]
      method <- gsva_ref_settings$kcdf
      ties <- gsva_ref_settings$ties
      cat("   OK: kcdf =", method, "| ties =", ties,
          "| max |diff| =", signif(gsva_ref_settings$max_abs_diff, 3), "\n")

      sens_folds <- which(fold_rep <= gsva_sensitivity_repeats)
      cat("   scoring", length(sens_folds), "folds with training-reference GSVA...\n")
      yb <- list()
      agree <- list()
      for (i in sens_folds) {
        tr <- folds[[i]]
        te <- setdiff(seq_len(n_obs), tr)
        tr_pairs <- rownames(x_baseline_matrix)[tr]
        ref_cols <- match(c(paste0(tr_pairs, "_NT"), paste0(tr_pairs, "_CET")),
                          colnames(expr_human))
        sc_all <- .ptpm_gsva_ref(expr_human, all_sets, ref_cols, method, ties, cl = cl)$scores
        sc <- sc_all[outcome_rows, , drop = FALSE]
        sp <- score_pairs(sc, pair_ids)
        y_f <- unname(sp$S[rownames(x_baseline_matrix)])
        b_f <- unname(sp$b[rownames(x_baseline_matrix)])
        x_f <- t(sc_all[colnames(x_baseline_matrix),
                        paste0(rownames(x_baseline_matrix), "_NT"), drop = FALSE])
        rownames(x_f) <- rownames(x_baseline_matrix)
        yb[[as.character(i)]] <- list(y = y_f, b = b_f, x = x_f)
        agree[[length(agree) + 1L]] <- data.frame(
          fold = i, repeat_id = fold_rep[i],
          cor_S_train = stats::cor(y_f[tr], y_main[tr]),
          cor_S_test = stats::cor(y_f[te], y_main[te]),
          cor_b_test = stats::cor(b_f[te], b_main[te]),
          shift_S_test = mean(y_f[te] - y_main[te]),
          sd_ratio_S_test = stats::sd(y_f[te]) / stats::sd(y_main[te]),
          cor_X_test = stats::cor(as.vector(x_f[te, ]), as.vector(x_baseline_matrix[te, ])))
      }
      gsva_score_agreement <- do.call(rbind, agree)

      sens_keys <- c("fold_internal", "fully_controlled")
      fr_gsva <- list()
      per_repeat_gsva <- do.call(rbind, lapply(sens_keys, function(key) {
        q <- which(regimes$key == key)
        cat("   regime:", key, "\n")
        fr <- eval_regime(regimes$fs_inside[q], regimes$coupling_correct[q],
                          sens_folds, yb = yb)
        fr_gsva[[key]] <<- fr
        pr <- per_repeat_metrics(fr)
        pr$regime <- key
        pr
      }))

      sens_reps <- seq_len(gsva_sensitivity_repeats)
      gsva_sens <- do.call(rbind, lapply(sens_keys, function(key) {
        bg <- bf <- NULL
        if (n_boot > 0) {
          bg <- boot_mean_r2(fold_results[[key]], sens_reps)
          bf <- boot_mean_r2(fr_gsva[[key]], sens_reps)
        }
        do.call(rbind, lapply(seq_len(nK), function(j) {
          K <- top_k_values[j]
          g <- per_repeat[[key]]
          g <- g[g$K == K & g$repeat_id <= gsva_sensitivity_repeats, ]
          f <- per_repeat_gsva[per_repeat_gsva$regime == key & per_repeat_gsva$K == K, ]
          f <- f[match(g$repeat_id, f$repeat_id), ]
          s <- summarise_vec(g$R2 - f$R2)
          bq <- if (!is.null(bg)) boot_q(bg[, j] - bf[, j]) else
            c(boot_lower = NA, boot_upper = NA, boot_se = NA)
          data.frame(regime = key, K = K,
                     R2_global_gsva = mean(g$R2), R2_train_ref_gsva = mean(f$R2),
                     cor2_global_gsva = mean(g$cor2), cor2_train_ref_gsva = mean(f$cor2),
                     diff_R2_mean = s["mean"], diff_R2_sd_repeats = s["sd"],
                     diff_boot_lower = bq["boot_lower"], diff_boot_upper = bq["boot_upper"],
                     prop_repeats_global_higher = s["prop_positive"], row.names = NULL)
        }))
      }))
    }
  }

  # 23b. ROBUSTNESS: KRAS CIRCULARITY AND PATHWAY-SPECIFIC ADJUSTMENT -----
  # (1) Pathway-specific coupling adjustment: S is residualised (within each
  #     outer training fold) on the baseline score of EACH pathway instead of
  #     their mean.
  # (2) KRAS-free phenotype: S and its baseline score are rebuilt from the
  #     non-KRAS pathways only, so the KRAS pathway is no longer part of the
  #     outcome.
  # (3) KRAS genes excluded from the baseline predictors (FS inside folds).
  # (4) Partial correlation of the baseline KRAS pathway score with the full
  #     and the KRAS-free phenotype, adjusted for the phenotype's own baseline
  #     score (in-sample, descriptive; patient-level bootstrap CI). The
  #     version adjusted only for the other pathways keeps the mechanical
  #     KRAS term shared by S and the predictor and serves as a coupled
  #     reference; the KRAS-free version is the circularity test.
  # All models use fold-internal feature selection and the same folds as the
  # main analysis; differences vs the main regimes are paired per repeat.
  robustness_r2 <- NULL
  robustness_partial_cor <- NULL
  rob_fr <- list()
  rob_pr <- list()
  if (robustness) {
    if (!kras_set %in% rownames(gsva_scores)) {
      warning("Gene set '", kras_set, "' not found among the GSVA scores; ",
              "robustness analyses skipped.")
    } else if (nrow(gsva_scores) < 2L) {
      warning("At least two gene sets are needed for the KRAS-free phenotype; ",
              "robustness analyses skipped.")
    } else {
      cat("\n-- Robustness analyses (pathway-specific adjustment, KRAS circularity)\n")
      xr <- rownames(x_baseline_matrix)
      other_sets <- setdiff(rownames(gsva_scores), kras_set)
      B_path <- t(gsva_scores[, paste0(xr, "_NT"), drop = FALSE])     # models x pathways
      rownames(B_path) <- xr
      b_kras <- unname(B_path[, kras_set])
      nk <- score_pairs(gsva_scores, pair_ids, sets = other_sets)
      y_nok <- unname(nk$S[xr])
      b_nok <- unname(nk$b[xr])
      keep_ov <- pathway_overlap$pathway[pathway_overlap$frac_shared_with_outcome <= overlap_threshold]
      x_nok <- x_baseline_matrix[, intersect(colnames(x_baseline_matrix), keep_ov), drop = FALSE]
      cat("   KRAS set:", kras_set, "| other sets:", paste(other_sets, collapse = ", "),
          "| predictor pathways kept after overlap filter (<=", overlap_threshold, "):",
          ncol(x_nok), "of", ncol(x_baseline_matrix), "\n")

      analyses <- list(
        list(id = "pathway_adjusted", ref = "fully_controlled", cc = TRUE,
             y = y_main, b = unname(B_path), x = x_baseline_matrix,
             desc = "Full S; adjustment on each baseline pathway score"),
        list(id = "kras_free_S_fold_internal", ref = NA_character_, cc = FALSE,
             y = y_nok, b = b_nok, x = x_baseline_matrix,
             desc = "KRAS-free S; no coupling adjustment"),
        list(id = "kras_free_S_fully_controlled", ref = NA_character_, cc = TRUE,
             y = y_nok, b = b_nok, x = x_baseline_matrix,
             desc = "KRAS-free S; adjustment on KRAS-free baseline score"),
        list(id = "overlap_excluded_fold_internal", ref = "fold_internal", cc = FALSE,
             y = y_main, b = b_main, x = x_nok,
             desc = "Full S; predictor pathways with high gene overlap with outcome sets excluded"),
        list(id = "overlap_excluded_fully_controlled", ref = "fully_controlled", cc = TRUE,
             y = y_main, b = b_main, x = x_nok,
             desc = "Full S; high-overlap predictor pathways excluded; mean-score adjustment")
      )

      rob_folds <- which(fold_rep <= robustness_repeats)
      rob_reps <- seq_len(robustness_repeats)
      ref_boot <- list()
      robustness_r2 <- do.call(rbind, lapply(analyses, function(a) {
        cat("   analysis:", a$id, "\n")
        fr <- eval_regime(TRUE, a$cc, rob_folds, y_vec = a$y, b_vec = a$b, x_mat = a$x)
        pr <- per_repeat_metrics(fr)
        rob_fr[[a$id]] <<- fr
        rob_pr[[a$id]] <<- pr
        bR <- if (n_boot > 0) boot_mean_r2(fr, rob_reps) else NULL
        if (!is.na(a$ref) && n_boot > 0 && is.null(ref_boot[[a$ref]])) {
          ref_boot[[a$ref]] <<- boot_mean_r2(fold_results[[a$ref]], rob_reps)
        }
        do.call(rbind, lapply(seq_len(nK), function(j) {
          K <- top_k_values[j]
          r2 <- pr$R2[pr$K == K]
          ci <- if (!is.null(bR)) boot_q(bR[, j]) else c(boot_lower = NA, boot_upper = NA, boot_se = NA)
          out <- data.frame(analysis = a$id, description = a$desc, K = K,
                            R2 = mean(r2), R2_boot_lower = ci["boot_lower"],
                            R2_boot_upper = ci["boot_upper"],
                            reference = a$ref, R2_reference = NA_real_,
                            diff_vs_reference = NA_real_,
                            diff_boot_lower = NA_real_, diff_boot_upper = NA_real_,
                            row.names = NULL)
          if (!is.na(a$ref)) {
            g <- per_repeat[[a$ref]]
            g <- g[g$K == K & g$repeat_id %in% rob_reps, ]
            f <- pr[pr$K == K, ]
            f <- f[match(g$repeat_id, f$repeat_id), ]
            out$R2_reference <- mean(g$R2)
            out$diff_vs_reference <- mean(f$R2 - g$R2)
            if (n_boot > 0) {
              dq <- boot_q(bR[, j] - ref_boot[[a$ref]][, j])
              out$diff_boot_lower <- dq["boot_lower"]
              out$diff_boot_upper <- dq["boot_upper"]
            }
          }
          out
        }))
      }))

      # Partial correlations with patient-level bootstrap
      pcor <- function(y, x, z) {
        z1 <- cbind(1, z)
        ry <- stats::lm.fit(z1, y)$residuals
        rx <- stats::lm.fit(z1, x)$residuals
        stats::cor(ry, rx)
      }
      pc_specs <- list(
        list(id = "KRAS_vs_full_S_adj_mean_baseline",
             desc = "baseline KRAS score vs full S | mean baseline score",
             y = y_main, z = b_main),
        list(id = "KRAS_vs_full_S_adj_other_pathways",
             desc = "baseline KRAS score vs full S | other pathways (coupled reference)",
             y = y_main, z = unname(B_path[, other_sets, drop = FALSE])),
        list(id = "KRAS_vs_KRAS_free_S_adj_own_baseline",
             desc = "baseline KRAS score vs KRAS-free S | KRAS-free baseline score",
             y = y_nok, z = b_nok),
        list(id = "KRAS_vs_KRAS_free_S_adj_each_pathway",
             desc = "baseline KRAS score vs KRAS-free S | each non-KRAS baseline pathway",
             y = y_nok, z = unname(B_path[, other_sets, drop = FALSE]))
      )
      patients <- unique(groups)
      rows_by_patient <- split(seq_len(n_obs), factor(groups, levels = patients))
      set.seed(seed + 80000)
      n_pb <- max(n_boot, 1000)
      boot_rows <- replicate(n_pb, unlist(rows_by_patient[sample(length(patients),
                                                               replace = TRUE)],
                                          use.names = FALSE), simplify = FALSE)
      robustness_partial_cor <- do.call(rbind, lapply(pc_specs, function(s) {
        z <- as.matrix(s$z)
        est <- pcor(s$y, b_kras, z)
        bt <- vapply(boot_rows, function(idx) pcor(s$y[idx], b_kras[idx],
                                                   z[idx, , drop = FALSE]), numeric(1))
        data.frame(analysis = s$id, description = s$desc, partial_r = est,
                   boot_lower = unname(stats::quantile(bt, alpha_ci, na.rm = TRUE)),
                   boot_upper = unname(stats::quantile(bt, 1 - alpha_ci, na.rm = TRUE)),
                   n_boot = n_pb, row.names = NULL)
      }))
    }
  }

  # 23c. SECONDARY SENSITIVITY: ssGSEA INSTEAD OF GSVA -------------------
  # ssGSEA without normalisation (normalize = FALSE) is a single-sample
  # method: each sample's score depends only on that sample, so no
  # information crosses between training and test patients by construction.
  # S and b are rebuilt from ssGSEA scores and the four main regimes (plus
  # the pathway-specific adjustment) are re-run on the same folds. R2 values
  # are compared with the GSVA-based results (paired per repeat, patient
  # bootstrap intervals). The ssGSEA scale differs from GSVA, so the
  # comparison is on R2, not on the scores themselves.
  ssgsea_results <- NULL
  ssgsea_contrasts <- NULL
  ssgsea_agreement <- NULL
  ssgsea_coupling_null <- NULL
  ssgsea_partial_cor <- NULL
  if (ssgsea_sensitivity) {
    cat("\n-- Secondary sensitivity: ssGSEA (normalize = FALSE) instead of GSVA\n")
    ss_param <- GSVA::ssgseaParam(expr_human, all_sets, alpha = ssgsea_alpha,
                                  normalize = FALSE)
    ss_all <- suppressMessages(GSVA::gsva(ss_param))
    ss_scores <- ss_all[outcome_rows, , drop = FALSE]
    xr <- rownames(x_baseline_matrix)
    X_ss <- t(ss_all[colnames(x_baseline_matrix), paste0(xr, "_NT"), drop = FALSE])
    rownames(X_ss) <- xr
    ss_order <- names(sort(apply(X_ss, 2, stats::var), decreasing = TRUE))
    ss <- score_pairs(ss_scores, pair_ids)
    y_ss <- unname(ss$S[xr])
    b_ss <- unname(ss$b[xr])
    B_ss <- t(ss_scores[, paste0(xr, "_NT"), drop = FALSE])
    ssgsea_agreement <- data.frame(
      cor_S_ssgsea_vs_gsva = stats::cor(y_ss, y_main),
      cor_b_ssgsea_vs_gsva = stats::cor(b_ss, b_main),
      cor_b_S_gsva = stats::cor(b_main, y_main),
      cor_b_S_ssgsea = stats::cor(b_ss, y_ss))
    print(round(ssgsea_agreement, 4), row.names = FALSE)

    ss_folds <- which(fold_rep <= ssgsea_repeats)
    ss_reps <- seq_len(ssgsea_repeats)
    ss_runs <- list(
      list(key = "standard", fs = FALSE, cc = FALSE, b = b_ss, ref_fr = fold_results$standard,
           ref_pr = per_repeat$standard),
      list(key = "coupling_adjusted", fs = FALSE, cc = TRUE, b = b_ss,
           ref_fr = fold_results$coupling_adjusted, ref_pr = per_repeat$coupling_adjusted),
      list(key = "fold_internal", fs = TRUE, cc = FALSE, b = b_ss,
           ref_fr = fold_results$fold_internal, ref_pr = per_repeat$fold_internal),
      list(key = "fully_controlled", fs = TRUE, cc = TRUE, b = b_ss,
           ref_fr = fold_results$fully_controlled, ref_pr = per_repeat$fully_controlled),
      list(key = "pathway_adjusted", fs = TRUE, cc = TRUE, b = unname(B_ss),
           ref_fr = rob_fr$pathway_adjusted, ref_pr = rob_pr$pathway_adjusted)
    )
    ss_boot <- list()
    ss_pr_all <- list()
    ssgsea_results <- do.call(rbind, lapply(ss_runs, function(a) {
      cat("   regime:", a$key, "\n")
      fr <- eval_regime(a$fs, a$cc, ss_folds, y_vec = y_ss, b_vec = a$b,
                        x_mat = X_ss, gorder = ss_order)
      pr <- per_repeat_metrics(fr)
      ss_pr_all[[a$key]] <<- pr
      bS <- if (n_boot > 0) boot_mean_r2(fr, ss_reps) else NULL
      bG <- if (n_boot > 0 && !is.null(a$ref_fr)) boot_mean_r2(a$ref_fr, ss_reps) else NULL
      ss_boot[[a$key]] <<- bS
      do.call(rbind, lapply(seq_len(nK), function(j) {
        K <- top_k_values[j]
        f <- pr[pr$K == K, ]
        ci <- if (!is.null(bS)) boot_q(bS[, j]) else c(boot_lower = NA, boot_upper = NA, boot_se = NA)
        out <- data.frame(regime = a$key, K = K, R2_ssgsea = mean(f$R2),
                          R2_ssgsea_boot_lower = ci["boot_lower"],
                          R2_ssgsea_boot_upper = ci["boot_upper"],
                          R2_gsva = NA_real_, diff_ssgsea_minus_gsva = NA_real_,
                          diff_boot_lower = NA_real_, diff_boot_upper = NA_real_,
                          row.names = NULL)
        if (!is.null(a$ref_pr)) {
          g <- a$ref_pr[a$ref_pr$K == K & a$ref_pr$repeat_id %in% ss_reps, ]
          f <- f[match(g$repeat_id, f$repeat_id), ]
          out$R2_gsva <- mean(g$R2)
          out$diff_ssgsea_minus_gsva <- mean(f$R2 - g$R2)
          if (!is.null(bS) && !is.null(bG)) {
            dq <- boot_q(bS[, j] - bG[, j])
            out$diff_boot_lower <- dq["boot_lower"]
            out$diff_boot_upper <- dq["boot_upper"]
          }
        }
        out
      }))
    }))

    # Coupling null with ssGSEA scores: observed cor(b, S) vs CET scores
    # shuffled across models (arithmetic coupling only)
    c_ss <- unname(ss$c[xr])
    set.seed(seed + 90000)
    r_null_ss <- if (n_perm_coupling > 0) {
      replicate(n_perm_coupling, stats::cor(b_ss, b_ss - c_ss[sample(n_obs)]))
    } else numeric(0)
    ssgsea_coupling_null <- data.frame(
      r_observed = stats::cor(b_ss, y_ss),
      r_null_mean = if (length(r_null_ss)) mean(r_null_ss) else NA_real_,
      r_null_q025 = if (length(r_null_ss)) unname(stats::quantile(r_null_ss, 0.025)) else NA_real_,
      r_null_q975 = if (length(r_null_ss)) unname(stats::quantile(r_null_ss, 0.975)) else NA_real_,
      n_perm = n_perm_coupling)

    # KRAS partial correlations with ssGSEA scores (same specifications as
    # the GSVA-based robustness analysis)
    if (kras_set %in% rownames(ss_scores) && nrow(ss_scores) >= 2L) {
      other_ss <- setdiff(rownames(ss_scores), kras_set)
      rownames(B_ss) <- xr
      bk_ss <- unname(B_ss[, kras_set])
      nk_ss <- score_pairs(ss_scores, pair_ids, sets = other_ss)
      y_nok_ss <- unname(nk_ss$S[xr])
      b_nok_ss <- unname(nk_ss$b[xr])
      pcor_ss <- function(y, x, z) {
        z1 <- cbind(1, z)
        stats::cor(stats::lm.fit(z1, y)$residuals, stats::lm.fit(z1, x)$residuals)
      }
      specs_ss <- list(
        list(desc = "baseline KRAS score vs full S | mean baseline score",
             y = y_ss, z = b_ss),
        list(desc = "baseline KRAS score vs full S | other pathways (coupled reference)",
             y = y_ss, z = unname(B_ss[, other_ss, drop = FALSE])),
        list(desc = "baseline KRAS score vs KRAS-free S | KRAS-free baseline score",
             y = y_nok_ss, z = b_nok_ss),
        list(desc = "baseline KRAS score vs KRAS-free S | each non-KRAS baseline pathway",
             y = y_nok_ss, z = unname(B_ss[, other_ss, drop = FALSE]))
      )
      pats <- unique(groups)
      rows_by_pat <- split(seq_len(n_obs), factor(groups, levels = pats))
      set.seed(seed + 95000)
      n_pb <- max(n_boot, 1000)
      bt_rows <- replicate(n_pb, unlist(rows_by_pat[sample(length(pats), replace = TRUE)],
                                        use.names = FALSE), simplify = FALSE)
      ssgsea_partial_cor <- do.call(rbind, lapply(specs_ss, function(sp) {
        z <- as.matrix(sp$z)
        bt <- vapply(bt_rows, function(idx) pcor_ss(sp$y[idx], bk_ss[idx],
                                                    z[idx, , drop = FALSE]), numeric(1))
        data.frame(description = sp$desc, partial_r = pcor_ss(sp$y, bk_ss, z),
                   boot_lower = unname(stats::quantile(bt, alpha_ci, na.rm = TRUE)),
                   boot_upper = unname(stats::quantile(bt, 1 - alpha_ci, na.rm = TRUE)),
                   n_boot = n_pb, row.names = NULL)
      }))
    }

    # Factorial decomposition with ssGSEA (same contrasts as the main analysis)
    ssgsea_contrasts <- do.call(rbind, lapply(names(contrasts), function(nm) {
      do.call(rbind, lapply(seq_len(nK), function(j) {
        K <- top_k_values[j]
        est_list <- lapply(stats::setNames(keys, keys), function(key) {
          pr <- ss_pr_all[[key]]
          mean(pr$R2[pr$K == K])
        })
        est <- contrasts[[nm]](est_list)
        bq <- c(boot_lower = NA, boot_upper = NA, boot_se = NA)
        if (n_boot > 0) {
          bl <- lapply(stats::setNames(keys, keys), function(key) ss_boot[[key]][, j])
          bq <- boot_q(contrasts[[nm]](bl))
        }
        data.frame(contrast = nm, K = K, estimate = est,
                   boot_lower = bq["boot_lower"], boot_upper = bq["boot_upper"],
                   row.names = NULL)
      }))
    }))
  }

  # 24. PERMUTATION NULLS OF THE NT/CET PAIRING ------------------------------
  set.seed(seed + 50000)
  r_obs <- stats::cor(b_main, y_main)
  r_null <- if (n_perm_coupling > 0) {
    replicate(n_perm_coupling, stats::cor(b_main, b_main - c_main[sample(n_obs)]))
  } else numeric(0)
  coupling_null <- data.frame(
    r_observed = r_obs,
    r_null_mean = if (length(r_null)) mean(r_null) else NA_real_,
    r_null_q025 = if (length(r_null)) unname(stats::quantile(r_null, 0.025)) else NA_real_,
    r_null_q975 = if (length(r_null)) unname(stats::quantile(r_null, 0.975)) else NA_real_,
    n_perm = n_perm_coupling
  )

  permutation_null <- NULL
  null_distribution <- NULL
  if (n_perm > 0) {
    perm_folds <- which(fold_rep <= perm_repeats)
    set.seed(seed + 60000)
    perms <- replicate(n_perm, sample(n_obs))
    tasks <- expand.grid(p = seq_len(n_perm), f = perm_folds)
    cat("\n-- Pairing null:", nrow(tasks), "permutation x fold tasks...\n")
    x_full <- x_baseline_matrix
    null_raw <- foreach::foreach(t = seq_len(nrow(tasks)), .packages = "glmnet") %dopar% {
      p <- tasks$p[t]; i <- tasks$f[t]
      y_perm <- b_main - c_main[perms[, p]]
      tr <- folds[[i]]
      te <- setdiff(seq_len(nrow(x_full)), tr)
      lapply(c(FALSE, TRUE), function(cc) {
        res <- fit_fold(x_full[tr, , drop = FALSE], x_full[te, , drop = FALSE],
                        y_perm[tr], y_perm[te], b_main[tr], b_main[te],
                        TRUE, global_gene_order, cc,
                        top_k_values, alpha_grid, inner_foldid_list[[i]])
        list(p = p, fold = i, coupling = cc, pred = res$pred, obs = res$obs)
      })
    }
    null_raw <- unlist(null_raw, recursive = FALSE)
    null_distribution <- do.call(rbind, lapply(seq_len(n_perm), function(p) {
      do.call(rbind, lapply(c(FALSE, TRUE), function(cc) {
        sel <- Filter(function(z) z$p == p && z$coupling == cc, null_raw)
        pr <- per_repeat_metrics(sel)
        agg <- stats::aggregate(cbind(R2, cor2) ~ K, pr, mean)
        data.frame(perm = p, regime = if (cc) "fully_controlled" else "fold_internal", agg)
      }))
    }))
    permutation_null <- do.call(rbind, lapply(c("fold_internal", "fully_controlled"), function(key) {
      do.call(rbind, lapply(top_k_values, function(K) {
        obs_r2 <- per_repeat[[key]]
        obs_r2 <- mean(obs_r2$R2[obs_r2$K == K & obs_r2$repeat_id <= perm_repeats])
        nd <- null_distribution$R2[null_distribution$regime == key & null_distribution$K == K]
        data.frame(regime = key, K = K, R2_observed = obs_r2,
                   R2_null_mean = mean(nd),
                   R2_null_q025 = unname(stats::quantile(nd, 0.025)),
                   R2_null_q975 = unname(stats::quantile(nd, 0.975)),
                   p_value = (1 + sum(nd >= obs_r2)) / (1 + length(nd)),
                   row.names = NULL)
      }))
    }))
  }

  # 25. LEAKAGE AUDIT ------------------------------------------------------
  leakage_audit <- data.frame(
    step = c("Outer train/test split",
             "Inner (tuning) split",
             "Gene ranking (feature selection)",
             "Elastic-net alpha / lambda",
             "Coupling regression (S ~ b)",
             "GSVA pathway scores (S, b)"),
    Standard = c("patient-grouped", "patient-grouped", "full dataset (pre-CV)",
                     "inner CV only", "not applied", "global (whole cohort)"),
    Coupling_adjusted = c("patient-grouped", "patient-grouped", "full dataset (pre-CV)",
                     "inner CV only", "outer and inner training folds", "global (whole cohort)"),
    Fold_internal = c("patient-grouped", "patient-grouped", "outer training fold",
                     "inner CV only", "not applied", "global (whole cohort)"),
    Fully_controlled = c("patient-grouped", "patient-grouped", "outer training fold",
                         "inner CV only", "outer and inner training folds", "global (whole cohort)"),
    note = c("",
             "inner folds built on patients of the outer training set",
             "variance ranking is unsupervised (does not use S)",
             "",
             "re-estimated within each inner training fold during tuning",
             "unsupervised, transductive; see gsva_sensitivity (training-reference GSVA)"),
    stringsAsFactors = FALSE
  )
  definitions <- c(
    fully_controlled = paste(
      "Every data-dependent step that could carry information about a held-out",
      "patient (outer and inner splits, gene ranking, hyperparameters, coupling",
      "regression) is estimated on training patients only. Exception: GSVA scores",
      "defining S and b are computed on the whole cohort (unsupervised); their",
      "influence is quantified by the training-reference GSVA sensitivity analysis."),
    R2 = paste(
      "Mean over repeats of R2_r = 1 - SSE_r/SST_r, where R2_r is computed on the",
      "pooled out-of-fold predictions of repeat r (not pooled over repeats,",
      "not averaged over folds)."),
    paired_difference = paste(
      "Difference between regimes within the same repeat (same folds);",
      "summaries describe partition variability, not sampling uncertainty."),
    confidence_intervals = paste(
      "Sampling uncertainty: patient-level cluster bootstrap of out-of-fold",
      "predictions (models held fixed) and corrected resampled t interval",
      "(Nadeau & Bengio 2003) on per-fold R2.")
  )

  # 26. PRINT & RETURN -----------------------------------------------------
  cat("\n============================================================\n")
  cat("FINAL RESULTS (R\u00B2 = mean of", repeats, "per-repeat R\u00B2, each 1 - SSE/SST on that repeat's out-of-fold predictions)\n")
  cat("============================================================\n")
  print(results_genome_wide, digits = 4)
  ci_lab <- paste0(round(100 * conf_level), "%")
  cat("\nR\u00B2 with", ci_lab, "intervals (boot = patient-level bootstrap; nb = corrected resampled t):\n")
  ci_print <- confidence_intervals
  ci_print[, -(1:2)] <- round(ci_print[, -(1:2)], 4)
  print(ci_print, row.names = FALSE)
  cat("\nContrasts between regimes with", ci_lab, "intervals:\n")
  ct_print <- contrast_intervals
  ct_print[, -(1:2)] <- round(ct_print[, -(1:2)], 4)
  print(ct_print, row.names = FALSE)
  if (!is.null(gsva_sens)) {
    cat("\nGSVA sensitivity (global GSVA vs training-reference GSVA; diff = global - reference):\n")
    print(gsva_sens, digits = 4, row.names = FALSE)
    cat("\nAgreement of scores, training-reference vs global GSVA (mean over folds):\n")
    print(round(colMeans(gsva_score_agreement[, -(1:2)]), 4))
  }
  if (!is.null(robustness_r2)) {
    cat("\nRobustness: pathway-specific adjustment and KRAS controls (fold-internal FS):\n")
    rr <- robustness_r2[, setdiff(names(robustness_r2), "description")]
    num <- vapply(rr, is.numeric, logical(1))
    rr[num] <- lapply(rr[num], round, 4)
    print(rr, row.names = FALSE)
    cat("\nPartial correlation of the baseline KRAS score (in-sample, patient bootstrap CI):\n")
    pc <- robustness_partial_cor
    pc[, c("partial_r", "boot_lower", "boot_upper")] <- round(pc[, c("partial_r", "boot_lower", "boot_upper")], 3)
    print(pc[, c("description", "partial_r", "boot_lower", "boot_upper")], row.names = FALSE)
  }
  if (!is.null(ssgsea_results)) {
    cat("\nSecondary sensitivity: ssGSEA (normalize = FALSE) vs GSVA:\n")
    print(round(ssgsea_agreement, 4), row.names = FALSE)
    sr <- ssgsea_results
    num <- vapply(sr, is.numeric, logical(1))
    sr[num] <- lapply(sr[num], round, 4)
    print(sr, row.names = FALSE)
    cat("\nFactorial decomposition with ssGSEA:\n")
    sc <- ssgsea_contrasts
    sc[, 3:5] <- round(sc[, 3:5], 4)
    print(sc, row.names = FALSE)
    cat("\nCoupling null with ssGSEA scores, cor(b, S):\n")
    print(ssgsea_coupling_null, digits = 4, row.names = FALSE)
    if (!is.null(ssgsea_partial_cor)) {
      cat("\nPartial correlation of the baseline KRAS score, ssGSEA scores (in-sample, patient bootstrap CI):\n")
      pc2 <- ssgsea_partial_cor
      pc2[, c("partial_r", "boot_lower", "boot_upper")] <- round(pc2[, c("partial_r", "boot_lower", "boot_upper")], 3)
      print(pc2[, c("description", "partial_r", "boot_lower", "boot_upper")], row.names = FALSE)
    }
  }
  cat("\nCoupling null, cor(b, S):\n")
  print(coupling_null, digits = 4, row.names = FALSE)
  if (!is.null(permutation_null)) {
    cat("\nPairing permutation null (cross-validated R\u00B2):\n")
    print(permutation_null, digits = 4, row.names = FALSE)
  }
  cat("\nDefinition - fully controlled:\n ", definitions[["fully_controlled"]], "\n")

  invisible(list(
    results_genome_wide = results_genome_wide,
    results_summary = results_summary,
    confidence_intervals = confidence_intervals,
    contrast_intervals = contrast_intervals,
    bootstrap_R2 = boot_R2,
    per_fold_R2 = per_fold_r2,
    per_repeat = per_repeat_all,
    paired_differences = paired_differences,
    paired_per_repeat = paired_per_repeat,
    gsva_sensitivity = gsva_sens,
    per_repeat_gsva_sensitivity = per_repeat_gsva,
    gsva_ref_validation = gsva_ref_validation,
    gsva_ref_settings = gsva_ref_settings,
    gsva_score_agreement = gsva_score_agreement,
    robustness_r2 = robustness_r2,
    robustness_partial_cor = robustness_partial_cor,
    ssgsea_results = ssgsea_results,
    ssgsea_contrasts = ssgsea_contrasts,
    ssgsea_agreement = ssgsea_agreement,
    ssgsea_coupling_null = ssgsea_coupling_null,
    ssgsea_partial_cor = ssgsea_partial_cor,
    coupling_null = coupling_null,
    coupling_null_distribution = r_null,
    permutation_null = permutation_null,
    permutation_null_distribution = null_distribution,
    leakage_audit = leakage_audit,
    definitions = definitions,
    R2_standard = stats::setNames(mean_r2("standard"), paste0("K_", top_k_values)),
    R2_coupling_adjusted = stats::setNames(mean_r2("coupling_adjusted"), paste0("K_", top_k_values)),
    R2_fold_internal = stats::setNames(mean_r2("fold_internal"), paste0("K_", top_k_values)),
    R2_fully_controlled = stats::setNames(mean_r2("fully_controlled"), paste0("K_", top_k_values)),
    fold_results = fold_results,
    global_gene_order = global_gene_order,
    folds = folds,
    inner_foldid = inner_foldid_list,
    groups = groups,
    baseline_score = b_main,
    cet_score = c_main,
    y = y_main,
    top_k_values = top_k_values,
    predictor_pathways = colnames(x_baseline_matrix),
    pathway_overlap = pathway_overlap,
    seed = seed,
    n_cores = n_cores
  ))
}
