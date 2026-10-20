#' Benchmark of alternative learners within the PTPM evaluation design
#'
#' Re-evaluates the suppression phenotype with learners other than the
#' elastic net (XGBoost, sparse PLS, linear SVR) using EXACTLY the same design
#' as run_ptpm(): same patient-grouped outer folds, same patient-grouped inner
#' folds, same 2 x 2 factorial of regimes, same fully nested coupling
#' adjustment, same metric and the same uncertainty summaries. The elastic-net
#' results are taken from the run_ptpm() output (not refitted), so every
#' comparison with the elastic net is paired: same folds, same repeats, same
#' bootstrap resamples of patients.
#'
#' REQUIREMENT
#' -----------
#' run_ptpm() must return the baseline predictor matrix. Add this element to
#' its returned list (section 26):
#'     x_baseline = x_baseline_matrix,
#' Alternatively pass the same matrix through the argument `x`. The function
#' checks that `x` is the matrix used by run_ptpm() (same models, same order,
#' same global variance ranking).
#'
#' EVALUATION DESIGN (identical to run_ptpm)
#' -----------------------------------------
#'   Regime (key)                          Feature selection (FS)        Coupling adjustment
#'   Standard           (standard)          pre-CV, on the full dataset   no
#'   Coupling-adjusted  (coupling_adjusted) pre-CV, on the full dataset   yes
#'   Fold-internal      (fold_internal)     inside outer training folds   no
#'   Fully controlled   (fully_controlled)  inside outer training folds   yes
#'
#' For every learner, hyperparameters are chosen ONLY by the patient-grouped
#' inner CV of run_ptpm() (`ptpm$inner_foldid`): the candidate with the lowest
#' mean squared error over all inner-validation models is refitted on the
#' whole outer training fold. With coupling adjustment, S is regressed on b in
#' each inner training fold and the coefficients are applied to the inner
#' validation fold; for the final fit the regression is estimated on the outer
#' training fold and applied unchanged to the held-out fold.
#'
#' Learner-specific choices:
#'   XGBoost : nrounds and max_depth tuned on a grid by the grouped inner CV
#'             (this replaces early stopping on random folds); one boosting
#'             run per max_depth, predictions at each nrounds from the same
#'             run. Tree models are invariant to feature scaling.
#'   sPLS    : keepX (same for all components) and ncomp tuned; mixOmics
#'             scales X with the training data of each fit only.
#'   SVR     : linear kernel, cost (and epsilon) tuned; X centred/scaled and
#'             y standardised with the training data of each fit only, so
#'             epsilon is on the standardised outcome scale.
#' Scaling is part of each fit in ALL regimes (as glmnet's internal
#' standardisation in run_ptpm), so regimes differ only in FS and coupling.
#'
#' METRICS AND UNCERTAINTY (identical to run_ptpm)
#' -----------------------------------------------
#' R2 = mean over repeats of 1 - SSE_r/SST_r on the pooled out-of-fold
#' predictions of repeat r. cor2 per repeat and the legacy pooled cor2 (the
#' metric of the earlier benchmark script) are also reported. Intervals:
#' (a) patient-level cluster bootstrap of out-of-fold predictions, with the
#' same resamples as run_ptpm (seed + 70000), primary; (b) corrected resampled
#' t interval (Nadeau & Bengio 2003) on per-fold R2, conservative check.
#'
#' @param ptpm Output of run_ptpm().
#' @param x Baseline predictor matrix used by run_ptpm (models x genes).
#' @param methods Learners to evaluate: any of "xgboost", "spls", "svr".
#' @param top_k_values Subset of ptpm$top_k_values to evaluate.
#' @param repeats Number of repeats (from the first) to evaluate; default all.
#' @param xgb_params Fixed XGBoost parameters.
#' @param xgb_nrounds_grid,xgb_max_depth_grid XGBoost tuning grid.
#' @param spls_keepX_grid,spls_ncomp_max sPLS tuning grid.
#' @param svr_cost_grid,svr_epsilon_grid SVR tuning grid.
#' @param n_boot Number of patient-level bootstrap resamples (0 = skip).
#' @param conf_level Confidence level of the intervals.
#' @param seed Random seed (default: the seed used by run_ptpm).
#' @param n_cores Number of parallel workers.
#' @param quick Logical; quick check for reviewers. Unless set explicitly,
#'   uses at most 2 repeats, smaller tuning grids (XGBoost nrounds 50/150/400,
#'   sPLS keepX 10/50/200, SVR cost 0.01/1/100) and n_boot = 200. Design,
#'   folds, regimes and metrics are unchanged; only precision is reduced.
#' @return A list (see the end of the function).
#' @export
run_ptpm_benchmark <- function(
    ptpm,
    x = ptpm$x_baseline,
    methods = c("xgboost", "spls", "svr"),
    top_k_values = ptpm$top_k_values,
    repeats = NULL,
    xgb_params = list(objective = "reg:squarederror", eta = 0.05,
                      subsample = 0.8, colsample_bytree = 0.8),
    xgb_nrounds_grid = c(25, 50, 100, 200, 400, 700, 1000),
    xgb_max_depth_grid = 3,
    spls_keepX_grid = c(5, 10, 25, 50, 100, 200, 300),
    spls_ncomp_max = 2,
    svr_cost_grid = c(0.01, 0.1, 1, 10, 100),
    svr_epsilon_grid = 0.1,
    n_boot = 2000,
    conf_level = 0.95,
    seed = ptpm$seed,
    n_cores = ptpm$n_cores,
    quick = FALSE
) {
  # 0. ARGUMENT CHECKS ---------------------------------------------------
  needed <- c("folds", "inner_foldid", "groups", "y", "baseline_score",
              "global_gene_order", "top_k_values", "fold_results", "seed")
  miss <- setdiff(needed, names(ptpm))
  if (length(miss)) stop("'ptpm' is not a run_ptpm() output; missing: ", paste(miss, collapse = ", "))
  if (is.null(x)) {
    stop("Baseline matrix not found. Add 'x_baseline = x_baseline_matrix' to the list ",
         "returned by run_ptpm() and re-run it, or pass the matrix through 'x'.")
  }
  methods <- unique(match.arg(methods, c("xgboost", "spls", "svr"), several.ok = TRUE))
  pkgs <- c(xgboost = "xgboost", spls = "mixOmics", svr = "e1071")[methods]
  for (pkg in c(unname(pkgs), "foreach", "doParallel")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.")
  }
  `%dopar%` <- foreach::`%dopar%`

  x <- as.matrix(x)
  y <- ptpm$y
  b <- ptpm$baseline_score
  if (!is.numeric(x) || anyNA(x)) stop("'x' must be a numeric matrix without missing values.")
  if (is.null(names(y)) || !identical(rownames(x), names(y)) || !identical(names(b), names(y))) {
    stop("Rows of 'x' must match names(ptpm$y) and names(ptpm$baseline_score), in the same order.")
  }
  y <- unname(y)
  b <- unname(b)
  global_gene_order <- ptpm$global_gene_order
  if (!identical(names(sort(apply(x, 2, stats::var), decreasing = TRUE)), global_gene_order)) {
    stop("'x' does not reproduce ptpm$global_gene_order: it is not the matrix used by run_ptpm().")
  }
  k_idx <- match(top_k_values, ptpm$top_k_values)
  if (length(top_k_values) == 0L || anyNA(k_idx)) {
    stop("'top_k_values' must be a subset of ptpm$top_k_values: ",
         paste(ptpm$top_k_values, collapse = ", "))
  }
  folds <- ptpm$folds
  inner_foldid_list <- ptpm$inner_foldid
  groups <- ptpm$groups
  n_obs <- nrow(x)
  nK <- length(top_k_values)
  fold_rep <- as.integer(sub(".*_Rep", "", names(folds)))
  max_rep <- max(fold_rep)
  if (!is.logical(quick) || length(quick) != 1L || is.na(quick)) stop("'quick' must be TRUE or FALSE.")
  if (quick) {
    # explicit arguments always win over the quick settings
    if (missing(repeats)) repeats <- min(2L, max_rep)
    if (missing(xgb_nrounds_grid)) xgb_nrounds_grid <- c(50, 150, 400)
    if (missing(spls_keepX_grid)) spls_keepX_grid <- c(10, 50, 200)
    if (missing(svr_cost_grid)) svr_cost_grid <- c(0.01, 1, 100)
    if (missing(n_boot)) n_boot <- 200
  }
  if (is.null(repeats)) repeats <- max_rep
  if (length(repeats) != 1L || repeats < 1 || repeats > max_rep || repeats != as.integer(repeats)) {
    stop("'repeats' must be an integer between 1 and ", max_rep, ".")
  }
  if (length(n_cores) != 1L || n_cores < 1) stop("'n_cores' must be >= 1.")
  reps <- seq_len(repeats)
  fold_idx <- which(fold_rep <= repeats)

  learners <- list(
    xgboost = function() .ptpm_learner_xgboost(xgb_params, xgb_nrounds_grid, xgb_max_depth_grid),
    spls    = function() .ptpm_learner_spls(spls_keepX_grid, spls_ncomp_max),
    svr     = function() .ptpm_learner_svr(svr_cost_grid, svr_epsilon_grid)
  )
  learners <- lapply(learners[methods], function(f) f())

  cat("\n============================================================\n")
  cat("PTPM benchmark of alternative learners", if (quick) "(QUICK CHECK: reduced precision)", "\n")
  cat("============================================================\n")
  cat("methods           :", paste(methods, collapse = ", "), "(reference: elastic net from run_ptpm)\n")
  cat("top_k_values      :", paste(top_k_values, collapse = ", "), "\n")
  cat("outer splits      :", length(fold_idx), "(", repeats, "repeats, same folds as run_ptpm )\n")
  cat("inner folds       : patient-grouped, same as run_ptpm\n")
  cat("seed / cores      :", seed, "/", n_cores, "\n")
  cat("============================================================\n\n")

  # 1. REGIMES -----------------------------------------------------------
  regimes <- data.frame(
    key   = c("standard", "coupling_adjusted", "fold_internal", "fully_controlled"),
    label = c("Standard R\u00B2", "Coupling-adjusted R\u00B2",
              "Fold-internal R\u00B2", "Fully controlled R\u00B2"),
    fs_inside        = c(FALSE, FALSE, TRUE, TRUE),
    coupling_correct = c(FALSE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  keys <- regimes$key

  # 2. PARALLEL BACKEND --------------------------------------------------
  if (n_cores > 1) {
    cat("Starting parallel backend with", n_cores, "cores...\n")
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit({ parallel::stopCluster(cl); foreach::registerDoSEQ() }, add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  eval_regime <- function(learner, fs_inside, coupling_correct) {
    # local copies: foreach auto-exports only objects of the calling frame
    x_full <- x; y_vec <- y; b_vec <- b
    folds_l <- folds; ifl <- inner_foldid_list
    gorder <- global_gene_order; kv <- top_k_values
    fitf <- .ptpm_fit_fold_learner
    base_seed <- seed
    foreach::foreach(i = fold_idx, .packages = learner$packages) %dopar% {
      tr <- folds_l[[i]]
      te <- setdiff(seq_len(nrow(x_full)), tr)
      res <- fitf(x_full[tr, , drop = FALSE], x_full[te, , drop = FALSE],
                  y_vec[tr], y_vec[te], b_vec[tr], b_vec[te],
                  fs_inside, gorder, coupling_correct, kv,
                  ifl[[i]], learner, base_seed + 1000L * i)
      list(fold = i, pred = res$pred, obs = res$obs,
           chosen = cbind(fold = i, repeat_id = as.integer(sub(".*_Rep", "", names(folds_l)[i])),
                          res$chosen))
    }
  }

  # 3. METRIC HELPERS (as in run_ptpm) ------------------------------------
  r2_sse <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
  r2_cor <- function(obs, pred) {
    if (stats::sd(pred) == 0) return(0)
    stats::cor(obs, pred)^2
  }
  per_repeat_metrics <- function(fr) {
    f_idx <- vapply(fr, `[[`, numeric(1), "fold")
    rr <- fold_rep[f_idx]
    out <- list()
    for (r in sort(unique(rr))) {
      sub <- fr[rr == r]
      obs <- unlist(lapply(sub, `[[`, "obs"))
      pred <- do.call(rbind, lapply(sub, `[[`, "pred"))
      for (j in seq_len(nK)) {
        out[[length(out) + 1L]] <- data.frame(
          repeat_id = r, K = top_k_values[j],
          R2 = r2_sse(obs, pred[, j]), cor2 = r2_cor(obs, pred[, j]))
      }
    }
    do.call(rbind, out)
  }
  pooled_cor2 <- function(fr) {
    obs <- unlist(lapply(fr, `[[`, "obs"))
    pred <- do.call(rbind, lapply(fr, `[[`, "pred"))
    vapply(seq_len(nK), function(j) r2_cor(obs, pred[, j]), numeric(1))
  }
  summarise_vec <- function(v) {
    c(mean = mean(v), sd = stats::sd(v),
      q025 = unname(stats::quantile(v, 0.025)), q975 = unname(stats::quantile(v, 0.975)),
      prop_positive = mean(v > 0))
  }

  # 4. ELASTIC NET REFERENCE (from run_ptpm, same folds) ---------------------
  fold_results <- list()
  fold_results$elastic_net <- lapply(stats::setNames(keys, keys), function(key) {
    fr <- ptpm$fold_results[[key]]
    if (is.null(fr)) stop("ptpm$fold_results$", key, " is missing.")
    fr <- Filter(function(z) z$fold %in% fold_idx, fr)
    lapply(fr, function(z) list(fold = z$fold, pred = z$pred[, k_idx, drop = FALSE], obs = z$obs))
  })

  # 5. ALTERNATIVE LEARNERS: FOUR REGIMES EACH ------------------------------
  tuning <- list()
  for (m in methods) {
    fold_results[[m]] <- list()
    for (q in seq_len(nrow(regimes))) {
      cat("-- ", m, " | regime: ", regimes$key[q],
          " | FS ", if (regimes$fs_inside[q]) "inside folds" else "pre-CV",
          " | coupling correction ", if (regimes$coupling_correct[q]) "yes" else "no", "\n", sep = "")
      fr <- eval_regime(learners[[m]], regimes$fs_inside[q], regimes$coupling_correct[q])
      fold_results[[m]][[regimes$key[q]]] <- lapply(fr, function(z) z[c("fold", "pred", "obs")])
      ch <- do.call(rbind, lapply(fr, `[[`, "chosen"))
      tuning[[length(tuning) + 1L]] <- data.frame(method = m, regime = regimes$key[q], ch,
                                                  row.names = NULL)
      pr <- per_repeat_metrics(fr)
      agg <- stats::aggregate(R2 ~ K, pr, mean)
      print(stats::setNames(round(agg$R2, 4), paste0("K_", agg$K)))
    }
  }
  all_methods <- c("elastic_net", methods)
  tuning <- do.call(rbind, tuning)

  per_repeat <- do.call(rbind, lapply(all_methods, function(m) {
    do.call(rbind, lapply(keys, function(key) {
      pr <- per_repeat_metrics(fold_results[[m]][[key]])
      cbind(method = m, regime = key, pr, stringsAsFactors = FALSE)
    }))
  }))
  rownames(per_repeat) <- NULL
  pr_get <- function(m, key, K) {
    p <- per_repeat[per_repeat$method == m & per_repeat$regime == key & per_repeat$K == K, ]
    p[order(p$repeat_id), ]
  }
  mean_r2 <- function(m, key) vapply(top_k_values, function(K) mean(pr_get(m, key, K)$R2), numeric(1))

  # 6. SUMMARY TABLES ----------------------------------------------------
  results_table <- do.call(rbind, lapply(all_methods, function(m) {
    d <- data.frame(Method = m, K = top_k_values,
                    a = mean_r2(m, "standard"), b = mean_r2(m, "coupling_adjusted"),
                    c = mean_r2(m, "fold_internal"), d = mean_r2(m, "fully_controlled"),
                    stringsAsFactors = FALSE)
    colnames(d) <- c("Method", "K", regimes$label)
    d
  }))
  rownames(results_table) <- NULL

  results_summary <- do.call(rbind, lapply(all_methods, function(m) {
    do.call(rbind, lapply(keys, function(key) {
      legacy <- pooled_cor2(fold_results[[m]][[key]])
      do.call(rbind, lapply(seq_len(nK), function(j) {
        p <- pr_get(m, key, top_k_values[j])
        s_r2 <- summarise_vec(p$R2)
        s_c2 <- summarise_vec(p$cor2)
        data.frame(method = m, regime = key, K = top_k_values[j],
                   R2_mean = s_r2["mean"], R2_sd = s_r2["sd"],
                   R2_q025 = s_r2["q025"], R2_q975 = s_r2["q975"],
                   cor2_mean = s_c2["mean"], cor2_sd = s_c2["sd"],
                   cor2_pooled_legacy = legacy[j], row.names = NULL)
      }))
    }))
  }))

  # 7. CONFIDENCE INTERVALS (as in run_ptpm) ------------------------------
  alpha_ci <- (1 - conf_level) / 2

  # (a) patient-level cluster bootstrap, same resamples as run_ptpm
  oof_arrays <- function(fr) {
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
  boot_R2 <- NULL
  if (n_boot > 0) {
    cat("\n-- Patient-level bootstrap:", n_boot, "resamples...\n")
    patients <- unique(groups)
    pat_idx <- match(groups, patients)
    set.seed(seed + 70000)
    W <- vapply(seq_len(n_boot), function(bb) {
      cnt <- tabulate(sample(length(patients), replace = TRUE), nbins = length(patients))
      as.numeric(cnt[pat_idx])
    }, numeric(n_obs))
    sw <- colSums(W)
    boot_mean_r2 <- function(fr) {
      arr <- oof_arrays(fr)
      O <- arr$O
      mo <- crossprod(W, O) / sw
      sst <- crossprod(W, O^2) - sw * mo^2
      out <- vapply(seq_len(nK), function(j) {
        Pj <- matrix(arr$P[, , j], nrow = n_obs)
        sse <- crossprod(W, (O - Pj)^2)
        rowMeans(1 - sse / sst)
      }, numeric(n_boot))
      matrix(out, nrow = n_boot)
    }
    boot_R2 <- lapply(stats::setNames(all_methods, all_methods), function(m) {
      lapply(stats::setNames(keys, keys), function(key) boot_mean_r2(fold_results[[m]][[key]]))
    })
  }
  boot_q <- function(v) {
    c(boot_lower = unname(stats::quantile(v, alpha_ci)),
      boot_upper = unname(stats::quantile(v, 1 - alpha_ci)),
      boot_se = stats::sd(v))
  }
  boot_na <- c(boot_lower = NA_real_, boot_upper = NA_real_, boot_se = NA_real_)

  # (b) corrected resampled t (Nadeau & Bengio 2003) on per-fold R2
  per_fold_r2 <- lapply(stats::setNames(all_methods, all_methods), function(m) {
    lapply(stats::setNames(keys, keys), function(key) {
      fr <- fold_results[[m]][[key]]
      fr <- fr[order(vapply(fr, `[[`, numeric(1), "fold"))]
      matrix(t(vapply(fr, function(z) {
        vapply(seq_len(nK), function(j) r2_sse(z$obs, z$pred[, j]), numeric(1))
      }, numeric(nK))), ncol = nK)
    })
  })
  ratio_te_tr <- mean(vapply(folds[fold_idx], function(tr) (n_obs - length(tr)) / length(tr),
                             numeric(1)))
  nb_interval <- function(d) {
    J <- length(d); mm <- mean(d)
    se <- sqrt((1 / J + ratio_te_tr) * stats::var(d))
    tq <- stats::qt(1 - alpha_ci, J - 1)
    c(nb_mean = mm, nb_lower = mm - tq * se, nb_upper = mm + tq * se,
      nb_p = if (se > 0) 2 * stats::pt(-abs(mm / se), J - 1) else NA_real_)
  }

  confidence_intervals <- do.call(rbind, lapply(all_methods, function(m) {
    do.call(rbind, lapply(keys, function(key) {
      do.call(rbind, lapply(seq_len(nK), function(j) {
        bq <- if (!is.null(boot_R2)) boot_q(boot_R2[[m]][[key]][, j]) else boot_na
        nb <- nb_interval(per_fold_r2[[m]][[key]][, j])
        data.frame(method = m, regime = key, K = top_k_values[j], R2 = mean_r2(m, key)[j],
                   t(bq), t(nb[c("nb_mean", "nb_lower", "nb_upper")]), row.names = NULL)
      }))
    }))
  }))

  # 8. FACTORIAL CONTRASTS WITHIN EACH METHOD (as in run_ptpm) --------------
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
  contrast_intervals <- do.call(rbind, lapply(all_methods, function(m) {
    do.call(rbind, lapply(names(contrasts), function(nm) {
      do.call(rbind, lapply(seq_len(nK), function(j) {
        est <- contrasts[[nm]](lapply(stats::setNames(keys, keys), function(key) mean_r2(m, key)[j]))
        pr_w <- lapply(stats::setNames(keys, keys), function(key) pr_get(m, key, top_k_values[j])$R2)
        s <- summarise_vec(contrasts[[nm]](pr_w))
        if (!is.null(boot_R2)) {
          bd <- contrasts[[nm]](lapply(stats::setNames(keys, keys),
                                       function(key) boot_R2[[m]][[key]][, j]))
          bq <- c(boot_q(bd), boot_prop_positive = mean(bd > 0))
        } else {
          bq <- c(boot_na, boot_prop_positive = NA_real_)
        }
        nb <- nb_interval(contrasts[[nm]](lapply(stats::setNames(keys, keys),
                                                 function(key) per_fold_r2[[m]][[key]][, j])))
        data.frame(method = m, contrast = nm, K = top_k_values[j], estimate = est,
                   sd_repeats = s["sd"], prop_repeats_positive = s["prop_positive"],
                   t(bq), t(nb), row.names = NULL)
      }))
    }))
  }))

  # 9. PAIRED COMPARISON WITH THE ELASTIC NET (same folds, same resamples) ---
  method_comparison <- do.call(rbind, lapply(methods, function(m) {
    do.call(rbind, lapply(keys, function(key) {
      do.call(rbind, lapply(seq_len(nK), function(j) {
        K <- top_k_values[j]
        pm <- pr_get(m, key, K)
        pe <- pr_get("elastic_net", key, K)
        pe <- pe[match(pm$repeat_id, pe$repeat_id), ]
        s <- summarise_vec(pm$R2 - pe$R2)
        if (!is.null(boot_R2)) {
          bd <- boot_R2[[m]][[key]][, j] - boot_R2$elastic_net[[key]][, j]
          bq <- c(boot_q(bd), boot_prop_positive = mean(bd > 0))
        } else {
          bq <- c(boot_na, boot_prop_positive = NA_real_)
        }
        nb <- nb_interval(per_fold_r2[[m]][[key]][, j] - per_fold_r2$elastic_net[[key]][, j])
        data.frame(method = m, regime = key, K = K,
                   R2_method = mean(pm$R2), R2_elastic_net = mean(pe$R2),
                   diff_mean = s["mean"], diff_sd_repeats = s["sd"],
                   prop_repeats_method_better = s["prop_positive"],
                   t(bq), t(nb[c("nb_lower", "nb_upper", "nb_p")]), row.names = NULL)
      }))
    }))
  }))

  # 10. LEAKAGE AUDIT ------------------------------------------------------
  leakage_audit <- do.call(rbind, lapply(all_methods, function(m) {
    tun <- if (m == "elastic_net") "alpha, lambda" else learners[[m]]$tuning
    sc <- if (m == "elastic_net") "glmnet internal standardisation (training data of each fit)"
          else learners[[m]]$scaling
    data.frame(
      method = m,
      step = c("Outer train/test split", "Inner (tuning) split", "Hyperparameters",
               "Scaling", "Gene ranking (feature selection)", "Coupling regression (S ~ b)",
               "GSVA pathway scores (S, b)"),
      Standard = c("patient-grouped", "patient-grouped (run_ptpm folds)", "inner CV only",
                   "within each fit", "full dataset (pre-CV)", "not applied", "global (whole cohort)"),
      Coupling_adjusted = c("patient-grouped", "patient-grouped (run_ptpm folds)", "inner CV only",
                            "within each fit", "full dataset (pre-CV)",
                            "outer and inner training folds", "global (whole cohort)"),
      Fold_internal = c("patient-grouped", "patient-grouped (run_ptpm folds)", "inner CV only",
                        "within each fit", "outer training fold", "not applied", "global (whole cohort)"),
      Fully_controlled = c("patient-grouped", "patient-grouped (run_ptpm folds)", "inner CV only",
                           "within each fit", "outer training fold",
                           "outer and inner training folds", "global (whole cohort)"),
      note = c("", "", tun, sc, "variance ranking is unsupervised (does not use S)",
               "re-estimated within each inner training fold during tuning",
               "unsupervised, transductive; see run_ptpm gsva_sensitivity"),
      stringsAsFactors = FALSE)
  }))

  # 11. PRINT & RETURN -----------------------------------------------------
  ci_lab <- paste0(round(100 * conf_level), "%")
  rnd <- function(d) { num <- vapply(d, is.numeric, logical(1)); d[num] <- lapply(d[num], round, 4); d }
  cat("\n============================================================\n")
  cat("BENCHMARK RESULTS (R\u00B2 = mean of", repeats, "per-repeat R\u00B2, 1 - SSE/SST)\n")
  cat("============================================================\n")
  print(rnd(results_table), row.names = FALSE)
  cat("\nR\u00B2 with", ci_lab, "intervals (boot = patient-level bootstrap; nb = corrected resampled t):\n")
  print(rnd(confidence_intervals), row.names = FALSE)
  cat("\nFactorial decomposition (Shapley) per method:\n")
  sh <- contrast_intervals[contrast_intervals$contrast %in%
                             c("Delta_FS_shapley", "Delta_coupling_shapley", "Interaction", "Total_gap"),
                           c("method", "contrast", "K", "estimate", "boot_lower", "boot_upper")]
  print(rnd(sh), row.names = FALSE)
  cat("\nMethod minus elastic net (paired: same folds, same bootstrap resamples):\n")
  print(rnd(method_comparison[, c("method", "regime", "K", "R2_method", "R2_elastic_net",
                                  "diff_mean", "boot_lower", "boot_upper",
                                  "prop_repeats_method_better")]), row.names = FALSE)

  invisible(list(
    results_table = results_table,
    results_summary = results_summary,
    confidence_intervals = confidence_intervals,
    contrast_intervals = contrast_intervals,
    method_comparison = method_comparison,
    per_repeat = per_repeat,
    per_fold_R2 = per_fold_r2,
    bootstrap_R2 = boot_R2,
    tuning = tuning,
    leakage_audit = leakage_audit,
    fold_results = fold_results,
    regimes = regimes,
    top_k_values = top_k_values,
    repeats = repeats,
    seed = seed
  ))
}

# ---------------------------------------------------------------------------
# Core fold fit for a generic learner (self-contained: it is shipped to the
# parallel workers and must not depend on other helpers).
# Mirrors fit_fold() of run_ptpm: gene ranking (pre-CV or inside the outer
# training fold), patient-grouped inner CV for tuning, coupling regression
# re-estimated within each inner training fold, final coupling regression
# estimated on the outer training fold and applied unchanged to the test fold.
# ---------------------------------------------------------------------------
.ptpm_fit_fold_learner <- function(x_train, x_test, y_train, y_test, b_train, b_test,
                                   fs_inside, global_gene_order, coupling_correct,
                                   top_k_values, inner_foldid, learner, seed) {
  ols <- function(B1, yy) {
    cf <- stats::lm.fit(B1, yy)$coefficients
    cf[is.na(cf)] <- 0
    cf
  }
  gene_order <- if (fs_inside) {
    names(sort(apply(x_train, 2, stats::var), decreasing = TRUE))
  } else {
    global_gene_order
  }
  Btr <- cbind(1, as.matrix(b_train))
  Bte <- cbind(1, as.matrix(b_test))
  if (coupling_correct) {
    cf <- ols(Btr, y_train)
    t_train <- as.numeric(y_train - Btr %*% cf)
    t_test <- as.numeric(y_test - Bte %*% cf)
  } else {
    t_train <- y_train
    t_test <- y_test
  }
  # Inner targets do not depend on K: computed once per outer fold
  inner <- lapply(sort(unique(inner_foldid)), function(f) {
    it <- inner_foldid != f
    iv <- !it
    if (coupling_correct) {
      cf_f <- ols(Btr[it, , drop = FALSE], y_train[it])
      list(it = it, iv = iv,
           t_it = as.numeric(y_train[it] - Btr[it, , drop = FALSE] %*% cf_f),
           t_iv = as.numeric(y_train[iv] - Btr[iv, , drop = FALSE] %*% cf_f))
    } else {
      list(it = it, iv = iv, t_it = y_train[it], t_iv = y_train[iv])
    }
  })
  pred <- matrix(NA_real_, nrow(x_test), length(top_k_values))
  chosen <- vector("list", length(top_k_values))
  for (j in seq_along(top_k_values)) {
    top <- gene_order[seq_len(min(top_k_values[j], length(gene_order)))]
    xtr <- x_train[, top, drop = FALSE]
    xte <- x_test[, top, drop = FALSE]
    grid <- learner$grid(length(top))
    err <- matrix(NA_real_, length(y_train), nrow(grid))
    for (q in seq_along(inner)) {
      p <- learner$fit_predict(xtr[inner[[q]]$it, , drop = FALSE], inner[[q]]$t_it,
                               xtr[inner[[q]]$iv, , drop = FALSE], grid, seed + q)
      err[inner[[q]]$iv, ] <- (inner[[q]]$t_iv - p)^2
    }
    cvm <- colMeans(err)
    cvm[!is.finite(cvm)] <- Inf
    best <- which.min(cvm)
    pred[, j] <- learner$fit_predict(xtr, t_train, xte, grid[best, , drop = FALSE], seed)[, 1]
    g <- grid[best, , drop = FALSE]
    chosen[[j]] <- data.frame(K = top_k_values[j],
                              hyperparameters = paste(names(g), unlist(g), sep = " = ",
                                                      collapse = ", "),
                              cv_mse = cvm[best], stringsAsFactors = FALSE)
  }
  list(pred = pred, obs = unname(t_test), chosen = do.call(rbind, chosen))
}

# ---------------------------------------------------------------------------
# Learners. Each returns a list with:
#   grid(K)                                   candidate hyperparameters
#   fit_predict(x_tr, y_tr, x_new, grid, seed) n_new x nrow(grid) predictions
# Arguments are forced so that the closures do not keep a reference to the
# calling frame (which would be serialised to every worker).
# ---------------------------------------------------------------------------
.ptpm_learner_xgboost <- function(params, nrounds_grid, max_depth_grid) {
  force(params); force(nrounds_grid); force(max_depth_grid)
  # iterationrange is [begin, end) in xgboost 1.x and [begin, end] in >= 3.x:
  # detect the convention once on a tiny deterministic model.
  inclusive <- local({
    xx <- matrix(sin(seq_len(200)), 50)
    d <- xgboost::xgb.DMatrix(xx, label = cos(seq_len(50)))
    m <- xgboost::xgb.train(params = list(objective = "reg:squarederror", max_depth = 2,
                                          eta = 0.3, nthread = 1),
                            data = d, nrounds = 3, verbose = 0)
    isTRUE(all.equal(stats::predict(m, d), stats::predict(m, d, iterationrange = c(1, 3))))
  })
  list(
    name = "xgboost",
    packages = "xgboost",
    tuning = paste0("nrounds {", paste(nrounds_grid, collapse = ", "), "} x max_depth {",
                    paste(max_depth_grid, collapse = ", "), "}"),
    scaling = "none needed (trees are invariant to monotone feature transformations)",
    grid = function(K) expand.grid(max_depth = max_depth_grid, nrounds = nrounds_grid,
                                   KEEP.OUT.ATTRS = FALSE),
    fit_predict = function(x_tr, y_tr, x_new, grid, seed) {
      out <- matrix(NA_real_, nrow(x_new), nrow(grid))
      dtr <- xgboost::xgb.DMatrix(x_tr, label = y_tr)
      dnew <- xgboost::xgb.DMatrix(x_new)
      for (d in unique(grid$max_depth)) {
        rows <- which(grid$max_depth == d)
        p <- utils::modifyList(params, list(max_depth = d, nthread = 1, seed = seed))
        set.seed(seed)
        m <- xgboost::xgb.train(params = p, data = dtr,
                                nrounds = max(grid$nrounds[rows]), verbose = 0)
        for (r in rows) {
          n <- grid$nrounds[r]
          out[, r] <- stats::predict(m, dnew, iterationrange = c(1, if (inclusive) n else n + 1))
        }
      }
      out
    }
  )
}

.ptpm_learner_spls <- function(keepX_grid, ncomp_max) {
  force(keepX_grid); force(ncomp_max)
  list(
    name = "spls",
    packages = "mixOmics",
    tuning = paste0("keepX {", paste(keepX_grid, collapse = ", "),
                    "} (<= K, same for all components) x ncomp 1..", ncomp_max),
    scaling = "mixOmics scale = TRUE (training data of each fit)",
    grid = function(K) {
      kx <- keepX_grid[keepX_grid <= K]
      if (length(kx) == 0L) kx <- K
      expand.grid(keepX = kx, ncomp = seq_len(ncomp_max), KEEP.OUT.ATTRS = FALSE)
    },
    fit_predict = function(x_tr, y_tr, x_new, grid, seed) {
      out <- matrix(NA_real_, nrow(x_new), nrow(grid))
      y_mat <- matrix(y_tr, ncol = 1L, dimnames = list(rownames(x_tr), "S"))
      for (kx in unique(grid$keepX)) {
        rows <- which(grid$keepX == kx)
        nc <- max(grid$ncomp[rows])
        fit <- mixOmics::spls(X = x_tr, Y = y_mat, ncomp = nc, keepX = rep(kx, nc),
                              mode = "regression", scale = TRUE)
        pr <- stats::predict(fit, newdata = x_new)$predict
        for (r in rows) out[, r] <- pr[, 1, grid$ncomp[r]]
      }
      out
    }
  )
}

.ptpm_learner_svr <- function(cost_grid, epsilon_grid) {
  force(cost_grid); force(epsilon_grid)
  list(
    name = "svr",
    packages = "e1071",
    tuning = paste0("cost {", paste(cost_grid, collapse = ", "), "} x epsilon {",
                    paste(epsilon_grid, collapse = ", "), "} (standardised y)"),
    scaling = "X centred/scaled and y standardised on the training data of each fit",
    grid = function(K) expand.grid(cost = cost_grid, epsilon = epsilon_grid,
                                   KEEP.OUT.ATTRS = FALSE),
    fit_predict = function(x_tr, y_tr, x_new, grid, seed) {
      mu <- colMeans(x_tr)
      s <- apply(x_tr, 2, stats::sd)
      s[!is.finite(s) | s == 0] <- 1
      xs_tr <- scale(x_tr, center = mu, scale = s)
      xs_new <- scale(x_new, center = mu, scale = s)
      my <- mean(y_tr)
      sy <- stats::sd(y_tr)
      if (!is.finite(sy) || sy == 0) sy <- 1
      ys <- (y_tr - my) / sy
      out <- matrix(NA_real_, nrow(x_new), nrow(grid))
      for (r in seq_len(nrow(grid))) {
        m <- e1071::svm(x = xs_tr, y = ys, type = "eps-regression", kernel = "linear",
                        cost = grid$cost[r], epsilon = grid$epsilon[r], scale = FALSE)
        out[, r] <- as.numeric(stats::predict(m, xs_new)) * sy + my
      }
      out
    }
  )
}
