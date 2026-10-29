# =============================================================================
# Shared internal helpers (not exported)
# Used by the training-reference GSVA sensitivity analysis of all run_ptpm*()
# functions.
# =============================================================================
# -----------------------------------------------------------------------------
# Reference-based GSVA (used by the GSVA sensitivity analysis of run_ptpm)
# -----------------------------------------------------------------------------
# GSVA has two steps:
#   (1) cross-sample: per-gene CDF of expression (Gaussian kernel, bandwidth
#       sd/4, or direct ECDF), the only step that uses other samples;
#   (2) per-sample: ranking of genes within each sample and random-walk
#       statistic on each gene set (tau = 1, maxDiff = TRUE).
# Step (1) is estimated on REFERENCE samples only (e.g. the training fold)
# and evaluated on ALL samples. With reference = all samples these functions
# reproduce GSVA::gsva(); run_ptpm() checks this on the user's data before
# using them.
# =============================================================================

.ptpm_cache <- new.env(parent = emptyenv())

# Gaussian-kernel CDF compiled with Rcpp (same arithmetic as GSVA's C code:
# bandwidth sd/4, pnorm lookup table with 10,001 entries on [0, 10]).
# Returns NULL when Rcpp or a compiler is not available.
.ptpm_kcdf_cpp <- function() {
  if (exists("f", envir = .ptpm_cache, inherits = FALSE)) return(.ptpm_cache$f)
  f <- NULL
  if (requireNamespace("Rcpp", quietly = TRUE)) {
    code <- '
NumericMatrix ptpm_kcdf_gauss(NumericMatrix X, NumericMatrix XR) {
  int G = X.nrow(), N = X.ncol(), M = XR.ncol();
  std::vector<double> tab(10001);
  for (int i = 0; i <= 10000; i++) tab[i] = R::pnorm(10.0 * i / 10000.0, 0.0, 1.0, 1, 0);
  NumericMatrix out(G, N);
  for (int g = 0; g < G; g++) {
    double mu = 0.0;
    for (int m = 0; m < M; m++) mu += XR(g, m);
    mu /= M;
    double ss = 0.0;
    for (int m = 0; m < M; m++) ss += (XR(g, m) - mu) * (XR(g, m) - mu);
    double bw = std::sqrt(ss / (M - 1)) / 4.0;
    for (int j = 0; j < N; j++) {
      double left = 0.0;
      for (int m = 0; m < M; m++) {
        double v = (X(g, j) - XR(g, m)) / bw;
        if (v < -10.0) continue;
        if (v > 10.0) { left += 1.0; continue; }
        double c = tab[(int)(std::fabs(v) / 10.0 * 10000.0)];
        left += (v < 0) ? 1.0 - c : c;
      }
      out(g, j) = left / M;
    }
  }
  return out;
}'
    f <- tryCatch(Rcpp::cppFunction(code, env = new.env()),
                  error = function(e) NULL)
  }
  assign("f", f, envir = .ptpm_cache)
  f
}

# Pure-R fallback of the Gaussian-kernel CDF (same result, slower).
.ptpm_kcdf_gauss_block <- function(X, XR) {
  tab <- stats::pnorm(10 * (0:10000) / 10000)
  bw <- apply(XR, 1, stats::sd) / 4
  acc <- matrix(0, nrow(X), ncol(X))
  for (m in seq_len(ncol(XR))) {
    v <- (X - XR[, m]) / bw
    cdf <- tab[floor(pmin(abs(v), 10) / 10 * 10000) + 1]
    # |v| > 10 gives pnorm(10) = 1 in double precision: 1 (v > 10), 0 (v < -10)
    acc <- acc + ifelse(v < 0, 1 - cdf, cdf)
  }
  acc / ncol(XR)
}

# Direct ECDF of each gene, estimated on reference samples.
.ptpm_ecdf_block <- function(X, XR) {
  out <- matrix(0, nrow(X), ncol(X))
  for (g in seq_len(nrow(X))) {
    out[g, ] <- findInterval(X[g, ], sort(XR[g, ])) / ncol(XR)
  }
  out
}

# Per-gene CDF matrix (genes x samples) with reference columns ref_idx.
.ptpm_gene_cdf <- function(expr, ref_idx, method = c("gaussian", "none"), cl = NULL) {
  method <- match.arg(method)
  XR <- expr[, ref_idx, drop = FALSE]
  if (method == "gaussian") {
    f_cpp <- .ptpm_kcdf_cpp()
    if (!is.null(f_cpp)) {
      out <- f_cpp(expr, XR)
      dimnames(out) <- dimnames(expr)
      return(out)
    }
    fun <- .ptpm_kcdf_gauss_block
  } else {
    fun <- .ptpm_ecdf_block
  }
  n_blocks <- if (is.null(cl)) 1L else length(cl)
  blocks <- if (n_blocks > 1L) {
    split(seq_len(nrow(expr)), cut(seq_len(nrow(expr)), n_blocks, labels = FALSE))
  } else list(seq_len(nrow(expr)))
  res <- if (n_blocks == 1L) {
    list(fun(expr, XR))
  } else {
    parallel::clusterMap(cl, function(X, XR, fun) fun(X, XR),
                         lapply(blocks, function(b) expr[b, , drop = FALSE]),
                         lapply(blocks, function(b) XR[b, , drop = FALSE]),
                         MoreArgs = list(fun = fun))
  }
  out <- do.call(rbind, res)
  dimnames(out) <- dimnames(expr)
  out
}

# Per-sample ranking + random walk (tau = 1, maxDiff = TRUE, absRanking = FALSE).
.ptpm_gsva_scores_from_cdf <- function(cdf, gset_idx, ties = c("last", "first"),
                                       tau = 1) {
  ties <- match.arg(ties)
  p <- nrow(cdf)
  sc <- matrix(NA_real_, length(gset_idx), ncol(cdf),
               dimnames = list(names(gset_idx), colnames(cdf)))
  for (j in seq_len(ncol(cdf))) {
    r <- rank(cdf[, j], ties.method = ties)
    dos <- p - r + 1
    srs <- abs(p / 2 - r)
    for (s in seq_along(gset_idx)) {
      g <- gset_idx[[s]]
      pos <- dos[g]
      step_in <- numeric(p)
      step_in[pos] <- srs[g]^tau
      step_in <- cumsum(step_in)
      step_out <- rep(1, p)
      step_out[pos] <- 0
      step_out <- cumsum(step_out)
      walk <- step_in / step_in[p] - step_out / step_out[p]
      sc[s, j] <- max(0, max(walk)) + min(0, min(walk))
    }
  }
  sc
}

# Reference-based GSVA. Genes constant in the reference samples are removed
# (GSVA removes constant genes of its input). 'cdf' can be supplied to reuse
# a CDF matrix computed with the same reference.
.ptpm_gsva_ref <- function(expr, gene_sets, ref_idx, method, ties, cl = NULL,
                           cdf = NULL) {
  keep <- apply(expr[, ref_idx, drop = FALSE], 1, stats::sd) > 0
  expr <- expr[keep, , drop = FALSE]
  gset_idx <- lapply(gene_sets, function(g) unique(stats::na.omit(match(g, rownames(expr)))))
  gset_idx <- gset_idx[lengths(gset_idx) > 0]
  if (is.null(cdf)) cdf <- .ptpm_gene_cdf(expr, ref_idx, method, cl)
  list(scores = .ptpm_gsva_scores_from_cdf(cdf, gset_idx, ties), cdf = cdf)
}
