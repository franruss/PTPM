## =====================================================================
## PTPM 0.1.6: R² on the scale of S and Oldham's test
## =====================================================================
## Reproduces Supplementary Table S2.2, the "Incremental R² on the scale
## of S" column of Table 1 and Oldham's test reported in the manuscript.
##
## Input (working directory):
##   ptpm_genome_wide.rds, ptpm_egfr_axis.rds, ptpm_pathways.rds,
##   ptpm_hybrid.rds, ptpm_hybrid_interaction.rds
##     saved by PTPM_0.1.6_complete_results.R (user guide, section 4)
##   ptpm_positive_control.rds (optional)
##     saved by PTPM_0.1.6_positive_control.R (user guide, section 6.1)
##
## Output:
##   console output (saved as PTPM_0.1.6_scale_of_S_and_Oldham.txt)
##   PTPM_0.1.6_scale_of_S.csv
##
## No model is refitted: the script only reuses the out-of-fold
## predictions stored in the run_ptpm*() results. Runtime: a few minutes.
## =====================================================================

options(width = 200)

# ---------------------------------------------------------------------
# R² = 1 - SSE/SST
# ---------------------------------------------------------------------
r2 <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)

# ---------------------------------------------------------------------
# Out-of-fold predictions of the fully controlled regime on the scale of S.
# In coupling-adjusted regimes z$obs is the test residual
#   R_j = S_j - beta0(T) - beta1(T) * b_j,
# with the coupling regression estimated on the outer training fold T.
# Hence S_j - R_j is the prediction from baseline pathway activity alone,
# and S_j - R_j + Rhat_j the prediction from baseline activity + predictors.
# ---------------------------------------------------------------------
s_scale_rows <- function(x, regime = "fully_controlled", K = 3000) {
  j <- match(K, x$top_k_values)
  if (is.na(j)) stop("K = ", K, " not found in top_k_values.")
  n <- length(x$y)
  do.call(rbind, lapply(x$fold_results[[regime]], function(z) {
    te    <- setdiff(seq_len(n), x$folds[[z$fold]])
    S     <- unname(x$y[te])
    fit_b <- S - z$obs
    data.frame(rep      = as.integer(sub(".*_Rep", "", names(x$folds)[z$fold])),
               patient  = x$groups[te],
               S        = S,
               b_only   = fit_b,
               b_plus_X = fit_b + z$pred[, j])
  }))
}

# ---------------------------------------------------------------------
# Mean over repeats of R² on S, with patient-level bootstrap 95% CI
# (patients resampled jointly across repeats, as in run_ptpm()).
# ---------------------------------------------------------------------
s_scale_summary <- function(x, regime = "fully_controlled", K = 3000,
                            B = 2000, seed = 123) {
  d <- s_scale_rows(x, regime, K)
  est <- function(dd) {
    v <- sapply(split(dd, dd$rep), function(z)
      c(r2(z$S, z$b_only), r2(z$S, z$b_plus_X)))
    m <- rowMeans(v)
    c(m[1], m[2], m[2] - m[1])
  }
  rows <- split(seq_len(nrow(d)), d$patient)
  set.seed(seed)
  bs <- replicate(B, est(d[unlist(rows[sample(names(rows), replace = TRUE)],
                                  use.names = FALSE), ]))
  out <- rbind(estimate = est(d),
               lower    = apply(bs, 1, quantile, 0.025),
               upper    = apply(bs, 1, quantile, 0.975))
  colnames(out) <- c("R2_b_only", "R2_b_plus_X", "incremental")
  out
}

# ---------------------------------------------------------------------
# 1. R² on the scale of S (Supplementary Table S2.2)
# ---------------------------------------------------------------------
models <- data.frame(
  file  = c("ptpm_genome_wide.rds", "ptpm_genome_wide.rds",
            "ptpm_positive_control.rds", "ptpm_egfr_axis.rds",
            "ptpm_pathways.rds", "ptpm_hybrid.rds",
            "ptpm_hybrid_interaction.rds"),
  model = c("fit_reg", "fit_reg", "fit_reg, supervised selection",
            "fit_EGFR_genes", "fit_pathways", "fit_hybrid",
            "fit_hybrid_interactions"),
  K     = c(3000, 300, 3000, 569, 47, 300, 300),
  stringsAsFactors = FALSE
)
# fit_weighted_score is not included: its outcome is rescaled within
# each training fold, so out-of-fold predictions are not on a common scale.

cache <- list()
tab <- list()
for (i in seq_len(nrow(models))) {
  f <- models$file[i]
  if (!file.exists(f)) {
    message("Skipping ", models$model[i], " (K = ", models$K[i], "): ", f, " not found.")
    next
  }
  if (is.null(cache[[f]])) cache[[f]] <- readRDS(f)
  s <- s_scale_summary(cache[[f]], K = models$K[i])
  tab[[length(tab) + 1]] <- data.frame(
    model             = models$model[i],
    K                 = models$K[i],
    R2_b_only         = s["estimate", "R2_b_only"],
    R2_b_only_lower   = s["lower",    "R2_b_only"],
    R2_b_only_upper   = s["upper",    "R2_b_only"],
    R2_b_plus_X       = s["estimate", "R2_b_plus_X"],
    R2_b_plus_X_lower = s["lower",    "R2_b_plus_X"],
    R2_b_plus_X_upper = s["upper",    "R2_b_plus_X"],
    incremental       = s["estimate", "incremental"],
    incremental_lower = s["lower",    "incremental"],
    incremental_upper = s["upper",    "incremental"]
  )
}
tab <- do.call(rbind, tab)

cat("\n============================================================\n")
cat("R² on the scale of S (fully controlled regime)\n")
cat("mean over 20 repeats; patient-level bootstrap 95% CI (2,000 resamples)\n")
cat("============================================================\n")
print(format(tab, digits = 3), row.names = FALSE)
write.csv(tab, "PTPM_0.1.6_scale_of_S.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# 2. Oldham's test: correlation between change (b - c) and the mean of
#    the baseline and post-treatment scores ((b + c) / 2)
# ---------------------------------------------------------------------
gw  <- if (!is.null(cache[["ptpm_genome_wide.rds"]])) cache[["ptpm_genome_wide.rds"]] else
         readRDS("ptpm_genome_wide.rds")
b   <- gw$baseline_score
cet <- gw$cet_score
change <- b - cet
avg    <- (b + cet) / 2

cat("\n============================================================\n")
cat("Oldham's test (n =", length(b), "paired models,",
    length(unique(gw$groups)), "patients)\n")
cat("============================================================\n")
print(cor.test(change, avg))

# Sensitivity: the 106 models derive from 72 patients, so the interval
# above treats models as independent. Patient-level bootstrap interval:
set.seed(123)
idx_by_patient <- split(seq_along(b), gw$groups)
r_boot <- replicate(2000, {
  i <- unlist(idx_by_patient[sample(names(idx_by_patient), replace = TRUE)],
              use.names = FALSE)
  cor(change[i], avg[i])
})
cat("Patient-level bootstrap 95% CI for Oldham's r (2,000 resamples):",
    sprintf("%.3f to %.3f", quantile(r_boot, 0.025), quantile(r_boot, 0.975)), "\n")

cat("\n")
print(sessionInfo())
