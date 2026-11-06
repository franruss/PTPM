# =============================================================================
# Supplementary File 3 - PTPM manuscript
# Transferability of the baseline KRAS pathway signal to a pharmacological
# endpoint in IRCC-PDX (Perron et al., Nat Commun 2024)
#
# Reproduces: end of Section 5 and Supplementary Note S9 (Tables S9.1, S9.2).
#
# Analysis: multivariable linear regression of tumour shrinkage after three
# weeks of cetuximab on baseline GSVA scores of the three EGFR-axis Hallmark
# pathways (KRAS, E2F, mTORC1). This is not a PTPM cross-validation: IRCC-PDX
# has baseline profiles only, so no paired suppression phenotype is built.
# IRCC-PDX and XENTURION derive from the same Candiolo PDX collection; the
# primary analysis is therefore restricted to patients not used in the
# XENTURION paired discovery analysis.
#
# Input data (Perron et al., doi:10.1038/s41467-024-53163-y)
#   - IRCC-PDX_rnaSeq.txt: RNA-seq counts, from the figshare link given in the
#     Data availability section of Perron et al.
#     (https://figshare.com/s/35d13c7e7cf8f4759334).
#   - cetuximab_response_corrected.txt: Supplementary Data 1, tab "Cetuximab
#     Response", columns A, H and I ("ircc_id", "Volume change in percentage at
#     3 weeks", "Cetuximab_Standard_3wks_cat"), saved as a tab-delimited file.
#   - kras_mutated.txt: Supplementary Data 1, tab "Mutations (110)", saved as a
#     tab-delimited file with one row per sample, including the columns
#     "ircc_id" and "KRAS_mut" (0 = wild-type, 1 = mutant).
#   - XENTURION metadata, loaded with PTPM::load_ptpm_example(), used only to
#     identify the 72 patients of the discovery analysis.
#
# Outputs (written to the working directory)
#   - TableS9_1_primary_nonoverlapping.txt  (upper block of Tables S9.1/S9.2)
#   - TableS9_1_all_patients.txt            (lower block of Tables S9.1/S9.2)
#
# Requirements: R >= 4.5, PTPM 0.1.5, GSVA, msigdbr, edgeR, dplyr, car.
# Runtime: a few minutes.
# =============================================================================

library(PTPM)
library(GSVA)
library(msigdbr)
library(edgeR)
library(dplyr)
library(car)

# -----------------------------------------------------------------------------
# 1. Patients of the XENTURION paired discovery analysis (72)
# -----------------------------------------------------------------------------

pheno_raw <- load_ptpm_example()$metadata
parse_treatment <- function(title) {
  if (grepl("cetux_6w_NT|cetuxi_h_NT|cetux_h_NT", title)) return("NT")
  else if (grepl("cetux|cetuxi", title))                   return("CET")
  else if (grepl("BASALE", title))                         return("BASALE")
  else                                                      return(NA)
}
pheno_raw$treatment_geo <- sapply(as.character(pheno_raw$title), parse_treatment)
pheno_raw$pair_id <- substr(pheno_raw$title, 1, 26)
pheno_raw$pair_id <- gsub("(LMO.*[0-9]{3})[0-9]R", "\\1XR", pheno_raw$pair_id)
pheno_raw$pair_id <- gsub("(LMX.*[0-9]{3})[0-9]{2}T", "\\1XXT", pheno_raw$pair_id)
pheno_raw$pair_id <- gsub("UMR[0-9]+$", "", pheno_raw$pair_id)
pair_table <- table(pheno_raw$pair_id, pheno_raw$treatment_geo)
paired_ids <- rownames(pair_table)[pair_table[, "CET"] > 0 & pair_table[, "NT"] > 0]

disc_ids <- unique(substr(paired_ids, 1, 7))
stopifnot(length(paired_ids) == 106, length(disc_ids) == 72)

# -----------------------------------------------------------------------------
# 2. IRCC-PDX expression: human genes, TMM normalisation, log-CPM
# -----------------------------------------------------------------------------

rnaseq <- read.table("IRCC-PDX_rnaSeq.txt", header = TRUE, sep = "\t", row.names = 1)
counts_human <- rnaseq[grepl("^H_", rownames(rnaseq)), ]
rownames(counts_human) <- sub("^H_", "", rownames(counts_human))

dge <- calcNormFactors(DGEList(counts = as.matrix(counts_human)), method = "TMM")
expr_tmm <- cpm(dge, log = TRUE, prior.count = 1)

# Average samples from the same patient (7-character ID)
core_ids <- substr(colnames(expr_tmm), 1, 7)
expr_avg <- sapply(unique(core_ids), function(id) rowMeans(expr_tmm[, core_ids == id, drop = FALSE]))
cat(sprintf("Expression: %d genes x %d samples -> %d patients\n",
            nrow(expr_tmm), ncol(expr_tmm), ncol(expr_avg)))
# Expected: 56,609 genes, 247 samples, 227 patients with response data (section 3)

# -----------------------------------------------------------------------------
# 3. Response data and alignment
# -----------------------------------------------------------------------------

meta_raw <- read.table("cetuximab_response_corrected.txt", header = TRUE, sep = "\t", row.names = 1)
meta_raw <- meta_raw[, !grepl("^X(\\.\\d+)?$", colnames(meta_raw))]   # drop empty columns, if any
stopifnot(ncol(meta_raw) == 2)
colnames(meta_raw) <- c("vol_change", "response_cat")
meta_raw$ircc_id <- substr(rownames(meta_raw), 1, 7)

common_ids <- intersect(colnames(expr_avg), meta_raw$ircc_id)
cat(sprintf("Patients with RNA-seq and response: %d\n", length(common_ids)))   # expected 227
expr_aligned <- expr_avg[, common_ids]
meta_aligned <- meta_raw[match(common_ids, meta_raw$ircc_id), ]

# -----------------------------------------------------------------------------
# 4. GSVA on the three EGFR-axis Hallmark gene sets
# -----------------------------------------------------------------------------

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
gene_sets_3 <- split(hallmark$gene_symbol, hallmark$gs_name)[c("HALLMARK_KRAS_SIGNALING_UP",
                                                               "HALLMARK_E2F_TARGETS",
                                                               "HALLMARK_MTORC1_SIGNALING")]
gsva_scores_3 <- gsva(gsvaParam(exprData = as.matrix(expr_aligned), geneSets = gene_sets_3))

# Outcome: y = tumour volume change (%); y_sup = -(y / 100) = tumour shrinkage,
# positive = response (same sign convention as the XENTURION suppression score).
df_vol <- data.frame(
  y    = meta_aligned$vol_change,
  KRAS = as.numeric(gsva_scores_3["HALLMARK_KRAS_SIGNALING_UP", common_ids]),
  E2F  = as.numeric(gsva_scores_3["HALLMARK_E2F_TARGETS",       common_ids]),
  MTOR = as.numeric(gsva_scores_3["HALLMARK_MTORC1_SIGNALING",  common_ids]),
  row.names = common_ids
)
stopifnot(!anyNA(df_vol))

# Variance of tumour-volume change explained by mean baseline EGFR-axis activity
baseline_mean <- (df_vol$KRAS + df_vol$E2F + df_vol$MTOR) / 3
cat(sprintf("Variance explained by mean baseline EGFR-axis activity: %.3f\n",
            summary(lm(df_vol$y ~ baseline_mean))$r.squared))
# Expected: 0.036 (compared with 0.383 for the uncorrected paired phenotype in XENTURION)

# -----------------------------------------------------------------------------
# 5. KRAS mutational status
# -----------------------------------------------------------------------------

kras_raw <- read.table("kras_mutated.txt", header = TRUE, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("ircc_id", "KRAS_mut") %in% colnames(kras_raw)))
kras_tbl <- kras_raw %>%
  mutate(core_id = substr(ircc_id, 1, 7)) %>%
  group_by(core_id) %>%
  summarise(KRAS_mut = max(KRAS_mut, na.rm = TRUE), .groups = "drop")
df_vol$KRAS_mut <- kras_tbl$KRAS_mut[match(rownames(df_vol), kras_tbl$core_id)]
stopifnot(!anyNA(df_vol$KRAS_mut))

# -----------------------------------------------------------------------------
# 6. Multivariable model with permutation test, by stratum
# -----------------------------------------------------------------------------

run_multivariate_permtest <- function(df, stratum_label, n_perm = 1000, seed = 42) {
  df$y_sup <- -(df$y / 100)
  fit <- lm(y_sup ~ KRAS + E2F + MTOR, data = df)
  cf  <- summary(fit)$coefficients
  # Partial correlation of each pathway, adjusted for the other two
  partial_r <- cf[-1, "t value"] / sqrt(cf[-1, "t value"]^2 + fit$df.residual)
  # Permutation test on the adjusted R2 of the whole model
  set.seed(seed)
  observed_r2 <- summary(fit)$adj.r.squared
  perm_r2 <- replicate(n_perm, {
    df_perm <- df
    df_perm$y_sup <- sample(df_perm$y_sup)
    summary(lm(y_sup ~ KRAS + E2F + MTOR, data = df_perm))$adj.r.squared
  })
  data.frame(
    Stratum        = stratum_label,
    n              = nrow(df),
    Adj_R2         = round(observed_r2, 4),
    KRAS_beta      = round(cf["KRAS", "Estimate"], 3),
    KRAS_partial_r = round(partial_r["KRAS"], 3),
    KRAS_p         = signif(cf["KRAS", "Pr(>|t|)"], 3),
    E2F_partial_r  = round(partial_r["E2F"], 3),
    E2F_p          = signif(cf["E2F", "Pr(>|t|)"], 3),
    MTOR_partial_r = round(partial_r["MTOR"], 3),
    MTOR_p         = signif(cf["MTOR", "Pr(>|t|)"], 3),
    Perm_emp_p     = signif(mean(perm_r2 >= observed_r2), 3),
    Perm_Z         = round((observed_r2 - mean(perm_r2)) / sd(perm_r2), 2),
    stringsAsFactors = FALSE
  )
}

# Three strata (all, KRAS wild-type, KRAS-mutant), BH adjustment of the KRAS
# P values within each analysis set, and KRAS x mutation interaction test.
analyse_set <- function(df, label) {
  tab <- rbind(
    run_multivariate_permtest(df, sprintf("%s: all (n=%d)", label, nrow(df))),
    run_multivariate_permtest(df[df$KRAS_mut == 0, ],
                              sprintf("%s: KRAS wild-type (n=%d)", label, sum(df$KRAS_mut == 0))),
    run_multivariate_permtest(df[df$KRAS_mut == 1, ],
                              sprintf("%s: KRAS-mutant (n=%d)", label, sum(df$KRAS_mut == 1)))
  )
  tab$KRAS_p_adj <- signif(p.adjust(tab$KRAS_p, method = "BH"), 3)

  df_int <- df
  df_int$y_sup    <- -(df_int$y / 100)
  df_int$KRAS_mut <- factor(df_int$KRAS_mut, levels = c(0, 1), labels = c("WT", "MUT"))
  fit_int <- lm(y_sup ~ KRAS * KRAS_mut + E2F + MTOR, data = df_int)
  ci <- confint(fit_int)["KRAS:KRAS_mutMUT", ]
  cat(sprintf("%s - KRAS x mutation interaction: beta = %+.3f [95%% CI %.3f; %.3f], P = %.3g\n",
              label, coef(fit_int)["KRAS:KRAS_mutMUT"], ci[1], ci[2],
              summary(fit_int)$coefficients["KRAS:KRAS_mutMUT", "Pr(>|t|)"]))
  tab
}

# -----------------------------------------------------------------------------
# 7. Primary analysis: patients not used in the XENTURION discovery set
# -----------------------------------------------------------------------------

overlap_disc <- intersect(rownames(df_vol), disc_ids)
cat(sprintf("IRCC-PDX patients also in the XENTURION discovery set: %d of %d\n",
            length(overlap_disc), nrow(df_vol)))            # expected 61 of 227
df_noov <- df_vol[!rownames(df_vol) %in% overlap_disc, ]
cat(sprintf("KRAS status of the overlapping patients: %d wild-type, %d mutant\n",
            sum(df_vol[overlap_disc, "KRAS_mut"] == 0), sum(df_vol[overlap_disc, "KRAS_mut"] == 1)))

table_primary <- analyse_set(df_noov, "Not used for discovery")
cat("\n=== Tables S9.1/S9.2, upper block (primary analysis) ===\n")
print(table_primary, row.names = FALSE)
# Expected (all, n = 166): KRAS beta -1.575, partial r -0.179, P 0.022, P_adj 0.065, permutation P 0.013

# -----------------------------------------------------------------------------
# 8. All patients (n = 227), including the 61 discovery patients
# -----------------------------------------------------------------------------

table_all <- analyse_set(df_vol, "All patients")
cat("\n=== Tables S9.1/S9.2, lower block (all patients) ===\n")
print(table_all, row.names = FALSE)
# Expected (all, n = 227): KRAS beta -2.084, partial r -0.230, P_adj 0.0015, permutation P < 0.001

# Same model on the original scale (tumour volume change, %): beta of KRAS
fit_vol <- lm(y ~ KRAS + E2F + MTOR, data = df_vol)
print(round(summary(fit_vol)$coefficients, 4))   # expected KRAS beta = +208
print(round(vif(fit_vol), 3))                     # all VIF < 3.1

write.table(table_primary, "TableS9_1_primary_nonoverlapping.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(table_all,     "TableS9_1_all_patients.txt",           sep = "\t", quote = FALSE, row.names = FALSE)

sessionInfo()
