# =============================================================================
# Supplementary File 2 - PTPM manuscript
# Validation of the composite suppression phenotype (CSP) against independent
# cetuximab response endpoints in XENTURION
#
# Reproduces: Section 5, Table 2, Supplementary Note S8 and Figure S8 of the
# manuscript, including the placebo specificity control.
#
# Input data
#   - XENTURION expression and metadata (GEO GSE204805), loaded with
#     PTPM::load_ptpm_example().
#   - Cetuximab_response_file_Xenturion.txt: Leto et al., Nat Commun 2024
#     (doi:10.1038/s41467-024-51909-2), Supplementary Data 6, first tab, saved
#     as a tab-delimited file. First column: patient ID (7 characters, e.g.
#     CRC0030). The next nine columns must be, in this order:
#       1. ATP ratio (cetuximab/untreated), 5,000 cells/well
#       2. ATP ratio (cetuximab/untreated), 1,250 cells/well
#       3. ATP ratio (cetuximab/untreated), 20,000 cells/well
#       4. Relative total area (day 7/day 0), 5,000 cells/well
#       5. Relative total area (day 7/day 0), 1,250 cells/well
#       6. Relative total area (day 7/day 0), 20,000 cells/well
#       7. PDX tumour volume change (%), cetuximab arm
#       8. PDX tumour volume change (%), placebo arm
#       9. Tumour Growth Inhibition (TGI) score
#     The columns are renamed by position (section 2), so their order matters.
#
# Outputs (written to the working directory)
#   - Table2_pearson_corrected.txt, Table2_spearman_corrected.txt
#     (Table 2: coupling-corrected CSP)
#   - Additional_pearson_naive.txt, Additional_spearman_naive.txt
#     (uncorrected CSP, reported in the text of Section 5)
#   - FigureS8_pharmacological_validation_corrected.pdf / .tif (Figure S8)
#   - Additional_figure_pharmacological_validation_naive.pdf / .tif
#
# Requirements: R >= 4.5, PTPM 0.1.5, GSVA, msigdbr, dplyr, ggplot2, ggpubr,
# patchwork. Runtime: a few minutes.
# =============================================================================

library(PTPM)
library(GSVA)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. XENTURION: paired models and composite suppression phenotype (CSP)
# -----------------------------------------------------------------------------

example_data <- load_ptpm_example()
rownames(example_data$expression) <- example_data$expression[, 1]   # gene identifiers
example_data$expression <- as.matrix(example_data$expression[, -1])
mode(example_data$expression) <- "numeric"

pheno_raw <- example_data$metadata

# Treatment label from the GEO title. "BASALE" samples are xenografts or
# organoids grown in normal conditions; they are not the untreated controls of
# the cetuximab experiments and are not used.
parse_treatment <- function(title) {
  if (grepl("cetux_6w_NT|cetuxi_h_NT|cetux_h_NT", title)) return("NT")
  else if (grepl("cetux|cetuxi", title))                   return("CET")
  else if (grepl("BASALE", title))                         return("BASALE")
  else                                                      return(NA)
}
pheno_raw$treatment_geo <- sapply(as.character(pheno_raw$title), parse_treatment)

# Model identifier shared by the NT and CET samples of the same model
pheno_raw$entire_codes <- substr(pheno_raw$title, 1, 26)
pheno_raw$pair_id <- pheno_raw$entire_codes
pheno_raw$pair_id <- gsub("(LMO.*[0-9]{3})[0-9]R", "\\1XR", pheno_raw$pair_id)
pheno_raw$pair_id <- gsub("(LMX.*[0-9]{3})[0-9]{2}T", "\\1XXT", pheno_raw$pair_id)
pheno_raw$pair_id <- gsub("UMR[0-9]+$", "", pheno_raw$pair_id)

pair_table <- table(pheno_raw$pair_id, pheno_raw$treatment_geo)
paired_ids <- rownames(pair_table)[pair_table[, "CET"] > 0 & pair_table[, "NT"] > 0]
stopifnot(length(paired_ids) == 106)                         # 78 LMX + 28 LMO
stopifnot(length(unique(substr(paired_ids, 1, 7))) == 72)    # 72 patients

# Average replicate samples per model and condition
pheno_paired <- pheno_raw[pheno_raw$pair_id %in% paired_ids, ]
pheno_paired$group <- paste(pheno_paired$pair_id, pheno_paired$treatment_geo, sep = "_")
group_levels <- unique(pheno_paired$group)
expr_avg <- sapply(group_levels, function(g) {
  cols <- pheno_paired$entire_codes[pheno_paired$group == g]
  rowMeans(example_data$expression[, cols, drop = FALSE])
})

# GSVA on the three EGFR-axis Hallmark gene sets, human genes only
hallmark  <- msigdbr(species = "Homo sapiens", collection = "H")
egfr_sets <- hallmark[hallmark$gs_name %in% c("HALLMARK_KRAS_SIGNALING_UP",
                                              "HALLMARK_MTORC1_SIGNALING",
                                              "HALLMARK_E2F_TARGETS"), ]
gene_sets <- split(egfr_sets$gene_symbol, egfr_sets$gs_name)

expr_human <- expr_avg[grepl("^H_", rownames(expr_avg)), ]
rownames(expr_human) <- sub("^H_", "", rownames(expr_human))
gsva_scores <- gsva(gsvaParam(expr_human, gene_sets))

gsva_nt  <- gsva_scores[, grepl("_NT$",  colnames(expr_avg))]
gsva_cet <- gsva_scores[, grepl("_CET$", colnames(expr_avg))]
pair_ids <- intersect(sub("_NT$", "", colnames(gsva_nt)), sub("_CET$", "", colnames(gsva_cet)))

# Suppression = NT - CET (positive = reduced pathway activity);
# CSP = unweighted mean of the three pathway-specific suppression scores
suppression_matrix <- sapply(pair_ids, function(pid) {
  gsva_scores[, paste0(pid, "_NT")] - gsva_scores[, paste0(pid, "_CET")]
})
suppression_score_naive <- colMeans(suppression_matrix)
print(summary(suppression_score_naive))
# Expected: median 0.173, mean 0.161 (Supplementary Note S4)

# -----------------------------------------------------------------------------
# 2. Coupling-corrected CSP
# -----------------------------------------------------------------------------
# Residual of the CSP regressed on the mean baseline (NT) score of the three
# pathways, estimated on all 106 models (descriptive analysis, no prediction).

baseline_score <- colMeans(gsva_nt)
names(baseline_score) <- sub("_NT$", "", names(baseline_score))
baseline_aligned <- baseline_score[names(suppression_score_naive)]
stopifnot(!anyNA(baseline_aligned))

m0 <- lm(suppression_score_naive ~ baseline_aligned)
suppression_score_corrected <- resid(m0)
names(suppression_score_corrected) <- names(suppression_score_naive)

cat(sprintf("Coupling model: R2 = %.3f (variance of the CSP explained by baseline)\n",
            summary(m0)$r.squared))
# Expected: 0.383

# -----------------------------------------------------------------------------
# 3. Pharmacological endpoints
# -----------------------------------------------------------------------------

pharma <- read.table("Cetuximab_response_file_Xenturion.txt",
                     row.names = 1, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
stopifnot(ncol(pharma) == 9)
colnames(pharma) <- c("ATP_5000", "ATP_1250", "ATP_20000",
                      "Area_5000", "Area_1250", "Area_20000",
                      "PDX_vol_cetux", "PDX_vol_placebo", "TGI")
pharma$model <- rownames(pharma)

# Check that patient IDs match before merging
overlap <- intersect(rownames(pharma), unique(substr(names(suppression_score_naive), 1, 7)))
cat(sprintf("Patients with both CSP and pharmacological data: %d\n", length(overlap)))
if (length(overlap) < 10) {
  stop("Too few matching patient IDs: check the ID format of the pharmacological file ",
       "(expected 7-character patient IDs, e.g. CRC0030).")
}

# Patient-level aggregation: models of the same patient are averaged, because
# the pharmacological endpoints are reported per patient.
build_merged_df <- function(score_vec) {
  data.frame(model = substr(names(score_vec), 1, 7),
             suppression_score = as.numeric(score_vec)) %>%
    group_by(model) %>%
    summarise(suppression_score = mean(suppression_score, na.rm = TRUE), .groups = "drop") %>%
    inner_join(pharma, by = "model")
}
df_naive     <- build_merged_df(suppression_score_naive)
df_corrected <- build_merged_df(suppression_score_corrected)

# -----------------------------------------------------------------------------
# 4. Correlations with the six endpoints (Table 2)
# -----------------------------------------------------------------------------

endpoints <- c("ATP_5000", "ATP_1250", "ATP_20000", "TGI", "PDX_vol_cetux", "Area_5000")

compute_correlations <- function(df, endpoints, method = "pearson") {
  out <- do.call(rbind, lapply(endpoints, function(ep) {
    d  <- df %>% filter(!is.na(.data[[ep]]))
    ct <- cor.test(d$suppression_score, d[[ep]], method = method,
                   exact = if (method == "spearman") FALSE else NULL)
    row <- data.frame(Endpoint = ep, n = nrow(d),
                      estimate = round(as.numeric(ct$estimate), 3),
                      p_value  = signif(ct$p.value, 3))
    if (method == "pearson") {
      row$CI_low  <- round(ct$conf.int[1], 3)
      row$CI_high <- round(ct$conf.int[2], 3)
    }
    row
  }))
  # Benjamini-Hochberg across the six endpoints, separately for each method
  out$p_adj  <- signif(p.adjust(out$p_value, method = "BH"), 3)
  out$signif <- ifelse(out$p_adj < 0.001, "***",
                ifelse(out$p_adj < 0.01,  "**",
                ifelse(out$p_adj < 0.05,  "*", "ns")))
  colnames(out)[colnames(out) == "estimate"] <- if (method == "pearson") "r" else "rho"
  out
}

corr_corrected_pearson  <- compute_correlations(df_corrected, endpoints, "pearson")
corr_corrected_spearman <- compute_correlations(df_corrected, endpoints, "spearman")
corr_naive_pearson      <- compute_correlations(df_naive,     endpoints, "pearson")
corr_naive_spearman     <- compute_correlations(df_naive,     endpoints, "spearman")

cat("\n=== Table 2: coupling-corrected CSP, Pearson ===\n");  print(corr_corrected_pearson)
cat("\n=== Table 2: coupling-corrected CSP, Spearman ===\n"); print(corr_corrected_spearman)
cat("\n=== Uncorrected CSP, Pearson (Section 5, text) ===\n"); print(corr_naive_pearson)
cat("\n=== Uncorrected CSP, Spearman (Section 5, text) ===\n"); print(corr_naive_spearman)

write.table(corr_corrected_pearson,  "Table2_pearson_corrected.txt",  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(corr_corrected_spearman, "Table2_spearman_corrected.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(corr_naive_pearson,      "Additional_pearson_naive.txt",  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(corr_naive_spearman,     "Additional_spearman_naive.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Specificity control: placebo arm
# -----------------------------------------------------------------------------
# The CSP should not correlate with tumour-volume change in placebo-treated
# PDXs; a correlation there would indicate a link with tumour growth in
# general rather than with the response to cetuximab.

d <- df_naive %>% filter(!is.na(PDX_vol_cetux), !is.na(PDX_vol_placebo))
ct_cetux   <- cor.test(d$suppression_score, d$PDX_vol_cetux,   method = "pearson")
ct_placebo <- cor.test(d$suppression_score, d$PDX_vol_placebo, method = "pearson")
cat(sprintf("\nMatched subset: n = %d models\n", nrow(d)))
cat(sprintf("Cetuximab arm: r = %+.3f, p = %.3g\n", ct_cetux$estimate,   ct_cetux$p.value))
cat(sprintf("Placebo arm:   r = %+.3f, p = %.3g\n", ct_placebo$estimate, ct_placebo$p.value))
# Expected: cetuximab r = -0.80, p < 0.001; placebo r = -0.12, p = 0.67 (n = 16)

# -----------------------------------------------------------------------------
# 6. Figures (Figure S8: coupling-corrected CSP)
# -----------------------------------------------------------------------------

theme_pub <- theme_classic(base_size = 11) +
  theme(axis.title    = element_text(size = 10),
        axis.text     = element_text(size = 9),
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))

scatter_cor <- function(data, y_var, y_lab, subtitle) {
  d <- data %>% filter(!is.na(.data[[y_var]]))
  ggplot(d, aes(x = suppression_score, y = .data[[y_var]])) +
    geom_point(size = 2, alpha = 0.7, color = "#2166AC") +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B", fill = "#FDDBC7", linewidth = 0.8) +
    stat_cor(method = "pearson",
             aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
             label.x.npc = "left", label.y.npc = "top", size = 3.2, color = "black") +
    labs(subtitle = subtitle, x = "PTPM suppression score", y = y_lab) +
    theme_pub
}

make_figure <- function(df_in, main_title = "") {
  pA <- scatter_cor(df_in, "ATP_5000",
                    "ATP ratio (cetuximab / untreated)\n5,000 cells/well", "Ex vivo PDXT response")
  pB <- scatter_cor(df_in, "TGI", "Tumor Growth Inhibition (TGI)", "In vivo PDX response")
  pC <- scatter_cor(df_in, "PDX_vol_cetux",
                    "Tumor volume change (%)\n3 weeks cetuximab", "In vivo PDX tumor volume")
  pD <- scatter_cor(df_in, "Area_5000",
                    "Relative total area (day 7 / day 0)\ncetuximab / untreated", "Ex vivo PDXT longitudinal")
  (pA | pB) / (pC | pD) +
    plot_annotation(title = main_title, tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
                    theme = theme(plot.title = element_text(size = 12, face = "bold"),
                                  plot.tag   = element_text(size = 11, face = "bold")))
}

fig_s8    <- make_figure(df_corrected)
fig_naive <- make_figure(df_naive, "Uncorrected PTPM suppression score")

ggsave("FigureS8_pharmacological_validation_corrected.pdf", fig_s8,
       width = 210, height = 140, units = "mm", device = cairo_pdf)
ggsave("FigureS8_pharmacological_validation_corrected.tif", fig_s8,
       width = 210, height = 140, units = "mm", dpi = 600, device = "tiff", compression = "lzw")
ggsave("Additional_figure_pharmacological_validation_naive.pdf", fig_naive,
       width = 170, height = 140, units = "mm", device = cairo_pdf)
ggsave("Additional_figure_pharmacological_validation_naive.tif", fig_naive,
       width = 170, height = 140, units = "mm", dpi = 600, device = "tiff", compression = "lzw")

sessionInfo()
