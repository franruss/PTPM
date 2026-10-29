# User Guide

This guide reproduces the PTPM analyses reported in the manuscript with the
installed `PTPM` package. Six functions apply the same protocol to different
predictor or outcome models:

| Function | Predictors | Outcome | Feature selection (K) | Default α grid |
|---|---|---|---|---|
| `run_ptpm()` | genome-wide baseline expression | mean pathway suppression S | top-K genes: 3000, 300 | 0–0.1 |
| `run_ptpm_egfr()` | a priori EGFR-axis genes (union of the Hallmark KRAS_SIGNALING_UP, MTORC1_SIGNALING and E2F_TARGETS gene sets; 569 genes) | mean pathway suppression S | top-K genes: 50, 100, 200, 300, 400, 500, all | 0–0.1 |
| `run_ptpm_pathways()` | baseline GSVA scores of the 47 non-outcome Hallmark pathways | mean pathway suppression S | top-K pathways: 10, 20, 30, 40, all | 0–1 |
| `run_ptpm_weighted_score()` | EGFR-axis genes (as `run_ptpm_egfr()`) | PC1-weighted pathway suppression, with the PCA fitted within training folds | top-K genes: 50, 100, 200, 300, 400, 500, all | 0–0.4 |
| `run_ptpm_hybrid()` | 47 non-outcome pathways (always included) + top-K EGFR-axis genes | mean pathway suppression S | top-K genes: 100, 200, 300, 400, all | 0–0.2 |
| `run_ptpm_hybrid_interactions()` | 3 non-outcome pathways (G2M_CHECKPOINT, GLYCOLYSIS, MYC_TARGETS_V1) and their pairwise interactions (always included) + top-K EGFR-axis genes | mean pathway suppression S | top-K genes: 100, 200, 300, 400, all | 0–0.2 |

All six functions share the same protocol:
- **Phenotype:** pathway suppression S = baseline (NT) minus cetuximab (CET) pathway activity, from GSVA scores of the three outcome Hallmark sets.
- **Cross-validation:** nested and patient-grouped, with 3 outer folds × 20 repeats. The inner 5-fold tuning CV is also grouped by patient.
- **2 × 2 design:**
  - feature selection on the full dataset (standard, coupling-adjusted) or within training folds (fold-internal, fully controlled);
  - with or without coupling adjustment.
- **Coupling adjustment:** the outcome is regressed on the baseline pathway score, and the regression is re-estimated within every outer and inner training fold.
- **Metric:** R² = mean over the 20 repeats of 1 − SSE/SST on each repeat's out-of-fold predictions.
- **Uncertainty:** patient-level bootstrap 95% confidence intervals, plus corrected resampled t intervals (Nadeau & Bengio) as a sensitivity analysis.
- **Diagnostics and sensitivity analyses:**
  - coupling null and partial correlations;
  - training-reference GSVA;
  - unnormalised ssGSEA;
  - pathway-specific adjustment;
  - KRAS-free phenotype;
  - exclusion of KRAS genes or of high-overlap predictor pathways.

A seventh function, `run_ptpm_benchmark()`, re-evaluates the genome-wide model
of `run_ptpm()` with three alternative learners (XGBoost, sparse PLS and linear
SVR) under the same protocol and on the same folds (section 8).

For K equal to the whole predictor pool, no feature selection is performed. The
standard and fold-internal regimes (and the coupling-adjusted and fully
controlled regimes) therefore coincide by construction. Residual differences of
about 0.001 reflect the numerical tolerance of the solver.

## 1. Install

From the directory containing `PTPM_0.1.5.tar.gz`:

```r
install.packages("PTPM_0.1.5.tar.gz", repos = NULL, type = "source")
```

Required packages: `glmnet`, `caret`, `GSVA`, `msigdbr`, `foreach`, `doParallel`.
`msigdbr` is also used to retrieve all 50 Hallmark gene sets for
`run_ptpm_pathways()`, `run_ptpm_hybrid()` and `run_ptpm_hybrid_interactions()`.
`Rcpp` is optional: it compiles a fast kernel for the training-reference GSVA
sensitivity analysis, and a slower pure-R version is used without it.

`GEOquery` (Bioconductor) is needed only to download other GEO series (section 9).

Optional packages for the benchmark of section 8: `xgboost`, `e1071` (CRAN) and
`mixOmics` (Bioconductor). Only the packages of the learners actually requested
are needed:

```r
install.packages(c("xgboost", "e1071"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("mixOmics")
```

## 2. Load

```r
library(PTPM)
example_data <- load_ptpm_example()
```

No `source()` call is required: the functions in the installed package are
loaded automatically by R.

## 3. Quick smoke test (a few minutes per model)

Every function accepts the same arguments. The settings below switch off the
heavier analyses and use only two repeats:

```r
smoke_args <- list(
  expression = example_data$expression,
  metadata = example_data$metadata,
  repeats = 2,
  n_boot = 200,
  gsva_sensitivity = FALSE,
  ssgsea_sensitivity = FALSE,
  robustness = FALSE,
  n_perm_coupling = 100,
  seed = 123,
  n_cores = parallel::detectCores() - 1
)
smoke <- do.call(run_ptpm_egfr, smoke_args)
smoke$results_genome_wide
```

## 4. Full analysis (manuscript settings)

The manuscript settings are `k = 3`, `repeats = 20`, `n_boot = 2000` and
`seed = 123`, with each function's default K grid. For `run_ptpm()` the K values
are `top_k_values = c(3000, 300)`.

```r
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
```

Each model can also be run on its own, for example:

```r
result <- run_ptpm_hybrid(
  expression = example_data$expression,
  metadata = example_data$metadata,
  k = 3, repeats = 20, n_boot = 2000, seed = 123,
  n_cores = parallel::detectCores() - 1
)
```

## 5. Inspect the principal results

The output has the same structure for all six models. The table
`results_genome_wide` keeps this name for every model.

```r
result <- results$egfr_axis       # any of the six models

result$results_genome_wide        # R² of the four regimes for each K
result$confidence_intervals       # R² with patient-bootstrap 95% CI
result$contrast_intervals         # paired differences between regimes, with 95% CI
result$leakage_audit              # where each data-dependent step is estimated
result$definitions                # definitions of "fully controlled" and of R²
```

A compact comparison of the fully controlled R² across models:

```r
do.call(rbind, lapply(names(results), function(nm) {
  ci <- results[[nm]]$confidence_intervals
  ci <- ci[ci$regime == "fully_controlled", c("K", "R2", "boot_lower", "boot_upper")]
  cbind(model = nm, ci)
}))
```

## 6. Diagnostics and sensitivity analyses

```r
result$coupling_null            # cor(b, S) vs permutation of the NT/CET pairing
result$gsva_ref_validation      # check that the reference GSVA reproduces GSVA::gsva()
result$gsva_sensitivity         # cohort-level vs training-reference GSVA
result$gsva_score_agreement
result$robustness_r2            # pathway-specific adjustment, KRAS-free phenotype and
                                # exclusion analyses (see below)
result$robustness_partial_cor   # baseline KRAS score vs phenotype (partial correlations)
result$ssgsea_results           # unnormalised ssGSEA instead of GSVA
result$ssgsea_contrasts
result$ssgsea_coupling_null
result$ssgsea_partial_cor
```

The exclusion analyses in `robustness_r2` depend on the predictors:

| Model | Exclusion analyses |
|---|---|
| `run_ptpm()`, `run_ptpm_egfr()`, `run_ptpm_weighted_score()` | KRAS_SIGNALING_UP genes excluded from the predictors |
| `run_ptpm_pathways()` | predictor pathways sharing more than 20% of genes with the outcome sets excluded (`overlap_threshold`) |
| `run_ptpm_hybrid()`, `run_ptpm_hybrid_interactions()` | both; interactions involving an excluded pathway are also removed |

Model-specific elements:

```r
results$egfr_axis$predictor_pool                    # EGFR-axis genes used as predictors
results$pathways$predictor_pathways                 # the 47 predictor pathways
results$pathways$pathway_overlap                    # gene overlap of each predictor pathway with the outcome sets
results$weighted_score$outcome_weights_global       # whole-cohort PC1 weights (descriptive only)
results$weighted_score$coupling_score               # coupling covariate ("weighted" or "mean")
results$hybrid$predictor_pathways                   # pathways always included
results$hybrid$predictor_gene_pool                  # EGFR-axis gene pool
results$hybrid_interaction$predictor_pathways       # pathways and interaction terms always included
results$hybrid_interaction$pathway_overlap
```

## 7. Inspect predictions

Out-of-fold predictions are stored per regime and per outer fold in
`result$fold_results`. Each element contains the fold index, the predictions
(one column per K) and the observed outcome.

In the coupling-adjusted regimes, the observed outcome is the training-fold
residual. In `run_ptpm_weighted_score()` it is the fold-specific PC1 score:
weights are estimated on each training fold, so outcomes of different folds
are on standardised but not identical scales.

```r
result <- results$egfr_axis
regime <- "fully_controlled"   # or "standard", "coupling_adjusted", "fold_internal"
fr <- result$fold_results[[regime]]
n  <- length(result$y)
preds <- do.call(rbind, lapply(fr, function(z) {
  test_idx <- setdiff(seq_len(n), result$folds[[z$fold]])
  data.frame(
    model     = names(result$y)[test_idx],
    repeat_id = as.integer(sub(".*_Rep", "", names(result$folds)[z$fold])),
    observed  = z$obs,
    setNames(as.data.frame(z$pred), paste0("pred_K", result$top_k_values))
  )
}))
head(preds)
```

## 8. Benchmark against other learners

`run_ptpm_benchmark()` takes the output of `run_ptpm()` and evaluates three
alternative learners on the same genome-wide baseline predictors and the same
phenotype S:

| Learner | Package | Hyperparameters tuned by the inner CV | Scaling |
|---|---|---|---|
| XGBoost | `xgboost` | `nrounds`, `max_depth` | not needed (trees) |
| Sparse PLS (sPLS) | `mixOmics` | `keepX` (≤ K, same for all components), `ncomp` | within each fit |
| Linear SVR | `e1071` | `cost`, `epsilon` (on standardised y) | X and y standardised within each fit |

The protocol is the one of `run_ptpm()`, with no re-estimation of the folds:
- **Folds:** the patient-grouped outer folds and inner tuning folds stored in the
  `run_ptpm()` output are reused unchanged.
- **2 × 2 design:** the same four regimes, with feature selection by baseline
  variance on the full dataset or within outer training folds.
- **Coupling adjustment:** re-estimated within every outer and inner training
  fold, as in `run_ptpm()`.
- **Tuning:** each learner's hyperparameters are chosen only by the
  patient-grouped inner CV (lowest mean squared error over the inner-validation
  models), then refitted on the whole outer training fold.
- **Metric and uncertainty:** the same R², the same patient-level bootstrap
  (same resampled patients as `run_ptpm()`) and the same corrected resampled t
  intervals.
- **Elastic-net reference:** taken from the `run_ptpm()` output, not refitted.
  Every comparison with the elastic net is therefore paired: same folds, same
  repeats, same resampled patients.

### 8.1 Quick check (a few minutes)

A reduced `run_ptpm()` run followed by the benchmark in quick mode:

```r
t0 <- Sys.time()
quick_ptpm <- run_ptpm(
  expression = example_data$expression,
  metadata = example_data$metadata,
  top_k_values = c(3000, 300),
  repeats = 2,
  n_boot = 200,
  gsva_sensitivity = FALSE,
  ssgsea_sensitivity = FALSE,
  robustness = FALSE,
  n_perm_coupling = 100,
  seed = 123,
  n_cores = parallel::detectCores() - 1
)
quick_bench <- run_ptpm_benchmark(quick_ptpm, quick = TRUE)
print(Sys.time() - t0)

quick_bench$results_table       # R² of the four regimes, one row per learner x K
quick_bench$method_comparison   # learner minus elastic net, paired, with 95% CI
```

`quick = TRUE` uses at most 2 repeats, smaller tuning grids (XGBoost `nrounds`
50/150/400; sPLS `keepX` 10/50/200; SVR `cost` 0.01/1/100) and 200 bootstrap
resamples. Any of these settings can still be given explicitly. Folds, regimes
and metric are unchanged, so only the precision is reduced: the values are
indicative and the intervals wide.

To make the check even faster:
- `methods = "svr"`: SVR is the fastest learner and XGBoost the slowest;
- `top_k_values = 300`: a single K instead of two.

```r
run_ptpm_benchmark(quick_ptpm, quick = TRUE, methods = "svr", top_k_values = 300)
```

### 8.2 Full benchmark (manuscript settings)

The benchmark is applied to the `run_ptpm()` result of section 4 (20 repeats),
with the full tuning grids and 2000 bootstrap resamples:

```r
genome_wide <- readRDS("ptpm_genome_wide.rds")   # saved in section 4
bench <- run_ptpm_benchmark(genome_wide, n_cores = parallel::detectCores() - 1)
saveRDS(bench, "ptpm_benchmark.rds")
```

The seed defaults to the one of the `run_ptpm()` run. The number of repeats
can be limited with `repeats` (the first repeats are used, on the same folds).

### 8.3 Inspect the benchmark

```r
bench$results_table           # R² of the four regimes, one row per learner x K
                              # (elastic net included as reference)
bench$confidence_intervals    # R² with patient-bootstrap 95% CI
bench$contrast_intervals      # factorial decomposition (feature selection vs
                              # coupling) for each learner
bench$method_comparison       # learner minus elastic net, paired, with 95% CI
bench$tuning                  # hyperparameters chosen in each outer fold
bench$leakage_audit           # where each data-dependent step is estimated
bench$results_summary         # includes cor2_pooled_legacy, the pooled squared
                              # correlation used by earlier versions
```

A compact comparison of the fully controlled R² across learners:

```r
ci <- bench$confidence_intervals
ci[ci$regime == "fully_controlled", c("method", "K", "R2", "boot_lower", "boot_upper")]
```

## 9. Other paired datasets

All seven functions can be applied to any dataset with a baseline and a
treated expression profile for each paired unit (model, tumour or patient).
The samples are described by a `design` table instead of the XENTURION
`metadata`:

| Column | Content |
|---|---|
| `sample` | expression column name |
| `pair` | paired unit (one baseline and one treated profile; replicates are averaged) |
| `condition` | `"baseline"` or `"treated"` (other values are ignored; labels can be changed with `baseline_label` / `treated_label`) |
| `patient` | optional; grouping of the cross-validation (defaults to `pair`) |

Expression must be on a log-like scale with gene symbols as identifiers.
Pairs without both profiles are excluded with a message. With `design`, no
gene prefix is expected (`gene_prefix = ""`); the `"H_"` prefix applies only
to the XENTURION data. With `design = NULL` the XENTURION parsing of
`metadata$title` is used, and results are identical to PTPM 0.1.4.

### 9.1 External example: GSE109756 (cetuximab, head and neck cancer)

[GEO] https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109756
GSE109756 (CHIRON window-of-opportunity study) profiles tumour biopsies taken
at baseline and at surgery after about two weeks of cetuximab (Affymetrix
HG-U133 Plus 2.0). `load_ptpm_geo("GSE109756")` (or its shortcut
`load_ptpm_gse109756()`) downloads the series with `GEOquery`, keeps one probe
per gene symbol and builds the `design` table.
The patients without cetuximab (control arm), if present, are returned
separately in `design_control`.

```r
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")

geo <- load_ptpm_geo("GSE109756")
geo$samples   # check that patient and time point of every sample are correct
geo$design

res_hnscc <- run_ptpm(
  expression = geo$expression,
  design = geo$design,
  top_k_values = c(300),
  k = 3, repeats = 2, n_boot = 2000, seed = 123,
  n_cores = parallel::detectCores() - 1
)
res_hnscc$results_genome_wide
res_hnscc$contrast_intervals
```

The examples of this section use `repeats = 2` and a single K for a quick
check; the manuscript settings are `repeats = 20` and `top_k_values = c(3000, 300)`.
With `design`, the outer folds are balanced by patient (`balanced_folds = TRUE`),
so that every repeat has exactly k folds of equal size also in small cohorts.

The cohort is small (about 15 complete pairs, so about 5 patients per test
fold): the analysis shows that the protocol transfers to an independent
cohort, tumour type and platform, but the intervals are wide. A control arm
of a few patients is too small for nested cross-validation; it can be used
descriptively (S close to 0 without the drug).

On microarray data, genes expressed in only one or two tumours (for example
cancer-testis antigens) can be selected by the whole-cohort variance ranking
while being almost constant in a training fold. The regimes with feature
selection on the full dataset then give unstable and strongly negative R².
The regimes with feature selection inside the folds are not affected. A
standard unsupervised expression filter, which does not use S, removes these
genes:

```r
ex   <- as.matrix(geo$expression[, -1])
keep <- rowSums(ex > quantile(ex, 0.25)) >= 0.2 * ncol(ex)   # expressed in >= 20% of arrays
res_hnscc_f <- run_ptpm(expression = geo$expression[keep, ], design = geo$design,
                        top_k_values = c(300), k = 3, repeats = 2, n_boot = 2000,
                        seed = 123, n_cores = parallel::detectCores() - 1)
```

If the automatic detection of the annotation fields fails, the function
stops and prints the available fields, which can then be given explicitly,
for example:

```r
geo <- load_ptpm_geo("GSE109756",
                     patient_field = "sample id:ch1",   # patient identifier in GSE109756
                     time_field = "treatment:ch1")      # "before treatment" / "after treatment"
```

### 9.2 Other treatments

The default outcome gene sets are those of the anti-EGFR analysis
(KRAS_SIGNALING_UP, MTORC1_SIGNALING, E2F_TARGETS). For another treatment,
pass the relevant gene sets and either name the set used by the KRAS controls
(`kras_set`) or switch the robustness block off (`robustness = FALSE`).
S is always baseline minus treated score, so it is positive when the
treatment suppresses the outcome gene sets.

[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20181)
GSE20181 (neoadjuvant letrozole in breast cancer; biopsies before treatment,
after 10-14 days and after about 3 months) identifies the patient only by the
number at the start of the title ("10A" = patient 10, pretreatment; "10B" =
patient 10, 10-14 days). The early time point is selected explicitly:

```r
geo_er <- load_ptpm_geo("GSE20181", treated_pattern = "10-14 days")
geo_er$samples   # check the classification; 3-month samples are left out
```

The patient and the time point are detected automatically; only the treated
time point has to be chosen, because the series has two.

GEO sample annotations are free text, so automatic detection cannot cover
every series. When it fails, the `design` table can be written directly from
the annotation, which is always possible. `load_ptpm_geo(gse, annotation_only
= TRUE)` returns the expression matrix and the sample annotation without
building the design. The same GSE20181 design written by hand:

```r
raw <- load_ptpm_geo("GSE20181", annotation_only = TRUE)   # expression + annotation
ann <- raw$samples
my_design <- data.frame(
  sample    = ann$sample,
  pair      = sub("^(\\d+)[A-Z].*", "\\1", ann$title),          # patient number
  condition = ifelse(grepl("pretreatment", ann$title), "baseline",
              ifelse(grepl("10-14 days", ann$title), "treated", NA))   # 3-month samples left out
)
my_design <- my_design[!is.na(my_design$condition), ]
table(my_design$condition)
```

The outcome gene sets for endocrine therapy:

```r
h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")  # category = "H" in msigdbr < 10
er_sets <- split(h$gene_symbol, h$gs_name)[c("HALLMARK_ESTROGEN_RESPONSE_EARLY",
                                             "HALLMARK_ESTROGEN_RESPONSE_LATE",
                                             "HALLMARK_E2F_TARGETS")]
res <- run_ptpm(expression = raw$expression, design = my_design, gene_sets = er_sets, 
  robustness = FALSE, top_k_values = c(300),
  k = 3, repeats = 2, n_boot = 2000, seed = 123,
  n_cores = parallel::detectCores() - 1)
```

### 9.3 Several models per patient: GSE108277 (cetuximab, CRC PDX)

GSE108277 ([GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108277))
profiles colorectal cancer PDX models untreated (NT) and after 24, 48 and
72 hours of cetuximab (Illumina HumanHT-12 v4). A patient can have more than
one PDX model, so the paired unit is the model and the cross-validation is
grouped by patient, as in the XENTURION analysis. The design is written by
hand; the 72-hour time point is used here.

The models belong to the Candiolo PDX collection, which is also the source of
the XENTURION data, so this series is a sensitivity analysis on another
platform and time point rather than an independent validation.

```r
raw <- load_ptpm_geo("GSE108277", annotation_only = TRUE)   # expression + annotation
s <- raw$samples

table(s$`treatment:ch1`)
table(sub("^\\S+\\s*", "", s$title))      # title suffixes: "", PHI, 24H, 48H, 72H

design <- data.frame(
  sample    = s$sample,
  pair      = substr(s$title, 1, 12),                  # PDX model, e.g. CRC0096LMX0a
  condition = ifelse(s$`treatment:ch1` == "NT", "baseline",
              ifelse(grepl("72H", s$title), "treated", NA)),
  patient   = s$`patient id:ch1`                       # CV grouped by patient
)
design <- design[!is.na(design$condition), ]
table(design$condition)

res_108 <- run_ptpm(expression = raw$expression, design = design,
                    top_k_values = c(300), k = 3, repeats = 2, n_boot = 2000,
                    seed = 123, n_cores = parallel::detectCores() - 1)
```

### 9.4 Settings used in the manuscript

The additional paired datasets reported in Supplementary Note S13 were analysed
with the commands of sections 9.1–9.3 and the manuscript settings
(`repeats = 20`, `top_k_values = c(3000, 300)`). Because all three series are
microarrays, the unsupervised expression filter of section 9.1 was applied to
each of them before analysis; the filter does not use S. For GSE108277, whose
models belong to the Candiolo collection, the overlap with the patients of the
XENTURION discovery analysis is also reported.

```r
# Unsupervised expression filter (section 9.1), applied to all three series
expr_filter <- function(expr) {
  ex   <- as.matrix(expr[, -1])
  keep <- rowSums(ex > quantile(ex, 0.25)) >= 0.2 * ncol(ex)   # expressed in >= 20% of arrays
  expr[keep, ]
}
```

GSE109756 (cetuximab, head and neck cancer):

```r
geo_hn <- load_ptpm_geo("GSE109756")
geo_hn$samples   # check that patient and time point of every sample are correct

res_hnscc <- run_ptpm(
  expression = expr_filter(geo_hn$expression),
  design = geo_hn$design,
  top_k_values = c(3000, 300),
  k = 3, repeats = 20, n_boot = 2000, seed = 123,
  n_cores = parallel::detectCores() - 1
)
res_hnscc$results_genome_wide
res_hnscc$confidence_intervals
res_hnscc$contrast_intervals
res_hnscc$coupling_null
```

GSE20181 (neoadjuvant letrozole, breast cancer; 10–14 days):

```r
raw_er <- load_ptpm_geo("GSE20181", annotation_only = TRUE)   # expression + annotation
ann <- raw_er$samples
design_er <- data.frame(
  sample    = ann$sample,
  pair      = sub("^(\\d+)[A-Z].*", "\\1", ann$title),                   # patient number
  condition = ifelse(grepl("pretreatment", ann$title), "baseline",
              ifelse(grepl("10-14 days", ann$title), "treated", NA))    # 3-month samples left out
)
design_er <- design_er[!is.na(design_er$condition), ]
table(design_er$condition)

h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")  # category = "H" in msigdbr < 10
er_sets <- split(h$gene_symbol, h$gs_name)[c("HALLMARK_ESTROGEN_RESPONSE_EARLY",
                                             "HALLMARK_ESTROGEN_RESPONSE_LATE",
                                             "HALLMARK_E2F_TARGETS")]

# robustness = FALSE: the KRAS controls do not apply to endocrine therapy, so the
# pathway-specific adjustment is not computed for this series
res_er <- run_ptpm(
  expression = expr_filter(raw_er$expression),
  design = design_er, gene_sets = er_sets, robustness = FALSE,
  top_k_values = c(3000, 300),
  k = 3, repeats = 20, n_boot = 2000, seed = 123,
  n_cores = parallel::detectCores() - 1
)
res_er$results_genome_wide
res_er$confidence_intervals
res_er$contrast_intervals
res_er$coupling_null
```

GSE108277 (cetuximab 72 h, CRC PDX, Candiolo collection):

```r
raw_108 <- load_ptpm_geo("GSE108277", annotation_only = TRUE)   # expression + annotation
s <- raw_108$samples
s
table(s$`treatment:ch1`)
table(sub("^\\S+\\s*", "", s$title))      # title suffixes: "", PHI, 24H, 48H, 72H

design_108 <- data.frame(
  sample    = s$sample,
  pair      = substr(s$title, 1, 12),
  condition = ifelse(s$`treatment:ch1` == "NT", "baseline",
              ifelse(s$`treatment:ch1` == "chronic", "treated", NA)),
  patient   = s$`patient id:ch1`
)
design_108 <- design_108[!is.na(design_108$condition), ]
table(design_108$condition)

# Overlap with the XENTURION discovery patients.
# Uses the result of section 4 if available; otherwise a quick run, since the
# paired models do not depend on the number of repeats or bootstrap resamples.
if (file.exists("ptpm_genome_wide.rds")) {
  genome_wide <- readRDS("ptpm_genome_wide.rds")
} else {
  example_data <- load_ptpm_example()
  genome_wide <- run_ptpm(
    expression = example_data$expression,
    metadata = example_data$metadata,
    top_k_values = 300, repeats = 2, n_boot = 200,
    gsva_sensitivity = FALSE, ssgsea_sensitivity = FALSE, robustness = FALSE,
    n_perm_coupling = 100, seed = 123,
    n_cores = parallel::detectCores() - 1
  )
}
disc_ids <- unique(substr(names(genome_wide$y), 1, 7))
length(disc_ids)                                                     # 72
length(intersect(unique(substr(design_108$patient, 1, 7)), disc_ids))

res_108 <- run_ptpm(expression = expr_filter(raw_108$expression), design = design_108,
                    top_k_values = c(3000, 300), k = 3, repeats = 20, n_boot = 2000,
                    seed = 123, n_cores = parallel::detectCores() - 1)
res_108$results_genome_wide
res_108$confidence_intervals
res_108$contrast_intervals
res_108$coupling_null
```

## Notes

- **Runtime:** reduce `repeats` first for a faster check. Each robustness and
  sensitivity block can be switched off (`gsva_sensitivity`, `ssgsea_sensitivity`,
  `robustness`) or limited to the first repeats (`gsva_sensitivity_repeats`,
  `ssgsea_repeats`, `robustness_repeats`).
- **Pairing permutation null:** the full null of the NT/CET pairing
  (`n_perm > 0`) is optional and computationally heavy. It is off by default.
- **Interaction pathways:** in `run_ptpm_hybrid_interactions()` the pathways used
  for interactions (`interaction_pathways`) are treated as fixed a priori. If
  they were chosen using the outcome on the same samples, the resulting
  estimates are optimistic.
- **Benchmark:** `run_ptpm_benchmark()` accepts only `run_ptpm()` output, which
  contains the baseline predictor matrix (`x_baseline`). The comparison with the
  elastic net is paired only when the benchmark is applied to that same output.
- **Reproducibility:** results are reproducible for a given `seed`, number of
  repeats and package versions (R, GSVA, glmnet; for the benchmark also xgboost,
  mixOmics and e1071). Report `sessionInfo()` together
  with the results.
