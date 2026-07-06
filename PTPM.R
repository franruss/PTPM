## PTPM_algoritm
## XENTURION is a population-level multidimensional resource of xenografts and tumoroids from metastatic colorectal cancer patients
## XENTURION is a bank of metastatic colorectal cancer (mCRC) Patient Derived Organoids (PDOs) obtained from Patient Derived Xenografts (PDXs).
## These models are annotated with different omics to advance our understanding of CRC.
# Clarification for columns “cell type”, “cell line” and “treatment”.
# The column “cell type” defines the type of samples
# LM: liver metastasis
# PR: primary (from colon)
# The column “cell line” defines the origin of samples
# H: patients
# X: xenografts
# O: organoids
# The column “treatment” describes the experiments done with cetuximab
# cetuxi: treated samples (72h (written as h) or 6w for xenografts)
# NT: non-treated control
# The word “BASALE” in the dataset indicates xenografts or organoids grown in normal conditions. They are NOT the non-treated control of the cetuximab experiments.

# È la tesi metodologica centrale, non un dettaglio tecnico. Il contributo di PTPM a una rivista di metodi non è "abbiamo predetto la risposta al cetuximab", 
# ma "abbiamo mostrato come stimare onestamente quella predicibilità in un disegno appaiato". La stima onesta è ciò che trasforma un'analisi di dati 
# esistenti in un contributo metodologico 

library(GEOquery)
library(R.utils)
library(affy)
library(data.table)
library(msigdbr)
library(GSVA)
library(glmnet)
library(caret)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pbapply)
library(caret)
set.seed(123)
here::here() 
setwd("C:/Users/Francesco Russo/OneDrive - CNR/machine learning")
# gse <- getGEO(filename = "GSE204805_series_matrix.txt.gz",GSEMatrix = TRUE)
# gse <- getGEO(filename = "C:/Users/Francesco Russo/Documents/GSE204805/GSE204805_merged_hs_mm.tsv.gz",GSEMatrix = TRUE)
#vsd_path <- "C:/Users/Francesco Russo/Documents/GSE204805/GSE204805_merged_hs_mm.tsv.gz"
vsd_path <- "GSE204805/GSE204805_vsd.tsv.gz"
expr_raw <- fread(vsd_path)
#gse <- getGEO("GSE204805", GSEMatrix = TRUE)
#expr_matrix <- exprs(gse[[1]])
expr_matrix <- as.data.frame(expr_raw)
dim(expr_matrix)
#[1] 34097  1226
expr_matrix[1:4,1:4]
#V1 CRC0018LMX0B02204TUMR14000 CRC0024LMH0000000000R10000 CRC0024LMX0A02201TUMR07000
#1     H_WASH7P                   7.393221                   7.540883                   7.433586
#2  H_MIR6859-1                   7.493525                   7.735542                   7.386376
#3 H_AL627309.7                   7.190104                   8.358940                   7.988427
rownames(expr_matrix) <- expr_matrix[, 1]# Remove first column
expr_matrix <- expr_matrix[, -1]# Convert to numeric matrix
expr_matrix <- as.matrix(expr_matrix)
mode(expr_matrix) <- "numeric"
dim(expr_matrix)
# 34097  1225
#saveRDS(expr_matrix, "GSE204805_expr_matrix.rds")
#expr_matrix <- readRDS("GSE204805_expr_matrix.rds")
expr_ids <- colnames(expr_matrix)
pheno <- data.frame(expr_id = expr_ids,stringsAsFactors = FALSE)
head(pheno)
dim(pheno)#[1] 1225    1
pheno$core_id <- substr(pheno$expr_id, 1, 7)
head(pheno$core_id)
#STEP 3 — Extract Model Type (LMX / LMO / LMH / PR etc.)
pheno$model_type <- substr(pheno$expr_id, 8, 10)
table(pheno$model_type)
#LMH LMO LMX PMX PRH PRX 
# 68 282 664   1  51 159 
####correggi errore di annotazione riscontrato
pheno$model_type[
  pheno$expr_id == "CRC1731PMX0A02002TUMR01000"
] <- "LMX"
table(pheno$model_type)
#LMH LMO LMX PRH PRX 
# 68 282 665  51 159 
# 1. Carica i metadati reali da GEO
gse <- getGEO("GSE204805", GSEMatrix = TRUE)
pheno_raw <- pData(gse[[1]])
# 2. Estrai sample_name (prima della virgola nel Title) e treatment
pheno_raw$sample_name <- trimws(sub(",.*", "", as.character(pheno_raw$title)))
parse_treatment <- function(title) {
  if (grepl("cetux_6w_NT|cetuxi_h_NT|cetux_h_NT", title)) return("NT")
  else if (grepl("cetux|cetuxi", title))                   return("CET")
  else if (grepl("BASALE", title))                         return("BASALE") #The word “basale” indicates xenografts or organoids grown in normal conditions. They are not the non-treated control of the cetuximab experiments
  else                                                      return(NA)
}
pheno_raw$treatment_geo <- sapply(as.character(pheno_raw$title), parse_treatment)

entire_codes <- substr(pheno_raw$title, 1, 26)
new_pheno = cbind(entire_codes ,pheno_raw[,c(1,2)])
head(new_pheno)

pheno_raw$entire_codes <- substr(pheno_raw$title, 1, 26)
pheno_raw$pair_id <- pheno_raw$entire_codes

pheno_raw$pair_id <- gsub("(LMO.*[0-9]{3})[0-9]R", "\\1XR", pheno_raw$pair_id)

# --- LMX: remove treatment-specific last digit block before T ---
pheno_raw$pair_id <- gsub("(LMX.*[0-9]{3})[0-9]{2}T", "\\1XXT", pheno_raw$pair_id)
pheno_raw$pair_id <- gsub("UMR[0-9]+$", "", pheno_raw$pair_id)
pair_table <- table(pheno_raw$pair_id, pheno_raw$treatment_geo)
# Recompute pairing
dim(pair_table) # 1077
#pair_table

pheno_raw$treatment_geo <- sapply(as.character(pheno_raw$title), parse_treatment)
pheno_raw$treatment_geo[1:10]
#[1] "BASALE" "BASALE" "BASALE" "BASALE" "BASALE" "BASALE" "BASALE" "BASALE" "BASALE" NA     
treatment_lookup <- setNames(pheno_raw$treatment_geo, pheno_raw$sample_name)
treatment_lookup[1:10]
#CRC0014LMX0B02202TUMR01000 CRC0017LMX0A02204TUMR02000 CRC0018LMO0B02004001R01000 CRC0018LMX0B02204TUMR14000 CRC0019LMX0B02202TUMR02000 CRC0021LMX0B02201TUMR02000 
#"BASALE"                   "BASALE"                   "BASALE"                   "BASALE"                   "BASALE"                   "BASALE" 
#CRC0022LMO0A04012002R01000 CRC0022LMO0A04013001R01000 CRC0022LMX0B02204TUMR03000 CRC0024LMH0000000000R10000 
#"BASALE"                   "BASALE"                   "BASALE"                         NA 
length(treatment_lookup) #1225
length(pheno$expr_id)#1225
pheno$expr_id[1:10]
#[1] "CRC0018LMX0B02204TUMR14000" "CRC0024LMH0000000000R10000" "CRC0024LMX0A02201TUMR07000" "CRC0024LMX0B02201TUMR01000" "CRC0028LMH0000000000R02000" "CRC0028LMX0B02204TUMR02000"
#[7] "CRC0029LMX0A02204TUMR05000" "CRC0050LMH0000000000R02000" "CRC0050LMX0A03210TUMR01000" "CRC0050PRX0A02201TUMR03000"
pheno$treatment_geo <- treatment_lookup[pheno$expr_id] # SUPERFLUO
#pheno$treatment_geo <- treatment_lookup
length(pheno$treatment_geo) # 1225
pheno$treatment_geo[1:10]
# [1] "BASALE" NA       "BASALE" "BASALE" NA       "BASALE" "BASALE" NA       "BASALE" "BASALE"
table(pheno$treatment_geo)
# BASALE    CET     NT 
# 859    126    121 
pheno_clean <- pheno[!is.na(pheno$treatment_geo), ]
dim(pheno_clean) # Keep only CET/NT/BASALE # 1106
pheno_na <- pheno[is.na(pheno$treatment_geo), ]
dim(pheno_na) # 119
#paired_ids <- rownames(table(pheno_raw$pair_id, pheno_raw$treatment_geo))[rowSums(table(pheno_raw$pair_id, pheno_raw$treatment_geo)[, c("CET","NT")] > 0) == 2]
paired_ids <- rownames(pair_table)[pair_table[, "CET"] > 0 & pair_table[, "NT"] > 0]
length(paired_ids)# 106
# Since we are using only the NT samples in the glmnet model
pheno_paired <- pheno_raw[pheno_raw$pair_id %in% paired_ids, ]
expr_paired  <- expr_matrix[, pheno_paired$entire_codes]
dim(expr_matrix)
dim(expr_paired )#[1] 34097 243
# Compute Response Score
# Then compute: EGFR Pathway Suppression Score 
# Use MSigDB: HALLMARK_EGF_SIGNALING, KRAS_SIGNALING_UP, MAPK_PATHWAY
# Compute GSVA before and after treatment. Then: EGFR_Suppression = Score_NT − Score_CET
# Large positive value = strong inhibition → Sensitive.
#  Define response based on downregulation of proliferation genes
# Modeling intrinsic transcriptional rewiring under EGFR blockade
# extend EDACA to:Paired drug perturbation modeling
pheno_paired$group <- paste(pheno_paired$pair_id, pheno_paired$treatment_geo,sep = "_")
group_levels <- unique(pheno_paired$group)
group_levels[1:12]
pheno_paired[1:12,1:2]
#table(pheno_paired$group)
#  average expression per group
expr_avg <- sapply(group_levels, function(g) {
  cols <- pheno_paired$entire_codes[pheno_paired$group == g]
  rowMeans(expr_matrix[, cols, drop = FALSE])
})
dim(expr_avg)#[1] 34097   243
# Step 1 — Compute GSVA on averaged matrix
hallmark <- msigdbr(species = "Homo sapiens",category = "H")
egfr_sets <- hallmark[hallmark$gs_name %in% c("HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_MTORC1_SIGNALING","HALLMARK_E2F_TARGETS"), ]
gene_sets <- split(egfr_sets$gene_symbol,egfr_sets$gs_name)
#gene_sets
#$HALLMARK_E2F_TARGETS
#[1] "AK2"      "ANP32E"   "ASF1A"    "ASF1B"    "ATAD2"    "AURKA"    "AURKB"   
# Now compute GSVA:
human_genes <- grepl("^H_", rownames(expr_avg))
expr_human <- expr_avg[human_genes, ]
rownames(expr_human) <- sub("^H_", "", rownames(expr_human))

# il leakage che introduce GSVA nella suppression matrix è minimo perché stai calcolando score su campioni paired (NT e CET dello stesso modello) 
# e la normalizzazione avviene a livello di coorte intera, non a livello di fold.
# Pathway activity scores were quantified using GSVA rather than ssGSEA. While ssGSEA generates enrichment scores independently for each sample, 
# GSVA estimates pathway activity relative to the expression distribution observed across the study cohort, providing a more stable representation 
# of coordinated pathway variation across tumor models. For PTPM, pathway suppression is defined as the paired difference between baseline (NT) 
# and post-treatment (CET) pathway activities. Because both measurements originate from the same tumor model and pathway definitions were 
# specified a priori using MSigDB Hallmark gene sets, the suppression phenotype reflects within-model transcriptional reprogramming rather 
# than information learned from drug-response labels. 
# We acknowledge that GSVA produces cohort-dependent pathway scores. However, in the present 
# study this property is advantageous because the predictive task aims to relate baseline pathway activity to differential therapeutic suppression 
# across models. In this setting, cohort-level normalization captures biologically meaningful between-model variation in pathway activity while 
# preserving the paired structure used to derive the suppression phenotype.

param <- gsvaParam(expr_human,gene_sets)
gsva_scores <- gsva(param) #una sola chiamata GSVA sulla matrice combinata a 212 campioni. 
# GSVA è stato calcolato congiuntamente su NT e CET — ottimo, è il punto che valida l'intero fenotipo di soppressione. 
# Se gsva() fosse stato chiamato separatamente su NT e su CET, la differenza NT−CET sarebbe un artefatto di normalizzazione.
gsva_scores[,1:10]
dim(gsva_scores) # 3 213

nt_cols  <- grepl("_NT$", colnames(expr_avg))
cet_cols <- grepl("_CET$", colnames(expr_avg))
gsva_nt  <- gsva_scores[, nt_cols]
gsva_cet <- gsva_scores[, cet_cols]
dim(gsva_nt)
#[1]   3 106
pair_ids <- intersect(sub("_NT$", "", colnames(gsva_nt)),sub("_CET$", "", colnames(gsva_cet)))

#The phenotype Suppression = GSVA(NT) − GSVA(CET) is constructed on pairs of the same tumor model.
#Becausethe pathways are defined a priori; NT and CET belong to the same patient/model; you are not using drug response to construct the scores;
#the difference between NT and CET represents a biological measure of transcriptional remodeling and not a supervised artifact.

suppression_matrix <- sapply(pair_ids, function(pid) {
  nt_col  <- paste0(pid, "_NT")
  cet_col <- paste0(pid, "_CET")
  gsva_scores[, nt_col] - gsva_scores[, cet_col] # IMPORTANT DEFINITION # Suppression = NT − CET
})

dim(suppression_matrix)  # 3 × 106
# Suppression = NT − CET. Because inhibition → pathway decreases.
# This gives: 3 pathways × 106 models. Now compute a single composite score. Stai modellando:"quanto collassa il programma EGFR" e non "il tumore risponde o no"
suppression_score <- colMeans(suppression_matrix) # IMPORTANT DEFINITION  which assumes w1 = w2 = w3 = 1/3 in the model: Suppression Score = w1·KRAS + w2·MTOR + w3·E2F
suppression_score[1:4]
#CRC0030LMX0B053XXT         CRC0031LMX0A050XXT         CRC0032LMX0A040XXT CRC0059LMO0A0300900XR01000 
#0.059941963               -0.002222302                0.032975742                0.053627182 
suppression_matrix[,1:4] # Positive values indicate pathway down-regulation following treatment
#HALLMARK_E2F_TARGETS               0.04669075        0.045915161         0.12794632                 0.70556695
#HALLMARK_KRAS_SIGNALING_UP        -0.01108370        0.018021896        -0.05769855                 0.01729233
#HALLMARK_MTORC1_SIGNALING          0.25135045        0.004106916         0.11729629                 0.34351016

summary(suppression_score) #a composite suppresion score per model was calculated as the mean suppression across the EGFR-axis pathway
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.37758   0.03267    0.17317   0.16148   0.28522   0.58088 

###########################################################################
# Adesso facciamo gli: istogrammi, densità, waterfall plot dei 106 modelli.
theme_ptpm <- theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),plot.title = element_text(face = "bold", hjust = 0.5)
  )
df <- data.frame(suppression_score = suppression_score)
# 1) Histogram
ggplot(df, aes(x = suppression_score)) +
  geom_histogram(bins = 20,color = "black",fill = "steelblue",alpha = 0.8) +
  geom_vline(xintercept = mean(df$suppression_score),linetype = 2,linewidth = 1) +
  labs(title = "Distribution of EGFR-axis Suppression Scores",x = "Suppression score",y = "Number of models") +
  theme_ptpm
# 2) Density
ggplot(df, aes(x = suppression_score)) +
  geom_density(fill = "steelblue",alpha = 0.4,linewidth = 1.2) +
  geom_vline(xintercept = mean(df$suppression_score),linetype = 2,linewidth = 1) +
  labs(title = "Density of EGFR-axis Suppression Scores",x = "Suppression score",y = "Density") +
  theme_ptpm
# 3) Waterfall plot. Ordina i modelli dal più resistente al più sensibile.
df_waterfall <- data.frame(model = names(suppression_score),score = suppression_score) %>%
  arrange(score) %>%
  mutate(model_rank = seq_len(n()) )
ggplot(df_waterfall,aes(x = reorder(model, score),y = score,fill = score > 0)) +
  geom_col(width = 0.9) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "firebrick","FALSE" = "steelblue"),guide = "none") +
  geom_hline(yintercept = 0,linewidth = 0.8) +
  labs(title = "Waterfall Plot of Cetuximab-Induced Pathway Suppression",x = "Models",y = "Suppression score") +
  theme_ptpm +
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

df_long <- as.data.frame(t(gsva_scores))
df_long$model_condition <- rownames(df_long)
df_long <- df_long %>%
  mutate(condition = ifelse(grepl("_NT$", model_condition),"NT","CET"),model = sub("_(NT|CET)$", "", model_condition)
  ) %>%
  pivot_longer(cols = starts_with("HALLMARK"),names_to = "pathway",values_to = "GSVA")
head(df_long)
df_long <- df_long %>%
  filter(grepl("_(NT|CET)$", model_condition)) %>%
  mutate(condition = sub("^.*_(NT|CET)$","\\1",model_condition),model = sub("_(NT|CET)$","",model_condition))
table(df_long$condition)
head(df_long)

ggplot(df_long,aes(x = condition,y = GSVA,fill = condition)) +
  geom_line(aes(group = model),alpha = 0.12,colour = "grey60") +
  geom_boxplot(width = 0.5,alpha = 0.8,outlier.shape = NA) +
  geom_jitter(width = 0.05,alpha = 0.25,size = 1) +
  stat_compare_means(paired = TRUE,method = "wilcox.test",label = "p.format") +
  facet_wrap( ~ pathway,nrow = 1) +
  scale_fill_manual(values = c(NT = "steelblue",CET = "firebrick")) +
  labs(title = "Cetuximab suppresses E2F, KRAS and mTORC1 signaling programs",x = NULL,y = "GSVA enrichment score") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",panel.grid = element_blank(),strip.text = element_text(face = "bold"))

########################################
# INIZIAMO con le prime definizioni
# Large positive value → CET strongly suppresses KRAS/mTOR/E2F → likely sensitive.
# Near zero → resistant.
# Negative value → Pathway paradoxically increases → strong resistance.
# Mechanistically defined cetuximab pathway suppression score Derived from paired perturbation RNA-seq
# Can baseline (NT) expression predict this suppression phenotype?
# STEP 1 — Extract Baseline Only
baseline_cols <- colnames(expr_avg)[grepl("_NT$", colnames(expr_avg))]
expr_baseline <- expr_avg[, baseline_cols]
colnames(expr_baseline) <- sub("_NT$", "", baseline_cols)
rownames(expr_baseline) <- sub("^H_", "", rownames(expr_baseline))
dim(expr_baseline) # 34097 × 106
# STEP 2 — Align with Response Label
x <- t(expr_baseline)
# STEP 3 — Reduce Dimensionality (Critical). We must avoid p >> n disaster.
# STEP 4 — Elastic Net (Robust Choice)

# So we model: Baseline transcriptome → Continuous suppression score.
# Align Continuous Outcome Correctly. Your alignment must be exact:

y_cont <- suppression_score[colnames(expr_baseline)] # defined as  colMeans(suppression_matrix) # IMPORTANT DEFINITION
# y_cont is derived from paired NT–CET data. Even though you're only using NT in x y_cont contains CET information
# If two samples from the same biological pair context (or related replicates) leak across folds. model indirectly "sees" test signal

sum(is.na(y_cont))  # must be 0
# Elastic Net Regression 
# Right now cross-validation is random at the sample level: trainControl(method = "repeatedcv", number = 3, repeats = 20)
# That means: Samples from the same patient (e.g., CRC0059) can appear both in training AND test folds
# But your data structure is: Same patient → multiple models (LMO, LMX, replicates) Same patient → highly correlated transcriptomes
# So your model is partially “seeing the test data” during training. This inflates performance. 
# Questa inflazione deriva da due forme di information leakage: situazioni in cui informazione che dovrebbe restare confinata 
# al test set influenza, anche indirettamente, la costruzione del modello. Due meccanismi sono particolarmente insidiosi nei disegni come quello qui adottato. 
# 1) Il primo è il feature-selection leakage: quando i geni più variabili vengono selezionati usando l'intero dataset prima della cross-validation, i fold 
# di test contribuiscono già alla scelta delle feature, e la performance stimata risulta ottimistica. Questo si può correggere nel codice.
# 2) Il secondo è il mathematical coupling: quando un fenotipo definito come differenza (baseline meno trattato) viene predetto a partire 
# dalle sue stesse componenti baseline, una porzione dell'associazione è garantita per costruzione algebrica, indipendentemente da qualsiasi segnale biologico.


pheno$core_id <- substr(pheno$expr_id, 1, 7)

groups <- substr(colnames(expr_baseline), 1, 7) # group the same patients together
x <- t(expr_baseline)
vars <- apply(x, 2, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:3000] # Reduce Dimensionality (Critical). We must avoid p >> n disaster.
x_reduced <- x[, top_genes] 
length(groups) == nrow(x_reduced)
# All samples from the same patient go into the same fold. So the model is tested on completely unseen patients. That's why we create group_folds.
# Quindi la mia definizione di analisi naïve è:
# Grouped patient-level CV (via groupKFold, corretto rispetto alla struttura raggruppata)
# Feature selection pre-CV (top 3000 geni selezionati una volta sul dataset intero)
# Nessuna coupling correction
# Nessuna correzione per FS leakage
# Corregge il grouping, ma non la FS leakage e non il coupling. Quindi Presento la naïve come "grouped ma senza altre correzioni" 
group_folds <- function(groups, k, repeats) {
  folds <- list()
  for (r in 1:repeats) {
    fold_ids <- groupKFold(groups, k = k)
    names(fold_ids) <- paste0("Fold", seq_along(fold_ids), "_Rep", r)
    folds <- c(folds, fold_ids)
  }
  return(folds)
}
folds <- group_folds(groups, k = 3, repeats = 20)
ctrl <- trainControl(method = "cv",index = folds,savePredictions = "final")
#Controllo allineamento y_cont: fai sempre:
stopifnot(all(names(y_cont) == rownames(x_reduced)))
fit_reg <- caret::train(x = x_reduced, y = y_cont,method = "glmnet",trControl = ctrl,tuneLength = 20, tuneGrid = expand.grid(alpha = seq(0, 0.1, length = 11),lambda = seq(0.001, 0.1, length = 20)))
#alpha spans from 0 to 0.1 since in the elastic net regression I prefer to keep all variables since the signal is dense and distributed across numerous genes
sd(y_cont)#  0.1802557
min(fit_reg$results$RMSE) # 0.141784
min(fit_reg$results$RMSE) < sd(y_cont) # If RMSE < SD → model captures signal.
subset(fit_reg$results,fit_reg$results$RMSE == min(fit_reg$results$RMSE))
#    alpha lambda      RMSE  Rsquared       MAE     RMSESD RsquaredSD      MAESD
#220   0.1    0.1 0.1323092 0.5067064 0.1019908 0.03849312  0.1269112 0.01836154
#controllo se il segnale deriva dal dal coupling
# Nel modello elastic net x_reduced contiene 3000 geni baseline, ma non contiene E2F_NT, KRAS_NT, mTOR_NT come covariate esplicite. Quindi il problema non è forte così come sembra. Tuttavia il coupling c'è
baseline_score <- colMeans(gsva_nt)
cor(baseline_score, suppression_score) #0.62

m0 <- lm(suppression_score ~ baseline_score) #Questa procedura rimuove l'effetto meccanico dovuto al fatto che NT compare in entrambi i lati dell'equazione.
resid_response <- resid(m0) #rappresenta la componente della risposta a cetuximab non spiegabile dal livello basale dei pathway EGFR.
#"Quali caratteristiche del trascrittoma baseline predicono una risposta a cetuximab maggiore o minore di quella attesa dato il livello basale dei pathway EGFR?"
stopifnot(all(names(resid_response ) == rownames(x_reduced)))
fit_reg_resid <- caret::train(x = x_reduced, y = resid_response ,method = "glmnet",trControl = ctrl,tuneLength = 20, tuneGrid = expand.grid(alpha = seq(0, 0.1, length = 11),lambda = seq(0.001, 0.1, length = 20)))
subset(fit_reg_resid$results, fit_reg_resid$results$RMSE == min(fit_reg_resid$results$RMSE))
#    alpha lambda      RMSE  Rsquared        MAE     RMSESD RsquaredSD      MAESD
#220   0.1    0.1 0.1235081 0.2727583 0.09815819 0.01694652   0.173505 0.01274438
cor(baseline_score, resid_response) #[1] -5.646499e-17

####################################################################################
# alpha = 1 → LASSO (sparse, feature selection)  When you suspect only a small subset of features are actually important.
# alpha = 0 → Ridge (dense, shrinkage) When many features are relevant or highly correlated.
# alpha ≈ 0.1 → mostly ridge
# Keep only best tuning parameters
# many correlated predictors (RNA-seq → huge collinearity)
# LASSO becomes unstable
# Ridge is much more stable
# Elastic net shifts toward low alpha
# So the model is basically telling you: “Don’t try to pick a few genes. The signal is distributed.”
pred <- fit_reg$pred  
best <- fit_reg$bestTune
pred_best <- pred[pred$alpha == best$alpha & pred$lambda == best$lambda, ]
# instead of computing R² per fold it’s better to compute it globally, not fold-by-fold.
# Global R² (across ALL held-out predictions)
R2_global <- cor(pred_best$obs, pred_best$pred)^2
R2_global # 0.444
coefs <- coef(fit_reg$finalModel, s = fit_reg$bestTune$lambda)
head(coefs)
coefs_df <- as.data.frame(as.matrix(coefs))
head(coefs_df)
coefs_df$Feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "Coefficient"
coefs_df <- coefs_df[coefs_df$Coefficient != 0, ]
top30 <- coefs_df[order(abs(coefs_df$Coefficient), decreasing = TRUE), ][1:30, ]
top30 <- top30[top30$Feature != "(Intercept)", ]
ggplot(top30,
       aes(x = reorder(Feature, Coefficient),
           y = Coefficient,
           fill = Coefficient > 0)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "",
    y = "Beta Coefficient"
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "firebrick",
      "FALSE" = "steelblue"
    ),
    guide = "none"
  ) +
  ggtitle("Top 30 out of the 3000 most variable genes") +
  theme_minimal(base_size = 14)
# In the plot I have
# MYBL2 → cell cycle regulator
# MCM4 → DNA replication
# PRKDC → DNA damage repair
# IGF1R → growth signaling
# If this gene increases (by 1 standardized unit), how does y change?” Mathematically:
# y = ∑j βj + xj + ...
# So:Positive coefficient → gene is associated with higher outcome Negative coefficient → gene is associated with lower outcome
# Example H_MYBL2 = +0.0036 → higher MYBL2 → higher predicted y
# H_AL365181.2 = -0.0050 → higher expression → lower predicted y
# These are not independent importance scores Because Genes are correlated
# Example: MYBL2, MCM4, PRKDC → all cell cycle / proliferation
# YBL2, MCM4, ATAD2 → cell cycle / proliferation
# PRKDC, GADD45G → DNA damage response
# IGF1R → growth signaling
# ENO1 → glycolysis
# Your model is capturing a proliferation / growth axis
# The predictive signal is weak at the level of individual genes (small coefficients), but consistent patterns emerge 
# at the level of biological programs, particularly proliferation, DNA repair, and growth signaling pathways.
pred <- predict(fit_reg,   x_reduced)
round(cor(y_cont, pred),2) #  0.95


# The baseline pathway score was moderately correlated with the suppression score (r = 0.62), indicating the presence of partial mathematical coupling between baseline pathway activity 
# and treatment-induced changes. To assess the robustness of the findings, we regressed the suppression score against the baseline pathway score and repeated all predictive 
# analyses using the residualized outcome. As expected, predictive performance decreased; however, the elastic-net model retained substantial predictive ability 
# (cross-validated R² = 0.27 versus 0.51 in the original analysis). These results suggest that although part of the predictive signal is attributable to baseline-to-change 
# coupling, a substantial fraction reflects independent transcriptional information present before treatment.
# e se volessi usare: CET ~ NT + predictors ? Dal punto di vista statistico questa è la soluzione più elegante.
# Invece di modellare NT − CET modelli direttamente: CET = β0 + β1*NT + β2*X1 + ... Questo è il classico approccio ANCOVA.
# Però qui nasce un problema interpretativo diverso. Con NT - CET un valore alto significa forte soppressione da cetuximab. cioè una misura diretta di sensibilità.
# Con CET sto invece predicendo l'attività residua dopo trattamento. Quindi: CET basso = sensibile, CET alto = resistente. L'interpretazione è meno immediata.
# Quindi posso modellare outcome = suppression score (biologicamente intuitivo) e in più fare un'analisi di robustezza:
# 1) outcome = CET score
# 2) aggiustamento per NT score
# e scriverei: To address the potential mathematical coupling between baseline pathway activity and treatment-induced changes, 
# we performed an additional ANCOVA analysis in which post-treatment pathway activity was modeled while adjusting for baseline pathway activity. 
# The results were consistent with those obtained using the suppression score, indicating that the observed associations were not solely driven by baseline-to-change coupling.
# Questa scelta conserva l'interpretabilità biologica del suppression score ma mostra anche che le conclusioni non dipendono dalla scelta dell'outcome come differenza NT−CET.
# Perché dimostra contemporaneamente: esiste una componente di coupling rimane comunque una sostanziale capacità predittiva indipendente.
# I already reduced the number of genes to 300 genes. And I got: R2 = 0.35. Elastic Net residualized R2 = 0.22. So I kept 3000 genes.
# Compared with the earlier 3000-gene model, the performance did not improve dramatically. Reducing to 300 genes reduced overfitting, but did not uncover a much stronger signal.
coef_orig  <- coef(fit_reg$finalModel, s = fit_reg$bestTune$lambda)
coef_resid <- coef(fit_reg_resid$finalModel,s = fit_reg_resid$bestTune$lambda)
# Ho confrontato l'intersezione tra i geni selezionati nei due modelli (the model with inflation and residualized one). 
# Se KRAS-related, cell-cycle o EGFR-related features sopravvivono alla residualizzazione, il messaggio biologico diventa più robusto.
co_orig <- as.matrix(coef_orig)
co_resid <- as.matrix(coef_resid)
sel_orig <- rownames(co_orig)[co_orig[,1] != 0]
sel_resid <- rownames(co_resid)[co_resid[,1] != 0]
length(sel_orig)
length(sel_resid)
common <- intersect(sel_orig, sel_resid)
length(common)
common
#[1] "(Intercept)" "PLA2G2A"     "NT5E"        "AC108865.1"  "PRNP"        "STOM"        "STC2"        "HOXB6"       "ETV4"        "FAM107B"     "WASF3"       "CRIP2"       "MEF2C"      
#[14] "TAC1"        "GUCA2A"      "C16orf62"    "AC108865.2"  "CGREF1"      "AFAP1-AS1"   "CAPN12"      "IER3"        "CDC42EP3"    "MLF1"        "NUDT12"      "NRSN2"       "VAT1"       
#[27] "PTGDR"       "CD109"       "TCEAL8"      "MYBL2"       "HLA-E"       "METTL7B"     "KREMEN1"     "SUSD3"       "NME4"        "GATA2"       "TLN2"        "MCM4"        "ADTRP"      
#[40] "GNAI2"       "SNHG7"       "TSEN2"       "KLF2"        "BATF"        "STEAP1"      "GPR87"       "ST5"         "SPATA20"     "PRKDC"       "ELL2"        "ABLIM3"      "TMEM25"     
#[53] "IFNGR2"      "GADD45G"     "DEGS1"       "PLA2G2F"     "CACNB2"      "ACTG1"       "SERPINF1"    "ZNF880"      "MANEAL"      "AL365181.2"  "MFNG"     
# Questo risultato è più interessante del semplice valore di R². contiene un' informazione aggiuntiva sulla risposta a cetuximab. Alcuni geni che saltano subito all'occhio sono:
#  | Gene            | Interesse biologico                                                    |
#  | --------------- | ---------------------------------------------------------------------- |
#  | **ETV4**        | target della cascata EGFR/MAPK, spesso associato a tumori EGFR-driven  |
#  | **MYBL2**       | regolatore del ciclo cellulare, funzionalmente vicino ai programmi E2F |
#  | **MCM4**        | replicazione del DNA, tipico target E2F                                |
#  | **IER3**        | gene immediato di risposta a MAPK/EGFR                                 |
#  | **KREMEN1**     | pathway WNT                                                            |
#  | **GATA2**       | regolatore trascrizionale importante in plasticità cellulare           |
#  | **PRKDC**       | DNA damage response                                                    |
#  | **GADD45G**     | stress response e arresto proliferativo                                |
#  | **NT5E (CD73)** | immunomodulazione e resistenza terapeutica                             |
# Quindi possiamo dire che esiste un mathematical coupling reale: cor = 0.62 tra baseline score e suppression score. Dopo aver rimosso completamente tale componente cor ≈ 0 e R2 passa da 0.51 a 0.27. 
# Il modello continua a predire significativamente la risposta. Un nucleo stabile di geni (tra cui ETV4, MYBL2, MCM4, IER3, ecc.) viene selezionato in entrambe le analisi.
# A questo punto non sostituirei il modello originale; lo terrei come analisi principale e aggiungerei il modello residualizzato come sensitivity analysis, 
# riportando che circa metà della capacità predittiva persiste dopo la rimozione della componente dovuta al baseline pathway score.
co_orig <- as.matrix(coef_orig)
co_resid <- as.matrix(coef_resid)
co_orig <- co_orig[common, , drop=FALSE]
co_resid <- co_resid[common, , drop=FALSE]
cor(abs(co_orig[,1]),abs(co_resid[,1])) #Se questa correlazione è alta (>0.5), allora i geni importanti rimangono importanti anche dopo aver eliminato il coupling.
# 0.91
#Se il segnale fosse stato quasi interamente dovuto al mathematical coupling, 
#ti aspetteresti: drastico cambiamento dei geni selezionati; coefficienti completamente diversi; correlazione vicina a 0.
#Invece osservi: molti geni comuni; ranking dei coefficienti molto simile; correlazione ≈ 0.91.
#Questo suggerisce che il modello residualizzato sta catturando sostanzialmente la stessa struttura biologica del modello originale.
#Questo è il tipo di controllo che distingue un lavoro esplorativo da un framework metodologico ben pensato.
#"Baseline-to-change coupling contributes to the signal, but does not fully account for the predictive performance or feature selection."
#The consistency of feature selection across the original and residualized models — evidenced by a coefficient correlation of 0.91 and the preservation of biologically coherent EGFR-axis 
#regulators including ETV4, MYBL2, MCM4, IER3, and PRKDC — indicates that the predictive signal cannot be attributed to baseline-to-change mathematical coupling alone, but reflects 
#genuine transcriptional encoding of cetuximab pathway dependency.
#The R² tells you how well the model predicts, but it doesn't tell you WHY it predicts.
#If coupling were the only driver of the signal, you would expect the model to select different genes after residualization, essentially noise. 
#Instead, the coefficient correlation of 0.91 between the two models and a set of common genes that tells a  biological story.
#ETV4 is a transcription factor of the PEA3 family, directly induced by the RAS/MAPK pathway downstream of EGFR. Its presence in the 
#residualized model means that its basal expression predicts a component of the cetuximab response that cannot be explained by the basal level of GSVA pathways—thus capturing 
#something more subtle than just aggregated pathway activity.
#MYBL2 and MCM4 together tell the same story: they are both direct transcriptional targets of E2F, linked to G1/S progression. Their co-selection suggests that the model is 
#capturing cell cycle status at a more granular level than the aggregated GSVA HALLMARK_E2F_TARGETS score.
#IER3 is an immediate-response gene to MAPK and NF-κB, often upregulated in the early stages of cellular stress responses, including 
#those induced by EGFR inhibitors. Its basal expression may reflect pre-existing response competence.
#PRKDC encodes DNA-PK, involved in double-strand repair. Its presence in the residualized model is biologically plausible because tumors with high replication stress—often 
#associated with high E2F/MYBL2/MCM4 activity—tend to be more dependent on EGFR for survival and therefore respond better to its inhibition.
#NT5E (CD73) is a marker of immunosuppression and therapeutic resistance, and suggests that the immunological 
#microenvironment—even in a PDX/PDXT system purified from human immune cells—leaves a transcriptional signature in tumor cells that correlates with resistance.

###############################
#Adesso riduciamo a 300 geni
top_genes <- names(sort(vars, decreasing = TRUE))[1:300] # Reduce Dimensionality (Critical). We must avoid p >> n disaster.
x_reduced <- x[, top_genes] 
length(groups) == nrow(x_reduced)
folds <- group_folds(groups, k = 3, repeats = 20)
ctrl <- trainControl(method = "cv",index = folds,savePredictions = "final")
#Controllo allineamento y_cont: fai sempre:
stopifnot(all(names(y_cont) == rownames(x_reduced)))
fit_reg <- caret::train(x = x_reduced, y = y_cont,method = "glmnet",trControl = ctrl,tuneLength = 20, tuneGrid = expand.grid(alpha = seq(0, 0.1, length = 11),lambda = seq(0.001, 0.1, length = 20)))
sd(y_cont)#  
min(fit_reg$results$RMSE) # 0.
min(fit_reg$results$RMSE) < sd(y_cont) # If RMSE < SD → model captures signal.
subset(fit_reg$results,fit_reg$results$RMSE == min(fit_reg$results$RMSE))
#   alpha      lambda      RMSE  Rsquared       MAE     RMSESD RsquaredSD      MAESD
#1      0 0.001000000 0.1506061 0.3237159 0.1217498 0.01694048  0.1189354 0.01347634

m0 <- lm(suppression_score ~ baseline_score) #Questa procedura rimuove l'effetto meccanico dovuto al fatto che NT compare in entrambi i lati dell'equazione.
resid_response <- resid(m0) #rappresenta la componente della risposta a cetuximab non spiegabile dal livello basale dei pathway EGFR.
#"Quali caratteristiche del trascrittoma baseline predicono una risposta a cetuximab maggiore o minore di quella attesa dato il livello basale dei pathway EGFR?"
stopifnot(all(names(resid_response ) == rownames(x_reduced)))
fit_reg_resid <- caret::train(x = x_reduced, y = resid_response ,method = "glmnet",trControl = ctrl,tuneLength = 20, tuneGrid = expand.grid(alpha = seq(0, 0.1, length = 11),lambda = seq(0.001, 0.1, length = 20)))
subset(fit_reg_resid$results, fit_reg_resid$results$RMSE == min(fit_reg_resid$results$RMSE))
#     alpha   lambda      RMSE  Rsquared        MAE     RMSESD RsquaredSD      MAESD
#1      0 0.001000000 0.1281993 0.2172025 0.1032625 0.01254216  0.1468001 0.01168498
cor(baseline_score, resid_response) #[1] -5.646499e-17

set.seed(123)
r2_distribution <- pbsapply(1:10, function(i) {
  folds_i <- group_folds(groups, k = 3, repeats = 5)
  ctrl_i <- trainControl(method = "cv", index = folds_i, savePredictions = "final")
  fit_i <- caret::train(x = x_reduced, y = y_cont,
                        method = "glmnet", trControl = ctrl_i,
                        tuneGrid = expand.grid(alpha = seq(0, 0.1, length = 11),lambda = seq(0.001, 0.1, length = 20)))
  best_i <- fit_i$bestTune
  pred_i <- fit_i$pred
  pred_best_i <- pred_i[pred_i$alpha == best_i$alpha & pred_i$lambda == best_i$lambda, ]
  cor(pred_best_i$obs, pred_best_i$pred)^2
})
# median CV R² = X (95% CI: Y–Z across 10 repeated cross-validation runs)
median(r2_distribution)  # 0.29
quantile(r2_distribution, c(0.025, 0.975))
#  2.5%     97.5% 
#  0.47      0.53 

set.seed(123)

r2_distribution_resid <- pbsapply(1:10, function(i) {
  folds_i <- group_folds(groups, k = 3, repeats = 5)
  ctrl_i <- trainControl(method = "cv", index = folds_i, savePredictions = "final")
  fit_i <- caret::train(x = x_reduced, y = resid_response, # resid_response
                        method = "glmnet", trControl = ctrl_i,
                        tuneLength = 20, tuneGrid = expand.grid(alpha = seq(0, 0.1, length = 11),lambda = seq(0.001, 0.1, length = 20)))
  best_i <- fit_i$bestTune
  pred_i <- fit_i$pred
  pred_best_i <- pred_i[pred_i$alpha == best_i$alpha & pred_i$lambda == best_i$lambda, ]
  cor(pred_best_i$obs, pred_best_i$pred)^2
})
# median CV R2 resid = X (95% CI: Y–Z across 10 repeated cross-validation runs)
median(r2_distribution_resid)  # 0.26
quantile(r2_distribution_resid, c(0.025, 0.975))
#  2.5%     97.5% 
#  0.23      0.29 

############# LEAKING REMOVAL
# Per eliminare completamente il leakage dovuto alla selezione dei geni devo fare la feature selection dentro ogni fold. 
# Devo spostare la selezione dei 3000 geni all'interno del ciclo CV.
##############################################################################
# fit_reg (genome-wide) — LEAKAGE-FREE + (opzionale) COUPLING-CORRECTED
#
# Predittori: top-K geni per varianza (K = 3000 o 300), selezionati DENTRO
# ogni training fold.
#
# Fix 1: groupKFold ritorna TRAINING indices
# Fix 2: alpha tuning via inner CV (non un vettore a cv.glmnet)
# Fix 3: coupling correction opzionale (m0 stimato sul training, applicato al test)
##############################################################################

# ---------- Funzione riusabile ----------
evaluate_genome_wide_leakfree <- function(x_full, y, folds, top_k,
                                          baseline_score  = NULL,
                                          alpha_grid      = seq(0, 0.1, length = 11),
                                          coupling_correct = FALSE) {
  
  # Sanity check
  if (coupling_correct) {
    if (is.null(baseline_score)) {
      stop("coupling_correct = TRUE richiede baseline_score in input.")
    }
    stopifnot(length(baseline_score) == nrow(x_full))
    stopifnot(!is.null(names(baseline_score)))
    stopifnot(identical(names(baseline_score), rownames(x_full)))
  }
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    # FIX: folds[[i]] sono TRAINING indices, non test
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(x_full)), train_idx)
    
    x_train <- x_full[train_idx, ]
    x_test  <- x_full[test_idx, ]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # ----- Feature selection SOLO sul training: top-K per varianza -----
    vars      <- apply(x_train, 2, var)
    top_genes <- names(sort(vars, decreasing = TRUE))[1:top_k]
    x_train   <- x_train[, top_genes, drop = FALSE]
    x_test    <- x_test[,  top_genes, drop = FALSE]
    
    # ----- Coupling correction opzionale (m0 stimato SOLO sul training) -----
    if (coupling_correct) {
      b_train <- baseline_score[train_idx]
      b_test  <- baseline_score[test_idx]
      
      m0 <- lm(y_train ~ b_train)
      y_train_corr <- resid(m0)
      
      # Applica la stessa correzione al test: y_test - predizione da baseline
      pred_from_baseline <- coef(m0)[1] + coef(m0)[2] * b_test
      y_test_corr <- y_test - pred_from_baseline
      
      y_train <- y_train_corr
      y_test  <- y_test_corr
    }
    
    # ----- Tuning di alpha via inner CV (5-fold) -----
    set.seed(123 + i)
    inner_foldid <- sample(rep(1:5, length.out = length(y_train)))
    
    cv_results <- lapply(alpha_grid, function(a) {
      cvfit <- cv.glmnet(x_train, y_train, alpha = a, foldid = inner_foldid)
      list(alpha = a, cvm_min = min(cvfit$cvm), cvfit = cvfit)
    })
    best_idx <- which.min(sapply(cv_results, function(r) r$cvm_min))
    best_fit <- cv_results[[best_idx]]$cvfit
    
    pred <- predict(best_fit, newx = x_test, s = "lambda.min")
    all_pred <- c(all_pred, as.numeric(pred))
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Preparazione dati ----------

x_baseline_matrix <- t(expr_baseline)   # 106 × 22741 geni human

# Pulizia rownames
rownames(x_baseline_matrix) <- sub("_NT$", "", rownames(x_baseline_matrix))

# baseline_score allineato per nome
baseline_score_clean <- baseline_score
names(baseline_score_clean) <- sub("_NT$", "", names(baseline_score_clean))
baseline_score_aligned <- baseline_score_clean[rownames(x_baseline_matrix)]

# y_cont allineato
y_cont_aligned <- y_cont[rownames(x_baseline_matrix)]

# Verifiche
stopifnot(identical(names(baseline_score_aligned), rownames(x_baseline_matrix)))
stopifnot(identical(names(y_cont_aligned), rownames(x_baseline_matrix)))
stopifnot(sum(is.na(baseline_score_aligned)) == 0)


# ---------- Run: 4 combinazioni (K ∈ {3000, 300} × {leakage-free, fully corrected}) ----------

cat("\n--- fit_reg (K=3000), LEAKAGE-FREE ONLY ---\n")
R2_3000_leakfree <- evaluate_genome_wide_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 3000,
  coupling_correct = FALSE
)
cat("  R² =", round(R2_3000_leakfree, 4), "\n")

cat("\n--- fit_reg (K=3000), FULLY CORRECTED ---\n")
R2_3000_fully_corrected <- evaluate_genome_wide_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 3000,
  baseline_score   = baseline_score_aligned,
  coupling_correct = TRUE
)
cat("  R² =", round(R2_3000_fully_corrected, 4), "\n")

cat("\n--- fit_reg_300 (K=300), LEAKAGE-FREE ONLY ---\n")
R2_300_leakfree <- evaluate_genome_wide_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 300,
  coupling_correct = FALSE
)
cat("  R² =", round(R2_300_leakfree, 4), "\n")

cat("\n--- fit_reg_300 (K=300), FULLY CORRECTED ---\n")
R2_300_fully_corrected <- evaluate_genome_wide_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 300,
  baseline_score   = baseline_score_aligned,
  coupling_correct = TRUE
)
cat("  R² =", round(R2_300_fully_corrected, 4), "\n")


# ---------- Riepilogo -----------
results_genome_wide <- data.frame(
  Model             = c("fit_reg (K=3000)", "fit_reg_300 (K=300)"),
  Leakage_free_only = c(R2_3000_leakfree, R2_300_leakfree),
  Fully_corrected   = c(R2_3000_fully_corrected, R2_300_fully_corrected),
  Delta_coupling    = c(R2_3000_leakfree - R2_3000_fully_corrected,
                        R2_300_leakfree - R2_300_fully_corrected)
)
cat("\n=== Genome-wide models (top-K variable genes) ===\n")
print(results_genome_wide, row.names = FALSE)
#              Model Leakage_free_only Fully_corrected Delta_coupling
#   fit_reg (K=3000)         0.4413085       0.2340480      0.2072604
#fit_reg_300 (K=300)         0.2580619       0.1295018      0.1285601

##################################
##### BENCHMARK against other methods
###############################

#### Proviamo a fare un confronto tra Elastic Net vs XGBoost
## (nested) Per ogni fold esterno fai una CV interna per scegliere il numero di alberi.
# In ML moderno molti userebbero: XGBoost invece di Elasic Net. Ma io ho privilegiato l'interpretabilità.
##############################################################################
# XGBoost — NAÏVE (grouped patient-level CV + FS pre-CV, no coupling correction)
##############################################################################

library(xgboost)

# ---------- Helper: estrazione robusta di best_iteration ----------
get_best_iter <- function(cv_fit) {
  if (!is.null(cv_fit$best_iteration)) return(cv_fit$best_iteration)
  eval_log <- cv_fit$evaluation_log
  test_col <- grep("^test_.*_mean$", names(eval_log), value = TRUE)[1]
  eval_log[which.min(get(test_col)), iter]
}

# ---------- Funzione naïve ----------
evaluate_xgboost_naive <- function(x_full, y, folds, top_k,
                                   xgb_params  = NULL,
                                   max_nrounds = 1000) {
  if (is.null(xgb_params)) {
    xgb_params <- list(
      objective        = "reg:squarederror",
      eval_metric      = "rmse",
      eta              = 0.05,
      max_depth        = 3,
      subsample        = 0.8,
      colsample_bytree = 0.8,
      nthread          = 4
    )
  }
  # NAÏVE: feature selection UNA VOLTA sull'intero dataset (leakage per costruzione)
  vars       <- apply(x_full, 2, var)
  top_genes  <- names(sort(vars, decreasing = TRUE))[1:top_k]
  x_selected <- x_full[, top_genes, drop = FALSE]
  all_pred <- c()
  all_obs  <- c()
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]                                      # training indices
    test_idx  <- setdiff(seq_len(nrow(x_selected)), train_idx)
    x_train <- x_selected[train_idx, ]
    x_test  <- x_selected[test_idx, ]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    dtrain <- xgb.DMatrix(as.matrix(x_train), label = y_train)
    dtest  <- xgb.DMatrix(as.matrix(x_test),  label = y_test)
    set.seed(123 + i)
    cv_fit <- xgb.cv(
      params                = xgb_params,
      data                  = dtrain,
      nrounds               = max_nrounds,
      nfold                 = 5,
      early_stopping_rounds = 20,
      verbose               = 0
    )
    best_nrounds <- get_best_iter(cv_fit)
    fit <- xgb.train(
      params  = xgb_params,
      data    = dtrain,
      nrounds = best_nrounds,
      verbose = 0
    )
    pred <- predict(fit, dtest)
    all_pred <- c(all_pred, pred)
    all_obs  <- c(all_obs,  y_test)
  }
  cor(all_pred, all_obs)^2
}
# ---------- Preparazione ----------
x_baseline_matrix <- t(expr_baseline)
rownames(x_baseline_matrix) <- sub("_NT$", "", rownames(x_baseline_matrix))
y_cont_aligned <- y_cont[rownames(x_baseline_matrix)]
stopifnot(identical(names(y_cont_aligned), rownames(x_baseline_matrix)))
# ---------- Run ----------
cat("\n--- XGBoost (K=3000), NAÏVE ---\n")
R2_xgb_3000_naive <- evaluate_xgboost_naive(x_full = x_baseline_matrix, y = y_cont_aligned, folds = folds, top_k = 3000)
cat("  R² =", round(R2_xgb_3000_naive, 4), "\n")
#R² = 0.29 
cat("\n--- XGBoost (K=300), NAÏVE ---\n")
R2_xgb_300_naive <- evaluate_xgboost_naive(x_full = x_baseline_matrix, y = y_cont_aligned, folds = folds, top_k = 300)
cat("  R² =", round(R2_xgb_300_naive, 4), "\n")
#  R² = 0.20
# XGBoost	fa peggio. Perché l'Elastic Net è particolarmente più adatto ai problemi genomici ad alta dimensionalità e basso numero di campioni.

##############################################################################
# XGBoost — LEAKAGE-FREE + (opzionale) COUPLING-CORRECTED
#
# Predittori: top-K geni per varianza, selezionati DENTRO ogni training fold
#
# Fix 1: groupKFold ritorna TRAINING indices (già presente nel tuo codice)
# Fix 2: eval_metric = "rmse" invece di "auc" (regressione, non classificazione)
# Fix 3: best_iteration estratto in modo robusto (workaround per xgboost 2.x
#        dove $best_iteration può ritornare NULL nonostante early stopping)
# Fix 4: coupling correction opzionale (m0 stimato sul training, applicata al test)
##############################################################################
library(xgboost)
# ---------- Helper: estrazione robusta di best_iteration ----------
# In xgboost 2.x, cv_fit$best_iteration può essere NULL anche quando l'early
# stopping ha funzionato. Estrazione difensiva dall'evaluation_log.
get_best_iter <- function(cv_fit) {
  if (!is.null(cv_fit$best_iteration)) {
    return(cv_fit$best_iteration)
  }
  eval_log <- cv_fit$evaluation_log
  test_col <- grep("^test_.*_mean$", names(eval_log), value = TRUE)[1]
  if (is.na(test_col)) stop("Nessuna colonna test_*_mean trovata in evaluation_log")
  eval_log[which.min(get(test_col)), iter]
}


# ---------- Funzione riusabile ----------
evaluate_xgboost_leakfree <- function(x_full, y, folds, top_k,
                                      baseline_score  = NULL,
                                      xgb_params      = NULL,
                                      max_nrounds     = 1000,
                                      coupling_correct = FALSE) {
  
  # Default hyperparameters
  if (is.null(xgb_params)) {
    xgb_params <- list(
      objective        = "reg:squarederror",
      eval_metric      = "rmse",
      eta              = 0.05,
      max_depth        = 3,
      subsample        = 0.8,
      colsample_bytree = 0.8,
      nthread          = 4
    )
  }
  
  # Sanity check
  if (coupling_correct) {
    if (is.null(baseline_score)) {
      stop("coupling_correct = TRUE richiede baseline_score in input.")
    }
    stopifnot(length(baseline_score) == nrow(x_full))
    stopifnot(!is.null(names(baseline_score)))
    stopifnot(identical(names(baseline_score), rownames(x_full)))
  }
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    # FIX: folds[[i]] sono TRAINING indices
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(x_full)), train_idx)
    
    x_train <- x_full[train_idx, ]
    x_test  <- x_full[test_idx, ]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # ----- Feature selection SOLO sul training (leakage-free) -----
    vars      <- apply(x_train, 2, var)
    top_genes <- names(sort(vars, decreasing = TRUE))[1:top_k]
    x_train   <- x_train[, top_genes, drop = FALSE]
    x_test    <- x_test[,  top_genes, drop = FALSE]
    
    # ----- Coupling correction opzionale (m0 stimato SOLO sul training) -----
    if (coupling_correct) {
      b_train <- baseline_score[train_idx]
      b_test  <- baseline_score[test_idx]
      
      m0 <- lm(y_train ~ b_train)
      y_train_corr <- resid(m0)
      
      # Applica la stessa correzione al test
      pred_from_baseline <- coef(m0)[1] + coef(m0)[2] * b_test
      y_test_corr <- y_test - pred_from_baseline
      
      y_train <- y_train_corr
      y_test  <- y_test_corr
    }
    
    # ----- Costruzione DMatrix -----
    dtrain <- xgb.DMatrix(as.matrix(x_train), label = y_train)
    dtest  <- xgb.DMatrix(as.matrix(x_test),  label = y_test)
    
    # ----- CV interna per selezione di nrounds -----
    set.seed(123 + i)
    cv_fit <- xgb.cv(
      params                = xgb_params,
      data                  = dtrain,
      nrounds               = max_nrounds,
      nfold                 = 5,
      early_stopping_rounds = 20,
      verbose               = 0
    )
    best_nrounds <- get_best_iter(cv_fit)
    
    # ----- Refit sul training intero con nrounds ottimo -----
    fit <- xgb.train(
      params  = xgb_params,
      data    = dtrain,
      nrounds = best_nrounds,
      verbose = 0
    )
    
    pred <- predict(fit, dtest)
    all_pred <- c(all_pred, pred)
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Preparazione dati ----------

x_baseline_matrix <- t(expr_baseline)
rownames(x_baseline_matrix) <- sub("_NT$", "", rownames(x_baseline_matrix))

baseline_score_clean <- baseline_score
names(baseline_score_clean) <- sub("_NT$", "", names(baseline_score_clean))
baseline_score_aligned <- baseline_score_clean[rownames(x_baseline_matrix)]

y_cont_aligned <- y_cont[rownames(x_baseline_matrix)]

# Verifiche
stopifnot(identical(names(baseline_score_aligned), rownames(x_baseline_matrix)))
stopifnot(identical(names(y_cont_aligned), rownames(x_baseline_matrix)))
stopifnot(sum(is.na(baseline_score_aligned)) == 0)


# ---------- Run: 4 combinazioni (K ∈ {3000, 300} × {leakage-free, fully corrected}) ----------

cat("\n--- XGBoost (K=3000), LEAKAGE-FREE ONLY ---\n")
R2_xgb_3000_leakfree <- evaluate_xgboost_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 3000,
  coupling_correct = FALSE
)
cat("  R² =", round(R2_xgb_3000_leakfree, 4), "\n")

cat("\n--- XGBoost (K=3000), FULLY CORRECTED ---\n")
R2_xgb_3000_fully_corrected <- evaluate_xgboost_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 3000,
  baseline_score   = baseline_score_aligned,
  coupling_correct = TRUE
)
cat("  R² =", round(R2_xgb_3000_fully_corrected, 4), "\n")

cat("\n--- XGBoost (K=300), LEAKAGE-FREE ONLY ---\n")
R2_xgb_300_leakfree <- evaluate_xgboost_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 300,
  coupling_correct = FALSE
)
cat("  R² =", round(R2_xgb_300_leakfree, 4), "\n")

cat("\n--- XGBoost (K=300), FULLY CORRECTED ---\n")
R2_xgb_300_fully_corrected <- evaluate_xgboost_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 300,
  baseline_score   = baseline_score_aligned,
  coupling_correct = TRUE
)
cat("  R² =", round(R2_xgb_300_fully_corrected, 4), "\n")


# ---------- Riepilogo -----------
results_xgb <- data.frame(
  K                 = c(3000, 300),
  Leakage_free_only = c(R2_xgb_3000_leakfree, R2_xgb_300_leakfree),
  Fully_corrected   = c(R2_xgb_3000_fully_corrected, R2_xgb_300_fully_corrected),
  Delta_coupling    = c(R2_xgb_3000_leakfree - R2_xgb_3000_fully_corrected,
                        R2_xgb_300_leakfree  - R2_xgb_300_fully_corrected)
)
cat("\n=== XGBoost (top-K variable genes) ===\n")
print(results_xgb, row.names = FALSE)
#   K Leakage_free_only Fully_corrected Delta_coupling
#3000         0.2812973       0.1237338      0.1575635
# 300         0.1989522       0.0604151      0.1385371

######################################################
####  Proviamo a fare un confronto tra Elastic Net vs sPLS
##############################################################################
# sPLS — NAÏVE (grouped patient-level CV + FS pre-CV, no coupling correction)
#
# Definizione naïve stretta:
#   - Feature selection top-K per varianza fatta UNA VOLTA sull'intero dataset
#   - CV grouped patient-level (via groupKFold, non random-per-sample)
#   - Nessuna coupling correction
#
# Nota: qui sposto la feature selection FUORI dal loop per ottenere la vera naïve.
##############################################################################
library(mixOmics)
# ---------- Funzione naïve ----------
evaluate_spls_naive <- function(x_full, y, folds, top_k,
                                keepX_grid = c(10, 25, 50, 100, 200),
                                ncomp      = 2,
                                inner_folds = 5) {
  
  # NAÏVE: feature selection UNA VOLTA sull'intero dataset (leakage per costruzione)
  vars       <- apply(x_full, 2, var)
  top_genes  <- names(sort(vars, decreasing = TRUE))[1:top_k]
  x_selected <- x_full[, top_genes, drop = FALSE]
  
  # Filtra keepX_grid ai valori sensati (non più di top_k)
  keepX_grid <- keepX_grid[keepX_grid <= top_k]
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]                                      # training indices
    test_idx  <- setdiff(seq_len(nrow(x_selected)), train_idx)
    
    x_train <- x_selected[train_idx, ]
    x_test  <- x_selected[test_idx, ]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # Tuning di keepX via CV interna sul training fold
    set.seed(123 + i)
    tune <- tune.spls(
      X          = x_train,
      Y          = y_train,
      test.keepX = keepX_grid,
      folds      = inner_folds,
      measure    = "R2",
      validation = "Mfold",
      ncomp      = ncomp
    )
    
    fit <- spls(
      X     = x_train,
      Y     = y_train,
      ncomp = ncomp,
      keepX = tune$choice.keepX
    )
    
    pred <- predict(fit, newdata = x_test)
    pred <- pred$predict[, 1, ncomp]   # predizione al componente finale
    
    all_pred <- c(all_pred, pred)
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Preparazione dati ----------
x_baseline_matrix <- t(expr_baseline)
rownames(x_baseline_matrix) <- sub("_NT$", "", rownames(x_baseline_matrix))
y_cont_aligned <- y_cont[rownames(x_baseline_matrix)]

stopifnot(identical(names(y_cont_aligned), rownames(x_baseline_matrix)))


# ---------- Run ----------
cat("\n--- sPLS (K=3000), NAÏVE ---\n")
R2_spls_3000_naive <- evaluate_spls_naive(
  x_full = x_baseline_matrix,
  y      = y_cont_aligned,
  folds  = folds,
  top_k  = 3000
)
cat("  R² =", round(R2_spls_3000_naive, 4), "\n")
R2_spls_3000_naive
# R² = 0.43
cat("\n--- sPLS (K=300), NAÏVE ---\n")
R2_spls_300_naive <- evaluate_spls_naive(
  x_full = x_baseline_matrix,
  y      = y_cont_aligned,
  folds  = folds,
  top_k  = 300,
  keepX_grid = c(5, 10, 25, 50, 100, 200, 300)  # grid adattato per K minore
)
cat("  R² =", round(R2_spls_300_naive, 4), "\n")
R2_spls_300_naive
# R² = 0.25

##############################################################################
# sPLS — LEAKAGE-FREE + (opzionale) COUPLING-CORRECTED
#
# Predittori: top-K geni per varianza, selezionati DENTRO ogni training fold
# Tuning di keepX via CV interna 5-fold sul training set
#
# Fix 1: groupKFold ritorna TRAINING indices
# Fix 2: coupling correction opzionale (m0 stimato sul training, applicata al test)
# Fix 3: set.seed per fold per determinismo del tuning interno
##############################################################################

library(mixOmics)

# ---------- Funzione riusabile ----------
evaluate_spls_leakfree <- function(x_full, y, folds, top_k,
                                   baseline_score  = NULL,
                                   keepX_grid      = c(10, 25, 50, 100, 200),
                                   ncomp           = 2,
                                   inner_folds     = 5,
                                   coupling_correct = FALSE) {
  
  # Sanity check
  if (coupling_correct) {
    if (is.null(baseline_score)) {
      stop("coupling_correct = TRUE richiede baseline_score in input.")
    }
    stopifnot(length(baseline_score) == nrow(x_full))
    stopifnot(!is.null(names(baseline_score)))
    stopifnot(identical(names(baseline_score), rownames(x_full)))
  }
  
  # Filtra keepX_grid ai valori sensati (non più di top_k)
  keepX_grid <- keepX_grid[keepX_grid <= top_k]
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    # FIX: folds[[i]] sono TRAINING indices
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(x_full)), train_idx)
    
    x_train <- x_full[train_idx, ]
    x_test  <- x_full[test_idx, ]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # ----- Feature selection SOLO sul training (leakage-free) -----
    vars      <- apply(x_train, 2, var)
    top_genes <- names(sort(vars, decreasing = TRUE))[1:top_k]
    x_train   <- x_train[, top_genes, drop = FALSE]
    x_test    <- x_test[,  top_genes, drop = FALSE]
    
    # ----- Coupling correction opzionale (m0 stimato SOLO sul training) -----
    if (coupling_correct) {
      b_train <- baseline_score[train_idx]
      b_test  <- baseline_score[test_idx]
      
      m0 <- lm(y_train ~ b_train)
      y_train_corr <- resid(m0)
      
      # Applica la stessa correzione al test
      pred_from_baseline <- coef(m0)[1] + coef(m0)[2] * b_test
      y_test_corr <- y_test - pred_from_baseline
      
      y_train <- y_train_corr
      y_test  <- y_test_corr
    }
    
    # ----- Tuning di keepX via CV interna sul training fold -----
    set.seed(123 + i)
    tune <- tune.spls(
      X          = x_train,
      Y          = y_train,
      test.keepX = keepX_grid,
      folds      = inner_folds,
      measure    = "R2",
      validation = "Mfold",
      ncomp      = ncomp
    )
    
    fit <- spls(
      X     = x_train,
      Y     = y_train,
      ncomp = ncomp,
      keepX = tune$choice.keepX
    )
    
    # Predizione al componente finale
    pred <- predict(fit, newdata = x_test)
    pred <- pred$predict[, 1, ncomp]
    
    all_pred <- c(all_pred, pred)
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Preparazione dati ----------

x_baseline_matrix <- t(expr_baseline)
rownames(x_baseline_matrix) <- sub("_NT$", "", rownames(x_baseline_matrix))

baseline_score_clean <- baseline_score
names(baseline_score_clean) <- sub("_NT$", "", names(baseline_score_clean))
baseline_score_aligned <- baseline_score_clean[rownames(x_baseline_matrix)]

y_cont_aligned <- y_cont[rownames(x_baseline_matrix)]

# Verifiche
stopifnot(identical(names(baseline_score_aligned), rownames(x_baseline_matrix)))
stopifnot(identical(names(y_cont_aligned), rownames(x_baseline_matrix)))
stopifnot(sum(is.na(baseline_score_aligned)) == 0)


# ---------- Run: 4 combinazioni (K ∈ {3000, 300} × {leakage-free, fully corrected}) ----------

cat("\n--- sPLS (K=3000), LEAKAGE-FREE ONLY ---\n")
R2_spls_3000_leakfree <- evaluate_spls_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 3000,
  coupling_correct = FALSE
)
cat("  R² =", round(R2_spls_3000_leakfree, 4), "\n")

cat("\n--- sPLS (K=3000), FULLY CORRECTED ---\n")
R2_spls_3000_fully_corrected <- evaluate_spls_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 3000,
  baseline_score   = baseline_score_aligned,
  coupling_correct = TRUE
)
cat("  R² =", round(R2_spls_3000_fully_corrected, 4), "\n")

cat("\n--- sPLS (K=300), LEAKAGE-FREE ONLY ---\n")
R2_spls_300_leakfree <- evaluate_spls_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 300,
  keepX_grid       = c(5, 10, 25, 50, 100, 200, 300),
  coupling_correct = FALSE
)
cat("  R² =", round(R2_spls_300_leakfree, 4), "\n")

cat("\n--- sPLS (K=300), FULLY CORRECTED ---\n")
R2_spls_300_fully_corrected <- evaluate_spls_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 300,
  keepX_grid       = c(5, 10, 25, 50, 100, 200, 300),
  baseline_score   = baseline_score_aligned,
  coupling_correct = TRUE
)
cat("  R² =", round(R2_spls_300_fully_corrected, 4), "\n")

# ---------- Riepilogo -----------
results_spls <- data.frame(
  K                 = c(3000, 300),
  Leakage_free_only = c(R2_spls_3000_leakfree, R2_spls_300_leakfree),
  Fully_corrected   = c(R2_spls_3000_fully_corrected, R2_spls_300_fully_corrected),
  Delta_coupling    = c(R2_spls_3000_leakfree - R2_spls_3000_fully_corrected,
                        R2_spls_300_leakfree  - R2_spls_300_fully_corrected)
)
cat("\n=== sPLS (top-K variable genes) ===\n")
print(results_spls, row.names = FALSE)
#   K Leakage_free_only Fully_corrected Delta_coupling
#3000         0.4463043       0.2425411      0.2037632
# 300         0.2530874       0.1377912      0.1152962
# sPLS non sta migliorando Elastic Net.


##############################################################################
# SVR (linear kernel) — NAÏVE (grouped patient-level CV + FS pre-CV,
# no coupling correction)
#
# Definizione naïve stretta:
#   - Feature selection top-K per varianza fatta UNA VOLTA sull'intero dataset
#   - Scaling calcolato sull'intero dataset (parte del regime naïve)
#   - CV grouped patient-level (via groupKFold)
#   - Nessuna coupling correction
#
# Nota: nel tuo codice originale scaling e FS erano DENTRO il loop = leakage-free.
# Qui li sposto FUORI per ottenere la vera naïve.
##############################################################################
library(e1071)
# ---------- Funzione naïve ----------
evaluate_svr_naive <- function(x_full, y, folds, top_k,
                               cost_grid = c(0.01, 0.1, 1, 10, 100)) {
  
  # NAÏVE step 1: feature selection UNA VOLTA sull'intero dataset (leakage)
  vars       <- apply(x_full, 2, var)
  top_genes  <- names(sort(vars, decreasing = TRUE))[1:top_k]
  x_selected <- x_full[, top_genes, drop = FALSE]
  
  # NAÏVE step 2: scaling calcolato UNA VOLTA sull'intero dataset (leakage)
  full_mean <- apply(x_selected, 2, mean)
  full_sd   <- apply(x_selected, 2, sd)
  full_sd[full_sd == 0] <- 1
  x_scaled  <- scale(x_selected, center = full_mean, scale = full_sd)
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    train_idx <- folds[[i]]                                      # training indices
    test_idx  <- setdiff(seq_len(nrow(x_scaled)), train_idx)
    
    x_train <- x_scaled[train_idx, ]
    x_test  <- x_scaled[test_idx, ]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # Tuning di cost sul training fold
    best_cost <- NULL
    best_rmse <- Inf
    for (cc in cost_grid) {
      set.seed(123 + i)
      tune_fit <- tune(
        svm,
        train.x = x_train,
        train.y = y_train,
        kernel  = "linear",
        ranges  = list(cost = cc)
      )
      rmse <- tune_fit$best.performance
      if (rmse < best_rmse) {
        best_rmse <- rmse
        best_cost <- cc
      }
    }
    
    fit <- svm(
      x      = x_train,
      y      = y_train,
      type   = "eps-regression",
      kernel = "linear",
      cost   = best_cost,
      scale  = FALSE   # già scalato fuori dal loop
    )
    
    pred <- predict(fit, x_test)
    all_pred <- c(all_pred, pred)
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Preparazione dati ----------

x_baseline_matrix <- t(expr_baseline)
rownames(x_baseline_matrix) <- sub("_NT$", "", rownames(x_baseline_matrix))
y_cont_aligned <- y_cont[rownames(x_baseline_matrix)]

stopifnot(identical(names(y_cont_aligned), rownames(x_baseline_matrix)))


# ---------- Run ----------

cat("\n--- SVR (K=3000), NAÏVE ---\n")
R2_svr_3000_naive <- evaluate_svr_naive(
  x_full = x_baseline_matrix,
  y      = y_cont_aligned,
  folds  = folds,
  top_k  = 3000
)
cat("  R² =", round(R2_svr_3000_naive, 4), "\n")
R2_svr_3000_naive # 0.38 3000 genes

cat("\n--- SVR (K=300), NAÏVE ---\n")
R2_svr_300_naive <- evaluate_svr_naive(
  x_full = x_baseline_matrix,
  y      = y_cont_aligned,
  folds  = folds,
  top_k  = 300
)
cat("  R² =", round(R2_svr_300_naive, 4), "\n")
R2_svr_300_naive # 0.24 300 genes

##############################################################################
# SVR (linear kernel) — LEAKAGE-FREE + (opzionale) COUPLING-CORRECTED
#
# Predittori: top-K geni per varianza, selezionati DENTRO ogni training fold
# Scaling: media/sd calcolate SOLO sul training fold, applicate al test
# Tuning di cost via CV interna sul training set (via tune() di e1071)
#
# Fix 1: groupKFold ritorna TRAINING indices
# Fix 2: coupling correction opzionale (m0 stimato sul training, applicata al test)
# Fix 3: set.seed per fold per riproducibilità del tuning interno
##############################################################################

library(e1071)

# ---------- Funzione riusabile ----------
evaluate_svr_leakfree <- function(x_full, y, folds, top_k,
                                  baseline_score  = NULL,
                                  cost_grid       = c(0.01, 0.1, 1, 10, 100),
                                  coupling_correct = FALSE) {
  
  # Sanity check
  if (coupling_correct) {
    if (is.null(baseline_score)) {
      stop("coupling_correct = TRUE richiede baseline_score in input.")
    }
    stopifnot(length(baseline_score) == nrow(x_full))
    stopifnot(!is.null(names(baseline_score)))
    stopifnot(identical(names(baseline_score), rownames(x_full)))
  }
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    # FIX: folds[[i]] sono TRAINING indices
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(x_full)), train_idx)
    
    x_train <- x_full[train_idx, ]
    x_test  <- x_full[test_idx, ]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # ----- Feature selection SOLO sul training (leakage-free) -----
    vars      <- apply(x_train, 2, var)
    top_genes <- names(sort(vars, decreasing = TRUE))[1:top_k]
    x_train   <- x_train[, top_genes, drop = FALSE]
    x_test    <- x_test[,  top_genes, drop = FALSE]
    
    # ----- Scaling SOLO sul training (leakage-free) -----
    train_mean <- apply(x_train, 2, mean)
    train_sd   <- apply(x_train, 2, sd)
    train_sd[train_sd == 0] <- 1
    x_train    <- scale(x_train, center = train_mean, scale = train_sd)
    x_test     <- scale(x_test,  center = train_mean, scale = train_sd)
    
    # ----- Coupling correction opzionale (m0 stimato SOLO sul training) -----
    if (coupling_correct) {
      b_train <- baseline_score[train_idx]
      b_test  <- baseline_score[test_idx]
      
      m0 <- lm(y_train ~ b_train)
      y_train_corr <- resid(m0)
      
      # Applica la stessa correzione al test
      pred_from_baseline <- coef(m0)[1] + coef(m0)[2] * b_test
      y_test_corr <- y_test - pred_from_baseline
      
      y_train <- y_train_corr
      y_test  <- y_test_corr
    }
    
    # ----- Tuning di cost sul training fold -----
    best_cost <- NULL
    best_rmse <- Inf
    for (cc in cost_grid) {
      set.seed(123 + i)
      tune_fit <- tune(
        svm,
        train.x = x_train,
        train.y = y_train,
        kernel  = "linear",
        ranges  = list(cost = cc)
      )
      rmse <- tune_fit$best.performance
      if (rmse < best_rmse) {
        best_rmse <- rmse
        best_cost <- cc
      }
    }
    
    fit <- svm(
      x      = x_train,
      y      = y_train,
      type   = "eps-regression",
      kernel = "linear",
      cost   = best_cost,
      scale  = FALSE   # già scalato manualmente
    )
    
    pred <- predict(fit, x_test)
    all_pred <- c(all_pred, pred)
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Preparazione dati ----------

x_baseline_matrix <- t(expr_baseline)
rownames(x_baseline_matrix) <- sub("_NT$", "", rownames(x_baseline_matrix))

baseline_score_clean <- baseline_score
names(baseline_score_clean) <- sub("_NT$", "", names(baseline_score_clean))
baseline_score_aligned <- baseline_score_clean[rownames(x_baseline_matrix)]

y_cont_aligned <- y_cont[rownames(x_baseline_matrix)]

# Verifiche
stopifnot(identical(names(baseline_score_aligned), rownames(x_baseline_matrix)))
stopifnot(identical(names(y_cont_aligned), rownames(x_baseline_matrix)))
stopifnot(sum(is.na(baseline_score_aligned)) == 0)


# ---------- Run: 4 combinazioni (K ∈ {3000, 300} × {leakage-free, fully corrected}) ----------

cat("\n--- SVR (K=3000), LEAKAGE-FREE ONLY ---\n")
R2_svr_3000_leakfree <- evaluate_svr_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 3000,
  coupling_correct = FALSE
)
cat("  R² =", round(R2_svr_3000_leakfree, 4), "\n")

cat("\n--- SVR (K=3000), FULLY CORRECTED ---\n")
R2_svr_3000_fully_corrected <- evaluate_svr_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 3000,
  baseline_score   = baseline_score_aligned,
  coupling_correct = TRUE
)
cat("  R² =", round(R2_svr_3000_fully_corrected, 4), "\n")

cat("\n--- SVR (K=300), LEAKAGE-FREE ONLY ---\n")
R2_svr_300_leakfree <- evaluate_svr_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 300,
  coupling_correct = FALSE
)
cat("  R² =", round(R2_svr_300_leakfree, 4), "\n")

cat("\n--- SVR (K=300), FULLY CORRECTED ---\n")
R2_svr_300_fully_corrected <- evaluate_svr_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  top_k            = 300,
  baseline_score   = baseline_score_aligned,
  coupling_correct = TRUE
)
cat("  R² =", round(R2_svr_300_fully_corrected, 4), "\n")


# ---------- Riepilogo -----------
results_svr <- data.frame(
  K                 = c(3000, 300),
  Leakage_free_only = c(R2_svr_3000_leakfree, R2_svr_300_leakfree),
  Fully_corrected   = c(R2_svr_3000_fully_corrected, R2_svr_300_fully_corrected),
  Delta_coupling    = c(R2_svr_3000_leakfree - R2_svr_3000_fully_corrected,
                        R2_svr_300_leakfree  - R2_svr_300_fully_corrected)
)
cat("\n=== SVR (top-K variable genes) ===\n")
print(results_svr, row.names = FALSE)
#   K Leakage_free_only Fully_corrected Delta_coupling
#3000         0.3690418       0.2209743      0.1480676
# 300         0.2283274       0.1032604      0.1250670

#SVR non sta migliorando Elastic Net.
#Probabilmente il limite non è l'algoritmo, ma è l'informazione contenuta nei geni baseline. non sembra esserci un modello non lineare nascosto che Elastic Net non riesce a catturare.
# Multiple machine-learning algorithms were benchmarked under grouped patient-level cross-validation, including Elastic Net, sparse Partial Least Squares, 
# Support Vector Regression, and XGBoost. Elastic Net consistently achieved the highest out-of-sample performance, while nonlinear approaches 
# such as XGBoost performed substantially worse (CV-R² ≈ 0.12–0.16). These findings indicate that cetuximab-induced pathway suppression is primarily encoded 
# through a distributed and approximately linear transcriptomic signal rather than complex nonlinear interactions.
###########################################################################################

#######################################################################################
####################### NEW MODEL
# Now, We focus only on 569  EGFR pathway genes
egfr_gene_list <- unique(unlist(gene_sets))
egfr_gene_list <- intersect(egfr_gene_list, rownames(expr_human)) # only human
dim(expr_baseline)
dim((expr_human))
#rownames(expr_baseline) <- sub("^H_", "", rownames(expr_baseline))
length(egfr_gene_list) # 569 genes
x_genes_in_pathways <- t(expr_baseline[egfr_gene_list,]) 
dim(x_genes_in_pathways)[2]  # 569 predictors
folds <- group_folds(groups, k = 3, repeats = 20)
ctrl <- trainControl(
  method = "cv",
  index = folds,
  savePredictions = "final",
)
#Controllo allineamento y_cont: fai sempre:
stopifnot(all(names(y_cont) == rownames(x_genes_in_pathways)))
fit_reg_genes_in_the_pathways <- caret::train(x = x_genes_in_pathways,y = y_cont,method = "glmnet",trControl = ctrl,tuneLength = 20,tuneGrid = expand.grid(alpha=seq(0,0.2,length=11),lambda = 0.1))
sd(y_cont)# 0.1802557
min(fit_reg_genes_in_the_pathways$results$RMSE)# 0.1256245
min(fit_reg_genes_in_the_pathways$results$RMSE) < sd(y_cont) # If RMSE < SD → model captures signal.
subset(fit_reg_genes_in_the_pathways$results, fit_reg_genes_in_the_pathways$results$RMSE == min(fit_reg_genes_in_the_pathways$results$RMSE))
#    alpha lambda      RMSE  Rsquared       MAE      RMSESD RsquaredSD       MAESD
#400  0.19    0.1 0.1256245 0.5471715 0.09892583 0.01514051 0.08768764 0.01061928

stopifnot(all(names(resid_response) == rownames(x_genes_in_pathways)))
fit_reg_genes_in_the_pathways_resid<-caret::train(x=x_genes_in_pathways,y=resid_response,method="glmnet",trControl=ctrl,tuneLength=20,tuneGrid=expand.grid(alpha=seq(0,0.2,length=11),lambda=0.1))
subset(fit_reg_genes_in_the_pathways_resid$results, fit_reg_genes_in_the_pathways_resid$results$RMSE == min(fit_reg_genes_in_the_pathways_resid$results$RMSE))
#  alpha lambda      RMSE  Rsquared        MAE     RMSESD RsquaredSD      MAESD
#8  0.14    0.1 0.1251917 0.2498665 0.09965073 0.01714054  0.1327762 0.01223854

############# solo i primi 300 out of 569 genes

x_genes_in_pathways <- t(expr_baseline[egfr_gene_list,]) 
vars <- apply(x_genes_in_pathways , 2, var)
top_300_genes = names(sort(vars, decreasing = TRUE))[1:300]
#Controllo allineamento y_cont: fai sempre:
stopifnot(all(names(y_cont) == rownames(top_300_genes)))
fit_reg_genes_top300  <- caret::train(x = x_genes_in_pathways[,top_300_genes],y = y_cont,method = "glmnet",trControl = ctrl,tuneLength = 20,tuneGrid = expand.grid(alpha=seq(0,0.2,length=11),lambda = 0.1))
sd(y_cont)# 0.1802557
min(fit_reg_genes_top300$results$RMSE)# 0.1256245
min(fit_reg_genes_top300$results$RMSE) < sd(y_cont) # If RMSE < SD → model captures signal.
subset(fit_reg_genes_top300$results, fit_reg_genes_top300$results$RMSE == min(fit_reg_genes_top300$results$RMSE))
#  alpha lambda      RMSE  Rsquared        MAE     RMSESD RsquaredSD      MAESD
#8  0.14    0.1 0.1266105 0.5270107 0.09942408 0.01586952  0.1036775 0.01298976

#Controllo allineamento y_cont: fai sempre:
stopifnot(all(names(resid_response) == rownames(top_300_genes)))
fit_reg_genes_in_the_pathways_resid <- caret::train(x = x_genes_in_pathways[,top_300_genes], y=resid_response,method="glmnet",trControl=ctrl,tuneLength=20,tuneGrid=expand.grid(alpha=seq(0,0.2,length=11),lambda=0.1))
subset(fit_reg_genes_in_the_pathways_resid$results, fit_reg_genes_in_the_pathways_resid$results$RMSE == min(fit_reg_genes_in_the_pathways_resid$results$RMSE))
#  alpha lambda      RMSE  Rsquared        MAE     RMSESD RsquaredSD     MAESD
#6   0.1    0.1 0.1232263 0.2666583 0.09781972 0.01346443  0.1276145 0.0130203

############# LEAKING REMOVAL
# NOTE: nessuna feature selection per questo modello — i 569 geni EGFR-axis sono definiti a priori da MSigDB Hallmark gene sets. Il loop CV serve solo
# a calcolare predizioni out-of-fold con tuning di (alpha, lambda) per fold. Non c'è "leakage removal" perché non c'è feature selection da cui rimuoverlo.
#For models with a priori predictor sets (EGFR-axis genes, n=569; GSVA Hallmark pathways, n=50; principal components, n=20), 
#feature selection is not performed and the leakage-free correction is not applicable. For these models, naïve and leakage-free 
#out-of-sample R² are identical by construction; only the coupling-corrected estimate is reported as a distinct regime."

# Se vuoi rispondere a un revisore che dice "ok, ma allora come fai a sapere che 569 geni EGFR-axis non sono troppo pochi/troppi?", 
# puoi aggiungere una feature selection interna anche dentro questi pathway-restricted models. 
# L'idea: dai 569 geni, seleziona i top-K più variabili dentro il fold di training e vedi se il segnale regge.
# Versione con feature selection INTERNA anche sui modelli pathway-restricted. Per sondare la robustezza del segnale a diversi K
##############################################################################
# fit_reg_genes (EGFR-axis, 569 geni) + SENSITIVITY ON K
# Due varianti: leakage-free-only e fully-corrected (leakage + coupling)
#
# Fix 1: groupKFold ritorna TRAINING indices
# Fix 2: alpha tuning via inner CV
# Fix 3: uso nrow(x_full) invece di nrow(x) (bug del codice originale)
# Fix 4: coupling correction opzionale (stima lm sul training, applicata al test)
##############################################################################

# ---------- Funzione riusabile ----------
evaluate_with_inner_fs <- function(x_full, y, folds, top_k,
                                   baseline_score  = NULL,
                                   alpha_grid      = seq(0, 0.2, length = 11),
                                   coupling_correct = FALSE) {
  
  # Sanity check
  if (coupling_correct) {
    if (is.null(baseline_score)) {
      stop("coupling_correct = TRUE richiede baseline_score in input.")
    }
    stopifnot(length(baseline_score) == nrow(x_full))
    stopifnot(!is.null(names(baseline_score)))
    stopifnot(identical(names(baseline_score), rownames(x_full)))
  }
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    # FIX: folds[[i]] sono TRAINING indices, non test
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(x_full)), train_idx)   # FIX: x_full, non x
    
    x_train <- x_full[train_idx, ]
    x_test  <- x_full[test_idx, ]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # ----- Feature selection INTERNA: top-K per varianza nel training fold -----
    if (top_k < ncol(x_train)) {
      vars      <- apply(x_train, 2, var)
      top_feats <- names(sort(vars, decreasing = TRUE))[1:top_k]
      x_train   <- x_train[, top_feats, drop = FALSE]
      x_test    <- x_test[,  top_feats, drop = FALSE]
    }
    
    # ----- Coupling correction opzionale (m0 stimato SOLO sul training) -----
    if (coupling_correct) {
      b_train <- baseline_score[train_idx]
      b_test  <- baseline_score[test_idx]
      
      m0 <- lm(y_train ~ b_train)
      y_train_corr <- resid(m0)
      
      # Applica la stessa correzione al test
      pred_from_baseline <- coef(m0)[1] + coef(m0)[2] * b_test
      y_test_corr <- y_test - pred_from_baseline
      
      y_train <- y_train_corr
      y_test  <- y_test_corr
    }
    
    # ----- Tuning di alpha via inner CV -----
    set.seed(123 + i)
    inner_foldid <- sample(rep(1:5, length.out = length(y_train)))
    
    cv_results <- lapply(alpha_grid, function(a) {
      cvfit <- cv.glmnet(x_train, y_train, alpha = a, foldid = inner_foldid)
      list(alpha = a, cvm_min = min(cvfit$cvm), cvfit = cvfit)
    })
    best_idx <- which.min(sapply(cv_results, function(r) r$cvm_min))
    best_fit <- cv_results[[best_idx]]$cvfit
    
    pred <- predict(best_fit, newx = x_test, s = "lambda.min")
    all_pred <- c(all_pred, as.numeric(pred))
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Preparazione dati ----------

# Matrice predittori: 569 geni EGFR-axis
x_genes_in_pathways <- t(expr_baseline[egfr_gene_list, ])

# Pulizia rownames (rimuove il suffisso _NT)
rownames(x_genes_in_pathways) <- sub("_NT$", "", rownames(x_genes_in_pathways))

# baseline_score con nomi puliti, allineato
baseline_score_clean <- baseline_score
names(baseline_score_clean) <- sub("_NT$", "", names(baseline_score_clean))
baseline_score_aligned <- baseline_score_clean[rownames(x_genes_in_pathways)]

# y_cont allineato
y_cont_aligned <- y_cont[rownames(x_genes_in_pathways)]

# Verifiche
stopifnot(identical(names(baseline_score_aligned), rownames(x_genes_in_pathways)))
stopifnot(identical(names(y_cont_aligned), rownames(x_genes_in_pathways)))
stopifnot(sum(is.na(baseline_score_aligned)) == 0)


# ---------- Sensitivity su K per il modello EGFR-axis: DUE varianti ----------

K_grid <- c(50, 100, 200, 300, 569)

# Variante 1: leakage-free only
cat("\n--- EGFR-axis, LEAKAGE-FREE ONLY ---\n")
R2_egfr_axis_leakfree <- sapply(K_grid, function(k) {
  cat("  K =", k, "...\n")
  evaluate_with_inner_fs(
    x_full           = x_genes_in_pathways,
    y                = y_cont_aligned,
    folds            = folds,
    top_k            = k,
    coupling_correct = FALSE
  )
})
names(R2_egfr_axis_leakfree) <- paste0("top_", K_grid)

# Variante 2: fully corrected (leakage + coupling)
cat("\n--- EGFR-axis, FULLY CORRECTED ---\n")
R2_egfr_axis_fully_corrected <- sapply(K_grid, function(k) {
  cat("  K =", k, "...\n")
  evaluate_with_inner_fs(
    x_full           = x_genes_in_pathways,
    y                = y_cont_aligned,
    folds            = folds,
    top_k            = k,
    baseline_score   = baseline_score_aligned,
    coupling_correct = TRUE
  )
})
names(R2_egfr_axis_fully_corrected) <- paste0("top_", K_grid)


# ---------- Riepilogo EGFR-axis -----------
results_egfr_axis <- data.frame(
  K                 = K_grid,
  Leakage_free_only = R2_egfr_axis_leakfree,
  Fully_corrected   = R2_egfr_axis_fully_corrected,
  Delta_coupling    = R2_egfr_axis_leakfree - R2_egfr_axis_fully_corrected
)
cat("\n=== EGFR-axis genes (569 a priori pool) ===\n")
print(results_egfr_axis, row.names = FALSE)
#  K Leakage_free_only Fully_corrected Delta_coupling
# 50         0.1766796       0.1446778     0.03200181
#100         0.3123345       0.1617299     0.15060463
#200         0.4693903       0.1961695     0.27322085
#300         0.4841986       0.1961092     0.28808938
#569         0.4778790       0.2230838     0.25479520

# ---------- Sensitivity aggiuntiva sui 50 GSVA pathway (senza esclusione) ----------
# NOTE: qui NON escludiamo E2F/KRAS/mTORC1, quindi il coupling matematico strutturale
#        è massimo. La correzione via baseline_score è particolarmente incisiva.

x_pathway_all <- t(gsva_nt_all)
rownames(x_pathway_all) <- sub("_NT$", "", rownames(x_pathway_all))
x_pathway_all <- x_pathway_all[rownames(x_genes_in_pathways), ]

stopifnot(identical(rownames(x_pathway_all), rownames(x_genes_in_pathways)))

K_grid_pw <- c(5, 10, 20, 50)

# Variante 1: leakage-free only
cat("\n--- 50 pathway (incl. coupled), LEAKAGE-FREE ONLY ---\n")
R2_pathway_leakfree <- sapply(K_grid_pw, function(k) {
  cat("  K =", k, "...\n")
  evaluate_with_inner_fs(
    x_full           = x_pathway_all,
    y                = y_cont_aligned,
    folds            = folds,
    top_k            = k,
    coupling_correct = FALSE
  )
})
names(R2_pathway_leakfree) <- paste0("top_", K_grid_pw)

# Variante 2: fully corrected
cat("\n--- 50 pathway (incl. coupled), FULLY CORRECTED ---\n")
R2_pathway_fully_corrected <- sapply(K_grid_pw, function(k) {
  cat("  K =", k, "...\n")
  evaluate_with_inner_fs(
    x_full           = x_pathway_all,
    y                = y_cont_aligned,
    folds            = folds,
    top_k            = k,
    baseline_score   = baseline_score_aligned,
    coupling_correct = TRUE
  )
})
names(R2_pathway_fully_corrected) <- paste0("top_", K_grid_pw)


# ---------- Riepilogo 50 pathway (incl. coupled) -----------
results_pathway_all <- data.frame(
  K                 = K_grid_pw,
  Leakage_free_only = R2_pathway_leakfree,
  Fully_corrected   = R2_pathway_fully_corrected,
  Delta_coupling    = R2_pathway_leakfree - R2_pathway_fully_corrected
)
cat("\n=== 50 GSVA pathways (E2F/KRAS/mTORC1 INCLUDED) ===\n")
print(results_pathway_all, row.names = FALSE)
# K Leakage_free_only Fully_corrected Delta_coupling
#5          0.3664593    0.0008822605      0.3655770
#10         0.3746250    0.0786594707      0.2959655
#20         0.4012112    0.0788055496      0.3224057
#50         0.4356299    0.1355571928      0.3000727

###########################################
pred <- fit_reg_genes_in_the_pathways$pred
# Keep only best tuning parameters
best <- fit_reg_genes_in_the_pathways$bestTune
pred_best <- pred[pred$alpha == best$alpha & pred$lambda == best$lambda, ]
# Global R² (across ALL held-out predictions)
R2_global <- cor(pred_best$obs, pred_best$pred)^2
R2_global # 0.50
coefs <- coef(fit_reg_genes_in_the_pathways$finalModel, s = fit_reg_genes_in_the_pathways$bestTune$lambda)
head(coefs)
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df$Feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "Coefficient"
coefs_df <- coefs_df[coefs_df$Coefficient != 0, ]
dim(coefs_df)
top <- coefs_df[order(abs(coefs_df$Coefficient), decreasing = TRUE), ][1:30, ]
top <- top[top$Feature != "(Intercept)", ]
ggplot(top, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "TRUE" = "firebrick",
      "FALSE" = "steelblue"
    ),
    guide = "none"
  ) +
  labs(x = "", y = "Beta Coefficient") +
  ggtitle("Top 30 predictors out of 569 EGFR pathway genes") +
  theme_minimal()
pred <- predict(fit_reg_genes_in_the_pathways,  x_genes_in_pathways)
round(cor(y_cont, pred),2) #  0.89



##########################
# But here I should use the gsva_nt values, not the expression values. As written in the paper.
# Now I use  50 GSVA pathways  as predictors
# NT pathway activity (3 features) → suppression score (derived from same pathways)
# Nel paper metterei il pathway model come analisi principale e il gene model come analisi esplorativa.
t_gsva_nt = t(gsva_nt)
#Correct way to use GSVA as predictors Use full pathway space  
hallmark_all <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets_all <- split(hallmark_all$gene_symbol, hallmark_all$gs_name)
param_all <-   gsvaParam(expr_human, gene_sets_all)
length(gene_sets_all) #50
dim(expr_human)
gsva_all <- gsva(param_all)
# baseline only
gsva_nt_all <- gsva_all[, grepl("_NT$", colnames(gsva_all))] # 50 pathways
dim(gsva_nt_all ) #50 predictors
colnames(gsva_nt_all) <- sub("_NT$", "", colnames(gsva_nt_all)) # 50 pathways
folds <- group_folds(groups, k = 3, repeats = 20)
ctrl <- trainControl(method = "cv",index = folds,savePredictions = "final",)
#Controllo allineamento y_cont: fai sempre:
stopifnot(all(names(y_cont) == rownames(t(gsva_nt_all))))
#Now Predictors are ~50 pathways Outcome = EGFR suppression (3 pathways). No explicit circularity
fit_reg_pathway <- caret::train(x = t(gsva_nt_all),y = y_cont,method = "glmnet",preProcess = c("center", "scale"),trControl = ctrl,
                                tuneLength = 20,tuneGrid = expand.grid(alpha = seq(0.9, 1, length = 11),lambda = seq(0.001, 0.1, length = 20)))
sum(is.na(fit_reg_pathway$results$RMSE))
sd(y_cont)# 0.1802557
min(fit_reg_pathway$results$RMSE)#  0.1302624
min(fit_reg_pathway$results$RMSE) < sd(y_cont) # If RMSE < SD → model captures signal.
subset(fit_reg_pathway$results,fit_reg_pathway$results$RMSE == min(fit_reg_pathway$results$RMSE))
#    alpha      lambda       RMSE  Rsquared        MAE      RMSESD RsquaredSD       MAESD
#3     0.9  0.01142105  0.1302624 0.5003283 0.1054797   0.01504382  0.1001408  0.01113432

stopifnot(all(names(resid_response) == rownames(t(gsva_nt_all))))
fit_reg_pathway_resid <- caret::train(x = t(gsva_nt_all),y = resid_response,method = "glmnet",preProcess = c("center", "scale"),trControl = ctrl,tuneLength = 20,
                                      tuneGrid = expand.grid(alpha = seq(0.9, 1, length = 11),lambda = seq(0.001, 0.1, length = 20)))
subset(fit_reg_pathway_resid$results,fit_reg_pathway_resid$results$RMSE == min(fit_reg_pathway_resid$results$RMSE))
#    alpha     lambda      RMSE  Rsquared       MAE    RMSESD RsquaredSD      MAESD
#205     1 0.02184211 0.1290255 0.2360594 0.1045494 0.0159466  0.1249254 0.01308539

###################################
# fit_reg_pathway_decoupled. only 47 pathways!
# Pathway che entrano nella definizione del suppression score (da escludere)
excluded_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_MTORC1_SIGNALING"
)

# Rimozione: 50 - 3 = 47 pathway candidati
x_full <- t(gsva_nt_all)
x_full <- x_full[, !colnames(x_full) %in% excluded_pathways]
cat("Pathway inclusi nella predizione:", ncol(x_full), "\n")
# Attesi: 47

folds <- group_folds(groups, k = 3, repeats = 20)
ctrl <- trainControl(method = "cv",index = folds,savePredictions = "final",)
#Controllo allineamento y_cont: fai sempre:
stopifnot(all(names(y_cont) == rownames(x_full)))
#Now Predictors are ~50 pathways Outcome = EGFR suppression (3 pathways). No explicit circularity
fit_reg_pathway_decoupled <- caret::train(x = x_full ,y = y_cont,method = "glmnet",preProcess = c("center", "scale"),trControl = ctrl,
                                tuneLength = 20,tuneGrid = expand.grid(alpha = seq(0.9, 1, length = 11),lambda = seq(0.001, 0.1, length = 20)))
sum(is.na(fit_reg_pathway_decoupled$results$RMSE))
sd(y_cont)# 0.1802557
min(fit_reg_pathway_decoupled$results$RMSE)#  0.1302624
min(fit_reg_pathway_decoupled$results$RMSE) < sd(y_cont) # If RMSE < SD → model captures signal.
subset(fit_reg_pathway_decoupled$results,fit_reg_pathway_decoupled$results$RMSE == min(fit_reg_pathway_decoupled$results$RMSE))
#     alpha    lambda      RMSE  Rsquared       MAE     RMSESD RsquaredSD     MAESD
#203     1 0.01142105 0.1349447 0.4547516 0.1079137 0.01344208  0.1258937 0.0105156

stopifnot(all(names(resid_response) == rownames(x_full )))
fit_reg_pathway_decoupled_resid <- caret::train(x = x_full,y = resid_response,method = "glmnet",preProcess = c("center", "scale"),trControl = ctrl,tuneLength = 20,
                                      tuneGrid = expand.grid(alpha = seq(0.9, 1, length = 11),lambda = seq(0.001, 0.1, length = 20)))
subset(fit_reg_pathway_decoupled_resid$results,fit_reg_pathway_decoupled_resid$results$RMSE == min(fit_reg_pathway_decoupled_resid$results$RMSE))
#    alpha     lambda      RMSE  Rsquared       MAE     RMSESD RsquaredSD      MAESD
#206     1 0.02705263 0.1288406 0.2356086 0.1042278 0.01219829  0.1194864 0.01101951

##############################################################################
# LEAKAGE-FREE + (opzionale) COUPLING-CORRECTED EVALUATION 
# OF 47 NON-COUPLED GSVA PATHWAYS + SENSITIVITY ON K
# 
# Predittori: 47 pathway Hallmark, esclusi i 3 pathway che definiscono 
# la suppression outcome (E2F, KRAS, mTORC1).
#
# Fix 1: groupKFold returns TRAINING indices
# Fix 2: alpha tuning via inner CV (not vector to cv.glmnet)
# Fix 3: exclude 3 coupled pathways → coupling STRUTTURALE già rimosso
# Fix 4: coupling correction opzionale via outcome residualization contro
#        il baseline_score (rimuove il coupling RESIDUO da collinearità 
#        fra pathway rimasti e pathway esclusi)
# 
# Sensitivity: K ∈ {10, 20, 30, 40, 47}
##############################################################################
#Il coupling matematico diretto è già stato rimosso escludendo E2F/KRAS/mTORC1 dai predittori. La suppression outcome è definita come mean(NT − CET) 
#sui 3 pathway esclusi, quindi i 47 pathway rimanenti non hanno relazione algebrica diretta con l'outcome.
#Il baseline_score, però, resta la media dei 3 pathway coupled (E2F, KRAS, mTORC1 in condizione NT). 
#Quando residualizzi contro baseline_score, stai rimuovendo la porzione dell'outcome spiegata da quei 3 pathway — non dai 47 predittori.

# ---------- Funzione riusabile per un dato K ----------
evaluate_pathway_leakfree <- function(x_full, y, folds, top_k,
                                      baseline_score  = NULL,
                                      alpha_grid      = seq(0, 1, length = 11),
                                      coupling_correct = FALSE) {
  
  # Sanity check
  if (coupling_correct) {
    if (is.null(baseline_score)) {
      stop("coupling_correct = TRUE richiede baseline_score in input.")
    }
    stopifnot(length(baseline_score) == nrow(x_full))
    stopifnot(!is.null(names(baseline_score)))
    stopifnot(identical(names(baseline_score), rownames(x_full)))
  }
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    # FIX: folds[[i]] sono TRAINING indices
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(x_full)), train_idx)
    
    x_train <- x_full[train_idx, ]
    x_test  <- x_full[test_idx, ]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # ----- Feature selection interna: top-K per varianza dentro il fold -----
    if (top_k < ncol(x_train)) {
      vars      <- apply(x_train, 2, var)
      top_paths <- names(sort(vars, decreasing = TRUE))[1:top_k]
      x_train   <- x_train[, top_paths, drop = FALSE]
      x_test    <- x_test[,  top_paths, drop = FALSE]
    }
    
    # ----- Coupling correction opzionale (m0 stimato SOLO sul training) -----
    if (coupling_correct) {
      b_train <- baseline_score[train_idx]
      b_test  <- baseline_score[test_idx]
      
      m0 <- lm(y_train ~ b_train)
      y_train_corr <- resid(m0)
      
      # Applica la stessa correzione al test
      pred_from_baseline <- coef(m0)[1] + coef(m0)[2] * b_test
      y_test_corr <- y_test - pred_from_baseline
      
      y_train <- y_train_corr
      y_test  <- y_test_corr
    }
    
    # ----- Tuning di alpha via inner CV -----
    set.seed(123 + i)
    inner_foldid <- sample(rep(1:5, length.out = length(y_train)))
    
    cv_results <- lapply(alpha_grid, function(a) {
      cvfit <- cv.glmnet(x_train, y_train, alpha = a, foldid = inner_foldid)
      list(alpha = a, cvm_min = min(cvfit$cvm), cvfit = cvfit)
    })
    best_idx <- which.min(sapply(cv_results, function(r) r$cvm_min))
    best_fit <- cv_results[[best_idx]]$cvfit
    
    pred <- predict(best_fit, newx = x_test, s = "lambda.min")
    all_pred <- c(all_pred, as.numeric(pred))
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Preparazione dati ----------

# Pathway che entrano nella definizione del suppression score (da escludere)
excluded_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_MTORC1_SIGNALING"
)

# 47 pathway rimasti come predittori
x_full <- t(gsva_nt_all)
rownames(x_full) <- sub("_NT$", "", rownames(x_full))
x_full <- x_full[, !colnames(x_full) %in% excluded_pathways]
cat("Pathway inclusi nella predizione:", ncol(x_full), "\n")   # atteso: 47

# baseline_score con nomi puliti, allineato a x_full
baseline_score_clean <- baseline_score
names(baseline_score_clean) <- sub("_NT$", "", names(baseline_score_clean))
baseline_score_aligned <- baseline_score_clean[rownames(x_full)]

# y_cont allineato
y_cont_aligned <- y_cont[rownames(x_full)]

# Verifiche
stopifnot(identical(names(baseline_score_aligned), rownames(x_full)))
stopifnot(identical(names(y_cont_aligned), rownames(x_full)))
stopifnot(sum(is.na(baseline_score_aligned)) == 0)


# ---------- Sensitivity su K: DUE varianti ----------

K_grid_pw <- c(10, 20, 30, 40, 47)

# Variante 1: leakage-free only (già calcolata prima)
cat("\n--- Variante 1: LEAKAGE-FREE ONLY ---\n")
R2_pathway_leakfree_only <- sapply(K_grid_pw, function(k) {
  cat("  K =", k, "...\n")
  evaluate_pathway_leakfree(
    x_full           = x_full,
    y                = y_cont_aligned,
    folds            = folds,
    top_k            = k,
    coupling_correct = FALSE
  )
})
names(R2_pathway_leakfree_only) <- paste0("top_", K_grid_pw)

# Variante 2: fully corrected (leakage-free + coupling-corrected)
cat("\n--- Variante 2: FULLY CORRECTED (leakage + coupling) ---\n")
R2_pathway_fully_corrected <- sapply(K_grid_pw, function(k) {
  cat("  K =", k, "...\n")
  evaluate_pathway_leakfree(
    x_full           = x_full,
    y                = y_cont_aligned,
    folds            = folds,
    top_k            = k,
    baseline_score   = baseline_score_aligned,
    coupling_correct = TRUE
  )
})
names(R2_pathway_fully_corrected) <- paste0("top_", K_grid_pw)


# ---------- Riepilogo -----------
results_table <- data.frame(
  K                   = K_grid_pw,
  Leakage_free_only   = R2_pathway_leakfree_only,
  Fully_corrected     = R2_pathway_fully_corrected,
  Delta_coupling      = R2_pathway_leakfree_only - R2_pathway_fully_corrected
)
print(results_table, row.names = FALSE)

# Ranking pathway per varianza sull'intero dataset (per capire quali entrano fra 30° e 40°)
vars_full <- apply(x_full, 2, var)
ranked_paths <- names(sort(vars_full, decreasing = TRUE))

paths_31_to_40 <- ranked_paths[31:40]
cat("Pathway aggiunti passando da top-30 a top-40:\n")
print(paths_31_to_40)
# print(paths_31_to_40)
#[1] "HALLMARK_MYOGENESIS"                 "HALLMARK_IL2_STAT5_SIGNALING"        "HALLMARK_P53_PATHWAY"                "HALLMARK_COMPLEMENT"                
#[5] "HALLMARK_APOPTOSIS"                  "HALLMARK_DNA_REPAIR"                 "HALLMARK_BILE_ACID_METABOLISM"       "HALLMARK_SPERMATOGENESIS"           
#[9] "HALLMARK_PI3K_AKT_MTOR_SIGNALING"    "HALLMARK_WNT_BETA_CATENIN_SIGNALING"

############


pred <- fit_reg_pathway$pred
# Keep only best tuning parameters
best <- fit_reg_pathway$bestTune
pred_best <- pred[pred$alpha == best$alpha & pred$lambda == best$lambda, ]
# Global R² (across ALL held-out predictions)
R2_global <- cor(pred_best$obs, pred_best$pred)^2
R2_global # 0.46
coefs <- coef(fit_reg_pathway$finalModel, s = fit_reg_pathway$bestTune$lambda)
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df$Feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "Coefficient"
coefs_df <- coefs_df[coefs_df$Coefficient != 0, ]
top30 <- coefs_df[order(abs(coefs_df$Coefficient), decreasing = TRUE), ]
top30 <- top30[top30$Feature != "(Intercept)", ]
top30 <- top30[abs(top30$Coefficient) > 0.0001, ]
ggplot(top30, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "TRUE" = "firebrick",
      "FALSE" = "steelblue"
    ),
    guide = "none"
  ) +
  labs(x = "", y = "Beta Coefficient") +
  ggtitle("50 GSVA pathways as predictors") +
  theme_minimal()
pred <- predict(fit_reg_pathway, t(gsva_nt_all))
round(cor(y_cont, pred),3) # 0.78
# This introduces mathematical coupling Because: GSVA_NT appears inside the definition of y_cont And also as a predictor. 
# predicting a function using one of its own components. This can inflate performance slightly.
# my model is: Baseline (NT) → suppression_score. Where: suppression_score = GSVA_NT − GSVA_CET. This means: Each y value requires a matched NT + CET pair
# my model is not: “predicting a phenotype from baseline”. It is predicting a treatment-induced change. that requires: paired design within-model comparison 


############################################################
#PCA First Often stronger signal because: Removes noise, Avoids p >> n, Stabilizes regression
pca <- prcomp(x)
x_pca <- pca$x[, 1:20]  # first 20 PCs # anche questa usa tutti i campioni. Le componenti principali vedono i dati di test.
folds <- group_folds(groups, k = 3, repeats = 20)
ctrl <- trainControl(method = "cv",index = folds,savePredictions = "final",)
stopifnot(all(names(y_cont) == rownames(x_pca)))
fit_reg_pca <- caret::train(x = x_pca,y = y_cont,method = "glmnet",trControl = ctrl,tuneLength = 20,
                            tuneGrid = expand.grid(alpha = seq(0, 0.1, length = 21),lambda = seq(0.001, 0.1, length = 20)))
sd(y_cont) #  0.1802557
min(fit_reg_pca$results$RMSE) #  0.1354931
min(fit_reg_pca$results$RMSE) < sd(y_cont) # If RMSE < SD → model captures signal.
subset(fit_reg_pca$results,fit_reg_pca$results$RMSE == min(fit_reg_pca$results$RMSE))
#    alpha     lambda      RMSE  Rsquared       MAE     RMSESD RsquaredSD      MAESD
#407   0.1 0.03226316 0.1360171 0.4482416 0.1097962 0.01662612   0.112969 0.0134707

m0 <- lm(suppression_score ~ baseline_score) #Questa procedura rimuove l'effetto meccanico dovuto al fatto che NT compare in entrambi i lati dell'equazione.
resid_response_pca <- resid(m0) 
stopifnot(all(names(resid_response_pca) == rownames(x_pca)))
fit_pca_resid<-caret::train(x=x_pca,y=resid_response_pca,method="glmnet",trControl = ctrl,tuneLength = 20,
                            tuneGrid = expand.grid(alpha = seq(0, 0.1, length = 21),lambda = seq(0.001, 0.1, length = 20)))
subset(fit_pca_resid$results,fit_pca_resid$results$RMSE == min(fit_pca_resid$results$RMSE))
#   alpha     lambda      RMSE  Rsquared       MAE     RMSESD RsquaredSD      MAESD
#51  0.01 0.05310526 0.1287657 0.2103124 0.1044869 0.01550952   0.134214 0.01203901

############# LEAKING REMOVAL fom fit_reg_pca
# Per eliminare completamente il leakage dovuto alla selezione dei geni devo fare la feature selection delle PCs dentro ogni fold. 
# Devo spostare la selezione delle Principal Components all'interno del ciclo CV.Rifai la PCA dentro ogni training fold e proiettare il test fold sui PC del training.
# Versione con feature selection INTERNA anche sui modelli pathway-restricted
# Per sondare la robustezza del segnale a diversi K

##############################################################################
# fit_reg_pca — FULLY CORRECTED (leakage-free + coupling-free)
#
# Fix 1: groupKFold ritorna TRAINING indices
# Fix 2: PCA calcolata SOLO sul training fold, test proiettato con center e
#        rotation del training
# Fix 3: tuning di alpha via inner CV
# Fix 4: coupling correction opzionale — residualizzazione dell'outcome
#        rispetto al baseline_score, con m0 stimato SOLO sul training
##############################################################################

evaluate_pca_leakfree <- function(x_full, y, folds,
                                  baseline_score  = NULL,
                                  n_pc            = 20,
                                  alpha_grid      = seq(0, 0.2, length = 11),
                                  coupling_correct = FALSE) {
  
  # Sanity check
  if (coupling_correct) {
    if (is.null(baseline_score)) {
      stop("coupling_correct = TRUE richiede baseline_score in input.")
    }
    stopifnot(length(baseline_score) == nrow(x_full))
    stopifnot(!is.null(names(baseline_score)))
    stopifnot(identical(names(baseline_score), rownames(x_full)))
  }
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    # FIX 1: TRAINING indices
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(x_full)), train_idx)
    
    x_train <- x_full[train_idx, ]
    x_test  <- x_full[test_idx, ]
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # ----- PCA SOLO sul training fold (leakage removal) -----
    pca_train <- prcomp(x_train, center = TRUE, scale. = FALSE)
    
    # Difesa: limita n_pc a quanti PC sono effettivamente disponibili
    n_pc_eff <- min(n_pc, ncol(pca_train$x))
    if (n_pc_eff < n_pc) {
      warning(sprintf("Fold %d: solo %d PC disponibili (richiesti %d). Uso %d.",
                      i, n_pc_eff, n_pc, n_pc_eff))
    }
    
    x_train_pc <- pca_train$x[, 1:n_pc_eff, drop = FALSE]
    
    # Proiezione del test sui PC del training
    x_test_centered <- scale(x_test, center = pca_train$center, scale = FALSE)
    x_test_pc       <- x_test_centered %*% pca_train$rotation[, 1:n_pc_eff, drop = FALSE]
    
    # ----- Coupling correction opzionale (m0 stimato sul training) -----
    if (coupling_correct) {
      b_train <- baseline_score[train_idx]
      b_test  <- baseline_score[test_idx]
      
      m0 <- lm(y_train ~ b_train)
      y_train_corr <- resid(m0)
      
      # Applica la stessa correzione al test: y_test - predizione da baseline
      pred_from_baseline <- coef(m0)[1] + coef(m0)[2] * b_test
      y_test_corr <- y_test - pred_from_baseline
      
      y_train <- y_train_corr
      y_test  <- y_test_corr
    }
    
    # ----- Tuning di alpha via inner CV -----
    set.seed(123 + i)
    inner_foldid <- sample(rep(1:5, length.out = length(y_train)))
    
    cv_results <- lapply(alpha_grid, function(a) {
      cvfit <- cv.glmnet(x_train_pc, y_train, alpha = a, foldid = inner_foldid)
      list(alpha = a, cvm_min = min(cvfit$cvm), cvfit = cvfit)
    })
    best_idx <- which.min(sapply(cv_results, function(r) r$cvm_min))
    best_fit <- cv_results[[best_idx]]$cvfit
    
    pred <- predict(best_fit, newx = x_test_pc, s = "lambda.min")
    all_pred <- c(all_pred, as.numeric(pred))
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Preparazione dati ----------

# La matrice di predittori: t(expr_baseline) — 106 × 22741 geni human
x_baseline_matrix <- t(expr_baseline)

# Assicurati che i rownames siano puliti dal suffisso _NT
rownames(x_baseline_matrix) <- sub("_NT$", "", rownames(x_baseline_matrix))

# baseline_score con nomi puliti e nello stesso ordine di x_baseline_matrix
baseline_score_clean <- baseline_score
names(baseline_score_clean) <- sub("_NT$", "", names(baseline_score_clean))
baseline_score_aligned <- baseline_score_clean[rownames(x_baseline_matrix)]

# Anche y_cont deve avere nomi coerenti
y_cont_aligned <- y_cont[rownames(x_baseline_matrix)]

# Verifiche
stopifnot(identical(names(baseline_score_aligned), rownames(x_baseline_matrix)))
stopifnot(identical(names(y_cont_aligned), rownames(x_baseline_matrix)))
stopifnot(sum(is.na(baseline_score_aligned)) == 0)


# ---------- Run: due varianti ----------

# Variante 1: leakage-free only (già calcolata prima, dovrebbe dare ~0.40)
R2_pca_leakfree_only <- evaluate_pca_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  baseline_score   = NULL,
  n_pc             = 20,
  coupling_correct = FALSE
)
R2_pca_leakfree_only #0.40

# Variante 2: fully corrected (leakage-free + coupling-corrected)
R2_pca_fully_corrected <- evaluate_pca_leakfree(
  x_full           = x_baseline_matrix,
  y                = y_cont_aligned,
  folds            = folds,
  baseline_score   = baseline_score_aligned,
  n_pc             = 20,
  coupling_correct = TRUE
)
R2_pca_fully_corrected #0.21


# ---------- Riepilogo -----------
cat("\nfit_reg_pca — evaluation regimes:\n")
cat(sprintf("  Leakage-free only (already reported): %.4f\n", R2_pca_leakfree_only))
cat(sprintf("  Fully corrected (leakage + coupling): %.4f\n", R2_pca_fully_corrected))



coefs <- coef(fit_reg_pca$finalModel, s = fit_reg_pca$bestTune$lambda)
head(coefs)
coefs_df <- as.data.frame(as.matrix(coefs))
head(coefs_df)
coefs_df$Feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "Coefficient"
coefs_df <- coefs_df[coefs_df$Coefficient != 0, ]
top30 <- coefs_df[order(abs(coefs_df$Coefficient), decreasing = TRUE), ]
top30 <- top30[top30$Feature != "(Intercept)", ]
top30 <- top30[abs(top30$Coefficient) > 0.0001, ]
ggplot(top30, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "TRUE" = "firebrick",
      "FALSE" = "steelblue"
    ),
    guide = "none"
  ) +
  labs(x = "", y = "Beta Coefficient") +
  ggtitle("Principal components predictors") +
  theme_minimal()
# Interpretation: Predictive signal is concentrated in EGFR-related pathways. We must rule out overfitting optimism. Compute simple correlation:
pred <- predict(fit_reg, x_reduced)
cor(pred, y_cont)^2 # 0.91
pred <- predict(fit_reg_genes_in_the_pathways, x_genes_in_pathways)
cor(pred, y_cont)^2 # 0.78   classic overfitting gap
# The model fits training data well but generalizes modestly — expected in p >> n settings.
pred <- predict(fit_reg_pathway, t(gsva_nt_all))
cor(pred, y_cont)^2 # 0.60 IF the gap is small so reduced overfitting.
pred <- predict(fit_reg_pca, x_pca)
cor(pred, y_cont)^2 # 0.67

#################################################
#################################################
#Now to improve this model: replace mean GSVA suppression with a weighted score
#Right now: suppression_score <- colMeans(suppression_matrix) This assumes: KRAS, MTORC1, E2F contribute equally. That’s biologically unlikely.
#Better idea: Suppression Score = w1*KRAS + w2*MTOR + w3*E2F. Where weights are data-driven.
#Adesso i predittori sono x_genes_in_pathways — cioè i 569 geni EGFR-axis — non i PC. 
#La "PCA" nel nome si riferisce all'outcome (y_cont_pca), che è il primo PC della suppression matrix usato come score pesato invece di colMeans. 
#Il modello predice quindi un outcome PCA-derivato, ma resta un regressore gene-livello.

pca_sup <- prcomp(t(suppression_matrix), scale. = TRUE)
# First principal component
suppression_score_pca <- pca_sup$x[,1]
#Sign ambiguity (important but simple) PCA signs are arbitrary. If you rerun PCA, you might get: -0.63, -0.36, -0.68. So to make it more robust:
# Sign check PRIMA di assegnare y_cont_pca
if(cor(suppression_score_pca, rowMeans(t(suppression_matrix))) < 0 ){suppression_score_pca <- -suppression_score_pca}
pca_sup$rotation #weights are data-driven.
#                                 PC1         PC2        PC3
#HALLMARK_E2F_TARGETS       0.6398000  0.39552697  0.6589494
#HALLMARK_KRAS_SIGNALING_UP 0.3426022 -0.91427962  0.2161401
#HALLMARK_MTORC1_SIGNALING  0.6879532  0.08747109 -0.7204645
# Ora assegna y_cont_pca con il segno già corretto
y_cont_pca <- suppression_score_pca[colnames(expr_baseline)]
# Why this works: Finds the dominant shared pattern of suppression. Automatically assigns weights. Avoids overfitting (only 3 variables)
#y_cont_pca <- suppression_score_pca[colnames(expr_baseline)]
folds <- group_folds(groups, k = 3, repeats = 20)
ctrl <- trainControl(method = "cv",index = folds,savePredictions = "final")
stopifnot(all(names(y_cont_pca) == rownames(x_genes_in_pathways)))
dim(x_genes_in_pathways)[2]
# 569
fit_grouped_pca <- caret::train(x = x_genes_in_pathways,y = y_cont_pca, method = "glmnet",trControl = ctrl,tuneLength = 20,
                                tuneGrid = expand.grid(alpha = seq(0, 0.4, length = 21),lambda = seq(0.01, 1, length = 20)))
sd(y_cont_pca) # 1.35
min(fit_grouped_pca$results$RMSE) #  0.99
min(fit_grouped_pca$results$RMSE) < sd(y_cont_pca) # If RMSE < SD → model captures signal.
subset(fit_grouped_pca$results,fit_grouped_pca$results$RMSE == min(fit_grouped_pca$results$RMSE))
#   alpha lambda      RMSE  Rsquared       MAE    RMSESD RsquaredSD      MAESD
#240  0.11      1 0.9902775 0.4887284 0.791832 0.1094606  0.1036493 0.08842566

m0_pca <- lm(suppression_score_pca ~ baseline_score) #Questa procedura rimuove l'effetto meccanico dovuto al fatto che NT compare in entrambi i lati dell'equazione.
resid_response_pca <- resid(m0_pca) 
stopifnot(all(names(resid_response_pca) == rownames( x_genes_in_pathways)))
fit_grouped_pca_resid<-caret::train(x=x_genes_in_pathways,y=resid_response_pca,method="glmnet",trControl = ctrl,tuneLength = 20,
                                    tuneGrid = expand.grid(alpha = seq(0, 0.2, length = 21),lambda = seq(0.01, 1, length = 20)))
subset(fit_grouped_pca_resid$results,fit_grouped_pca_resid$results$RMSE == min(fit_grouped_pca_resid$results$RMSE))
#   alpha lambda      RMSE  Rsquared        MAE     RMSESD RsquaredSD      MAESD
#200  0.09      1 0.9426405 0.2754918 0.7657408 0.09555809  0.1407425 0.07843237
# sanity check: Rsquared ≈ 1 − RMSE^2/SD(y)^2

##############################################################################
# fit_grouped_pca — LEAKAGE-FREE + COUPLING-CORRECTED
#
# Fix 1: groupKFold returns TRAINING indices
# Fix 2: alpha tuning via inner CV (not hardcoded)
# Fix 3: PCA weights for the outcome computed INSIDE each training fold,
#        then applied to test fold (avoids outcome-side leakage)
# Fix 4: sign alignment of PC1 enforced per fold using training-only signal
# Fix 5: (optional) coupling correction via baseline-score residualization,
#        also computed inside the training fold
##############################################################################

evaluate_grouped_pca_leakfree <- function(x_genes,               # 106 × 569
                                          suppression_matrix,    # 3 × 106 (pathway × model)
                                          baseline_score,        # 106 vettore, per coupling
                                          folds,
                                          alpha_grid = seq(0, 0.4, length = 11),
                                          coupling_correct = FALSE) {
  
  # Trasponi la suppression matrix a (model × pathway) come si aspetta prcomp
  supp_mat_t <- t(suppression_matrix)   # 106 × 3
  
  # Sanity check
  stopifnot(all(rownames(x_genes) == rownames(supp_mat_t)))
  stopifnot(all(rownames(x_genes) == names(baseline_score)))
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    # FIX: TRAINING indices
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(x_genes)), train_idx)
    
    # -------- 1. PCA dell'outcome SOLO sul training --------
    supp_train <- supp_mat_t[train_idx, ]
    supp_test  <- supp_mat_t[test_idx, ]
    
    # Salvo center e scale per applicarli identici al test
    pca_train  <- prcomp(supp_train, center = TRUE, scale. = TRUE)
    
    # Sign alignment: PC1 deve correlare positivamente con la MEDIA
    # dei pathway di training (proxy della direzione biologica "sopprime")
    train_mean_supp <- rowMeans(supp_train)
    if (cor(pca_train$x[, 1], train_mean_supp) < 0) {
      pca_train$rotation[, 1] <- -pca_train$rotation[, 1]
      pca_train$x[, 1]        <- -pca_train$x[, 1]
    }
    
    y_train <- pca_train$x[, 1]  # outcome training: PC1 sui suppression del training
    
    # Proietta il test sui PC del training (usando gli stessi center e scale)
    supp_test_scaled <- scale(supp_test,
                              center = pca_train$center,
                              scale  = pca_train$scale)
    y_test_pc1 <- supp_test_scaled %*% pca_train$rotation[, 1]
    y_test     <- as.numeric(y_test_pc1)
    
    # -------- 2. Coupling correction (opzionale) --------
    if (coupling_correct) {
      # Residualizza sia training che test rispetto al baseline_score
      # Il modello di residualizzazione si stima SOLO sul training
      m0     <- lm(y_train ~ baseline_score[train_idx])
      y_train_corr <- resid(m0)
      # Applica il modello al test
      pred_test_from_baseline <- coef(m0)[1] + coef(m0)[2] * baseline_score[test_idx]
      y_test_corr  <- y_test - pred_test_from_baseline
      
      y_train <- y_train_corr
      y_test  <- y_test_corr
    }
    
    # -------- 3. Predittori: 569 geni EGFR-axis (a priori, nessuna FS) --------
    x_train <- x_genes[train_idx, ]
    x_test  <- x_genes[test_idx, ]
    
    # -------- 4. Elastic net con tuning di alpha (inner CV) --------
    set.seed(123 + i)
    inner_foldid <- sample(rep(1:5, length.out = length(y_train)))
    
    cv_results <- lapply(alpha_grid, function(a) {
      cvfit <- cv.glmnet(x_train, y_train, alpha = a, foldid = inner_foldid)
      list(alpha = a, cvm_min = min(cvfit$cvm), cvfit = cvfit)
    })
    best_idx <- which.min(sapply(cv_results, function(r) r$cvm_min))
    best_fit <- cv_results[[best_idx]]$cvfit
    
    pred <- predict(best_fit, newx = x_test, s = "lambda.min")
    
    all_pred <- c(all_pred, as.numeric(pred))
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Sanity: allineamento nomi ----------
# rownames(x_genes_in_pathways) e colnames(suppression_matrix) devono
# corrispondere per riga/modello. Se non è così, riallinea prima.

# Assicurati che baseline_score sia nominato coerentemente
baseline_score_aligned <- baseline_score[rownames(x_genes_in_pathways)]

# La suppression matrix ha modelli sulle colonne → trasponi per allineare
supp_matrix_aligned <- suppression_matrix[, rownames(x_genes_in_pathways)]

stopifnot(identical(colnames(supp_matrix_aligned), rownames(x_genes_in_pathways)))

# Riallinea baseline_score rimuovendo il suffisso _NT dai nomi
baseline_score_clean <- baseline_score
names(baseline_score_clean) <- sub("_NT$", "", names(baseline_score_clean))

# Ora il subsetting per nome funziona
baseline_score_aligned <- baseline_score_clean[rownames(x_genes_in_pathways)]

# Verifica
sum(is.na(baseline_score_aligned))                                       # atteso: 0
identical(names(baseline_score_aligned), rownames(x_genes_in_pathways))  # atteso: TRUE
stopifnot(identical(names(baseline_score_aligned), rownames(x_genes_in_pathways)))

#leakage-free
R2_grouped_pca_leakage-free <- evaluate_grouped_pca_leakfree(
  x_genes            = x_genes_in_pathways,
  suppression_matrix = supp_matrix_aligned,
  baseline_score     = baseline_score_aligned,
  folds              = folds,
  coupling_correct   = FALSE
)
R2_grouped_pca_leakage-free #0.45

#leakage-free + coupling-corrected (stima onesta)
R2_grouped_pca_full_corrected <- evaluate_grouped_pca_leakfree(
  x_genes            = x_genes_in_pathways,
  suppression_matrix = supp_matrix_aligned,
  baseline_score     = baseline_score_aligned,
  folds              = folds,
  coupling_correct   = TRUE
)
R2_grouped_pca_full_corrected #0.20

coefs <- coef(fit_grouped_pca$finalModel, s = fit_grouped_pca$bestTune$lambda)
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df$Feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "Coefficient"
coefs_df <- coefs_df[coefs_df$Coefficient != 0, ]
top30 <- coefs_df[order(abs(coefs_df$Coefficient), decreasing = TRUE), ][1:20, ]
top30 <- top30[top30$Feature != "(Intercept)", ]
ggplot(top30, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "TRUE" = "firebrick",
      "FALSE" = "steelblue"
    ),
    guide = "none"
  ) +
  labs(x = "", y = "Beta Coefficient") +
  ggtitle("Top 20 EGFR-axis gene predictors among\n 569 total predictors (PCA-weighted outcome)") +
  theme_minimal()
pred <- predict(fit_grouped_pca,  x =  t(expr_baseline))
cor(y_cont_pca, pred)^2 # 0.75


#########################################################
########################################### HYBRID MODEL
# Now gene + pathway genes = hybrid model. Combine gene-level (x_reduced) and pathway-level (PCA)
x_genes_in_pathways[1:5,1:5]
#                                   AK2    ANP32E    ASF1A    ASF1B     ATAD2
#CRC0030LMX0B053XXT_NT         11.37128 10.259300 9.156231 10.06935  9.297617
#CRC0031LMX0A050XXT_NT         11.28961 11.082842 9.158245 10.13628 11.243454
#CRC0032LMX0A040XXT_NT         10.96308  9.391762 8.344355 10.19022 10.041282
#CRC0059LMO0A0300900XR01000_NT 11.40280  9.282717 9.428570 10.59960 10.550726
#CRC0059LMO0A0301700XR01000_NT 11.56707  9.755493 9.580645 10.61686 10.012349

rownames(x_genes_in_pathways) <- sub("_NT$", "", rownames(x_genes_in_pathways))
dim(x_genes_in_pathways)
gsva_nt  <- gsva_scores[, nt_cols]
t_gsva_nt_all = t(gsva_nt_all)
rownames(t_gsva_nt_all) <- sub("_NT$", "", rownames(t_gsva_nt_all))
all(rownames(x_genes_in_pathways) == rownames(t_gsva_nt_all))

excluded_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_MTORC1_SIGNALING"
)

# Rimozione: 50 - 3 = 47 pathway candidati
x_full <- t(gsva_nt_all)
x_full <- x_full[, !colnames(x_full) %in% excluded_pathways]
cat("Pathway inclusi nella predizione:", ncol(x_full), "\n")
# Attesi: 47

vars <- apply(x_genes_in_pathways , 2, var)
top_300_genes = names(sort(vars, decreasing = TRUE))[1:300]
top_300_x_genes_in_pathways = x_genes_in_pathways[,top_300_genes]
dim(top_300_x_genes_in_pathways)

X_combined <- cbind(x_full, top_300_x_genes_in_pathways) #Combine gene-level (x_reduced) and pathway-level (PCA)
#X_combined <- X_combined[, !duplicated(colnames(X_combined))]
dim(X_combined)[2] # 347
stopifnot(all(names(y_cont) == rownames(X_combined)))
fit_hybrid <- caret::train(x = X_combined,y = y_cont,method = "glmnet",trControl = ctrl, tuneGrid = expand.grid(alpha =  seq(0, 0.2, length = 11),lambda = seq(0.001, 0.1, length = 20)))
subset(fit_hybrid$results,fit_hybrid$results$RMSE == min(fit_hybrid$results$RMSE))
#    alpha     lambda      RMSE  Rsquared        MAE     RMSESD RsquaredSD      MAESD
#216   0.2 0.07915789 0.1263187 0.5198025 0.09838936 0.01255792  0.1126867 0.01128182

m0 <- lm(suppression_score ~ baseline_score) #Questa procedura rimuove l'effetto meccanico dovuto al fatto che NT compare in entrambi i lati dell'equazione.
resid_response <- resid(m0) 
stopifnot(all(names(resid_response ) == rownames(X_combined)))
fit_hybrid_resid <- caret::train(x = X_combined,y = resid_response ,method = "glmnet",trControl = ctrl, tuneGrid = expand.grid(alpha =  seq(0, 0.2, length = 11),lambda = seq(0.001, 0.1, length = 20)))
subset(fit_hybrid_resid$results,fit_hybrid_resid$results$RMSE == min(fit_hybrid_resid$results$RMSE))
#    alpha     lambda     RMSE  Rsquared        MAE     RMSESD RsquaredSD      MAESD
# 120   0.1    0.1 0.1252063 0.2522464 0.09998193 0.01383137  0.1117236 0.0125718

############# LEAKING REMOVAL
# Per eliminare completamente il leakage dovuto alla selezione dei geni devo fare la feature selection dentro ogni fold. 
# Devo spostare la selezione dei 300 geni + 47 pathways all'interno del ciclo CV.
##############################################################################
# fit_hybrid — LEAKAGE-FREE + (opzionale) COUPLING-CORRECTED
#
# Predittori:
#   - 47 GSVA pathways (Hallmark, esclusi E2F/KRAS/mTORC1 che entrano nell'outcome)
#   - Top-K geni EGFR-axis (569) selezionati per varianza DENTRO ogni fold
#
# Fix 1: groupKFold ritorna TRAINING indices
# Fix 2: selezione top-K geni fatta DENTRO il training fold
# Fix 3: tuning di alpha via inner CV
# Fix 4: coupling correction opzionale (stima lm sul training, applicata al test)
##############################################################################

evaluate_hybrid_leakfree <- function(x_pathways,           # 106 × 47 (pathway esclusi 3 coupled)
                                     x_egfr_genes,         # 106 × 569 (EGFR-axis)
                                     y,                    # suppression score
                                     baseline_score,       # per coupling correction
                                     folds,
                                     top_k_genes = 300,
                                     alpha_grid  = seq(0, 0.2, length = 11),
                                     coupling_correct = FALSE) {
  
  # Sanity check allineamento
  stopifnot(identical(rownames(x_pathways), rownames(x_egfr_genes)))
  stopifnot(identical(rownames(x_pathways), names(y)))
  if (coupling_correct) {
    stopifnot(identical(rownames(x_pathways), names(baseline_score)))
  }
  
  all_pred <- c()
  all_obs  <- c()
  
  for (i in seq_along(folds)) {
    # FIX: TRAINING indices
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(x_pathways)), train_idx)
    
    # ----- Split pathway (nessuna FS: sono già 47 predittori a priori) -----
    x_train_pw <- x_pathways[train_idx, , drop = FALSE]
    x_test_pw  <- x_pathways[test_idx,  , drop = FALSE]
    
    # ----- Split geni EGFR-axis (569) + FS interna al fold -----
    x_train_g <- x_egfr_genes[train_idx, , drop = FALSE]
    x_test_g  <- x_egfr_genes[test_idx,  , drop = FALSE]
    
    # Feature selection SOLO sul training: top-K per varianza fra i 569
    if (top_k_genes < ncol(x_train_g)) {
      vars       <- apply(x_train_g, 2, var)
      top_genes  <- names(sort(vars, decreasing = TRUE))[1:top_k_genes]
      x_train_g  <- x_train_g[, top_genes, drop = FALSE]
      x_test_g   <- x_test_g[,  top_genes, drop = FALSE]
    }
    
    # ----- Combina pathway + geni -----
    x_train <- cbind(x_train_pw, x_train_g)
    x_test  <- cbind(x_test_pw,  x_test_g)
    
    # ----- Outcome -----
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # ----- Coupling correction opzionale (stima solo sul training) -----
    if (coupling_correct) {
      m0 <- lm(y_train ~ baseline_score[train_idx])
      y_train_corr <- resid(m0)
      # Applica il modello al test: y_test - predizione da baseline
      pred_from_baseline <- coef(m0)[1] + coef(m0)[2] * baseline_score[test_idx]
      y_test_corr <- y_test - pred_from_baseline
      
      y_train <- y_train_corr
      y_test  <- y_test_corr
    }
    
    # ----- Tuning di alpha via inner CV -----
    set.seed(123 + i)
    inner_foldid <- sample(rep(1:5, length.out = length(y_train)))
    
    cv_results <- lapply(alpha_grid, function(a) {
      cvfit <- cv.glmnet(x_train, y_train, alpha = a, foldid = inner_foldid)
      list(alpha = a, cvm_min = min(cvfit$cvm), cvfit = cvfit)
    })
    best_idx <- which.min(sapply(cv_results, function(r) r$cvm_min))
    best_fit <- cv_results[[best_idx]]$cvfit
    
    pred <- predict(best_fit, newx = x_test, s = "lambda.min")
    all_pred <- c(all_pred, as.numeric(pred))
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}


# ---------- Preparazione dati ----------

# Pathway: 50 Hallmark → escludi i 3 coupled → 47
excluded_pathways <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_MTORC1_SIGNALING"
)
x_pathways_47 <- t(gsva_nt_all)
rownames(x_pathways_47) <- sub("_NT$", "", rownames(x_pathways_47))
x_pathways_47 <- x_pathways_47[, !colnames(x_pathways_47) %in% excluded_pathways]
cat("Pathway inclusi:", ncol(x_pathways_47), "\n")   # atteso: 47

# Geni EGFR-axis: 569 predittori (nessuna pre-selezione qui, la FS avviene dentro il fold)
# rownames devono essere puliti dal suffisso _NT
if (any(grepl("_NT$", rownames(x_genes_in_pathways)))) {
  rownames(x_genes_in_pathways) <- sub("_NT$", "", rownames(x_genes_in_pathways))
}

# baseline_score con nomi puliti
baseline_score_clean <- baseline_score
names(baseline_score_clean) <- sub("_NT$", "", names(baseline_score_clean))

# Allineamento su ordine comune
common_ids <- rownames(x_pathways_47)
x_pathways_47        <- x_pathways_47[common_ids, ]
x_egfr_axis          <- x_genes_in_pathways[common_ids, ]
y_cont_aligned       <- y_cont[common_ids]
baseline_score_aligned <- baseline_score_clean[common_ids]

# Verifiche
stopifnot(identical(rownames(x_pathways_47), rownames(x_egfr_axis)))
stopifnot(identical(rownames(x_pathways_47), names(y_cont_aligned)))
stopifnot(sum(is.na(baseline_score_aligned)) == 0)

#leakage-free 
R2_hybrid_leakage_free <- evaluate_hybrid_leakfree(
  x_pathways       = x_pathways_47,
  x_egfr_genes     = x_egfr_axis,
  y                = y_cont_aligned,
  baseline_score   = baseline_score_aligned,
  folds            = folds,
  top_k_genes      = 300,
  coupling_correct = FALSE
)
R2_hybrid_leakage_free # 0.48

#leakage-free + coupling-corrected (stima onesta)
R2_hybrid_fully_corrected <- evaluate_hybrid_leakfree(
  x_pathways       = x_pathways_47,
  x_egfr_genes     = x_egfr_axis,
  y                = y_cont_aligned,
  baseline_score   = baseline_score_aligned,
  folds            = folds,
  top_k_genes      = 300,
  coupling_correct = TRUE
)
R2_hybrid_fully_corrected # 0.20


# ---------- Sensitivity: variazione di top_k_genes ----------
# Utile per Supplementary: quanti geni EGFR-axis "bastano" nell'hybrid?
K_grid <- c(100, 200, 300, 400, 569)
R2_hybrid_sensitivity <- sapply(K_grid, function(k) {
  cat("K =", k, "...\n")
  evaluate_hybrid_leakfree(
    x_pathways       = x_pathways_47,
    x_egfr_genes     = x_egfr_axis,
    y                = y_cont_aligned,
    baseline_score   = baseline_score_aligned,
    folds            = folds,
    top_k_genes      = k,
    coupling_correct = FALSE
  )
})
names(R2_hybrid_sensitivity) <- paste0("top_", K_grid)
R2_hybrid_sensitivity 
#  top_100   top_200   top_300   top_400   top_569 
#0.3891443 0.4698640 0.4827910 0.4895248 0.4764793


coefs <- coef(fit_hybrid$finalModel, s = fit_hybrid$bestTune$lambda)
head(coefs)
coefs_df <- as.data.frame(as.matrix(coefs))
head(coefs_df)
coefs_df$Feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "Coefficient"
coefs_df <- coefs_df[coefs_df$Coefficient != 0, ]
top30 <- coefs_df[order(abs(coefs_df$Coefficient), decreasing = TRUE), ][1:40,]
top30 <- top30[top30$Feature != "(Intercept)", ]
ggplot(top30, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "TRUE" = "firebrick",
      "FALSE" = "steelblue"
    ),
    guide = "none"
  ) +
  labs(x = "", y = "Beta Coefficient") +
  ggtitle(paste("Top 40 predictors in Hybrid model ",dim(X_combined)[2]," predictors", sep="")) +
  theme_minimal()
pred <- predict(fit_hybrid,  x =  X_combined)
cor(y_cont, pred)^2 # 0.79


################################### HYBRID MODEL + INTERACTIONS between the 3 pathways!
#ONE global linear function. same coefficients for all samples. same relationship for all tumors. same biology assumed for everyone
#It does NOT learn distinct subpopulations. 3 Pathway activity + interactions
##############################################################################
# fit_hybrid_interactions — LEAKAGE-FREE con due strategie di decoupling
# Strategia : predittori a interazione da pathway NON coupled
##############################################################################
# ---------- Preparazione dati puliti ----------
# Pulisci rownames
rownames(x_genes_in_pathways) <- sub("_NT$", "", rownames(x_genes_in_pathways))
# Estrai i 3 pathway coupled per l'interaction matrix
coupled_paths <- c("HALLMARK_E2F_TARGETS",
                   "HALLMARK_KRAS_SIGNALING_UP",
                   "HALLMARK_MTORC1_SIGNALING")
gsva_nt_named <- gsva_nt_all
colnames(gsva_nt_named) <- sub("_NT$", "", colnames(gsva_nt_named))
# Pathway coupled: usati per costruire le interazioni "classiche"
gsva_nt_coupled <- t(gsva_nt_named[coupled_paths, ])   # 106 × 3
# Pathway NON coupled: candidati per l'Opzione 3
gsva_nt_decoupled <- t(gsva_nt_named[!rownames(gsva_nt_named) %in% coupled_paths, ])  # 106 × 47
# baseline_score pulito
baseline_score_clean <- baseline_score
names(baseline_score_clean) <- sub("_NT$", "", names(baseline_score_clean))
# Allineamento
common_ids <- rownames(x_genes_in_pathways)
gsva_nt_coupled       <- gsva_nt_coupled[common_ids, ]
gsva_nt_decoupled     <- gsva_nt_decoupled[common_ids, ]
baseline_score_aligned <- baseline_score_clean[common_ids]
y_cont_aligned        <- y_cont[common_ids]
stopifnot(identical(rownames(x_genes_in_pathways), rownames(gsva_nt_coupled)))
stopifnot(identical(rownames(x_genes_in_pathways), rownames(gsva_nt_decoupled)))
stopifnot(sum(is.na(baseline_score_aligned)) == 0)

# ---------- Helper: costruzione matrice interazioni fra 3 pathway ----------
build_interaction_matrix <- function(x3) {
  # x3: matrice n × 3 (pathway attività basali)
  # Ritorna: n × 6 (3 principali + 3 interazioni pairwise)
  stopifnot(ncol(x3) == 3)
  df <- as.data.frame(x3)
  X_int <- model.matrix(~ (.)^2, data = df)[, -1]  # rimuove intercetta
  # Naming pulito
  pnames <- colnames(x3)
  colnames(X_int) <- c(pnames,
                       paste0("int_", pnames[1], "_x_", pnames[2]),
                       paste0("int_", pnames[1], "_x_", pnames[3]),
                       paste0("int_", pnames[2], "_x_", pnames[3]))
  X_int
}

selected_decoupled <- c("HALLMARK_G2M_CHECKPOINT","HALLMARK_GLYCOLYSIS","HALLMARK_MYC_TARGETS_V1")
# Verifica che esistano
stopifnot(all(selected_decoupled %in% colnames(gsva_nt_decoupled)))
x_paths_for_int <- gsva_nt_decoupled[, selected_decoupled]
x_pathways_and_int <- build_interaction_matrix(x_paths_for_int)
dim(x_pathways_and_int)

vars <- apply(x_genes_in_pathways , 2, var)
top_300_genes = names(sort(vars, decreasing = TRUE))[1:300]
top_300_x_genes_in_pathways = x_genes_in_pathways[,top_300_genes]
dim(top_300_x_genes_in_pathways)

X_full <- cbind(x_pathways_and_int,top_300_x_genes_in_pathways)
dim(X_full) #[1] 106 306
X_full[1:4,1:10]
stopifnot(all(names(y_cont) == rownames(X_full)))
fit_hybrid_interactions <- train(x = X_full,y = y_cont, method = "glmnet",trControl = ctrl,tuneLength = 20,tuneGrid = expand.grid(alpha =  seq(0, 0.2, length = 21),lambda = seq(0.001, 0.1, length = 20)))
subset(fit_hybrid_interactions$results, fit_hybrid_interactions$results$RMSE == min(fit_hybrid_interactions$results$RMSE))
#    alpha     lambda     RMSE  Rsquared        MAE     RMSESD RsquaredSD      MAESD
#416   0.2 0.07915789 0.125718 0.5247892 0.09807988 0.01230647  0.1138572 0.01101761

m0 <- lm(suppression_score ~ baseline_score) #Questa procedura rimuove l'effetto meccanico dovuto al fatto che NT compare in entrambi i lati dell'equazione.
resid_response <- resid(m0) 
stopifnot(all(names(resid_response ) == rownames(X_full)))
fit_hybrid_interactions_resid <- train(x = X_full,y = resid_response, method = "glmnet",trControl = ctrl,tuneLength = 20,tuneGrid = expand.grid(alpha =  seq(0, 0.2, length = 11),lambda = seq(0.001, 0.1, length = 20)))
subset(fit_hybrid_interactions_resid$results, fit_hybrid_interactions_resid$results$RMSE == min(fit_hybrid_interactions_resid$results$RMSE))
#    alpha lambda     RMSE  Rsquared        MAE     RMSESD RsquaredSD      MAESD
#120   0.1    0.1 0.1251802 0.2506879 0.0991027 0.01452767  0.1143999 0.01297233

# ---------- Funzione riusabile ----------
evaluate_hybrid_int_leakfree <- function(x_genes,             # 106 × 569 (geni EGFR-axis)
                                         x_paths_for_int,     # 106 × 3 (pathway per interazioni)
                                         y,                   # suppression score
                                         baseline_score,      # per coupling correction
                                         folds,
                                         top_k_genes = 300,
                                         alpha_grid  = seq(0, 0.2, length = 11),
                                         coupling_correct = FALSE) {
  
  stopifnot(identical(rownames(x_genes), rownames(x_paths_for_int)))
  stopifnot(identical(rownames(x_genes), names(y)))
  all_pred <- c()
  all_obs  <- c()
  for (i in seq_along(folds)) {
    # FIX: TRAINING indices
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(x_genes)), train_idx)
    # ----- Split geni EGFR-axis + FS interna al fold -----
    x_train_g <- x_genes[train_idx, , drop = FALSE]
    x_test_g  <- x_genes[test_idx,  , drop = FALSE]
    if (top_k_genes < ncol(x_train_g)) {
      vars      <- apply(x_train_g, 2, var)
      top_genes <- names(sort(vars, decreasing = TRUE))[1:top_k_genes]
      x_train_g <- x_train_g[, top_genes, drop = FALSE]
      x_test_g  <- x_test_g[,  top_genes, drop = FALSE]
    }
    # ----- Costruzione interaction matrix (nessuna FS: 6 predittori a priori) -----
    x_train_int <- build_interaction_matrix(x_paths_for_int[train_idx, , drop = FALSE])
    x_test_int  <- build_interaction_matrix(x_paths_for_int[test_idx,  , drop = FALSE])
    
    # ----- Combina -----
    x_train <- cbind(x_train_g, x_train_int)
    dim(x_train) # 70, 306 = 300 + 3 pathways + 3 interactions
    x_test  <- cbind(x_test_g,  x_test_int)
    dim(x_test) #36 306
    
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    # ----- Coupling correction opzionale -----
    if (coupling_correct) {
      m0 <- lm(y_train ~ baseline_score[train_idx])
      y_train_corr <- resid(m0)
      pred_from_baseline <- coef(m0)[1] + coef(m0)[2] * baseline_score[test_idx]
      y_test_corr <- y_test - pred_from_baseline
      y_train <- y_train_corr
      y_test  <- y_test_corr
    }
    
    # ----- Tuning di alpha via inner CV -----
    set.seed(123 + i)
    inner_foldid <- sample(rep(1:5, length.out = length(y_train)))
    
    cv_results <- lapply(alpha_grid, function(a) {
      cvfit <- cv.glmnet(x_train, y_train, alpha = a, foldid = inner_foldid)
      list(alpha = a, cvm_min = min(cvfit$cvm), cvfit = cvfit)
    })
    best_idx <- which.min(sapply(cv_results, function(r) r$cvm_min))
    best_fit <- cv_results[[best_idx]]$cvfit
    
    pred <- predict(best_fit, newx = x_test, s = "lambda.min")
    all_pred <- c(all_pred, as.numeric(pred))
    all_obs  <- c(all_obs,  y_test)
  }
  
  cor(all_pred, all_obs)^2
}

# interazioni fra 3 pathway NON coupled (rigore massimo) ----------
# Sceglie 3 pathway biologicamente candidati: G2M, GLYCOLYSIS, MYC_TARGETS_V1
# (dominanti nel modello bayesiano gerarchico)
selected_decoupled <- c("HALLMARK_G2M_CHECKPOINT",
                        "HALLMARK_GLYCOLYSIS",
                        "HALLMARK_MYC_TARGETS_V1")

# Verifica che esistano
stopifnot(all(selected_decoupled %in% colnames(gsva_nt_decoupled)))

gsva_nt_int_B <- gsva_nt_decoupled[, selected_decoupled]

R2_hybrid_int_B <- evaluate_hybrid_int_leakfree(
  x_genes          = x_genes_in_pathways,
  x_paths_for_int  = gsva_nt_int_B,                # G2M, Glycolysis, MYC
  y                = y_cont_aligned,
  baseline_score   = baseline_score_aligned,
  folds            = folds,
  coupling_correct = FALSE                          # non necessaria: pathway già decoupled
)
R2_hybrid_int_B #0.49

dim(x_genes_in_pathways)
# Per completezza, anche B con residualization
R2_hybrid_int_B_corrected <- evaluate_hybrid_int_leakfree(
  x_genes          = x_genes_in_pathways,
  x_paths_for_int  = gsva_nt_int_B,
  y                = y_cont_aligned,
  baseline_score   = baseline_score_aligned,
  folds            = folds,
  coupling_correct = TRUE
)
R2_hybrid_int_B_corrected # 0.20

# ---------- Riepilogo ----------
results_table <- data.frame(
  Configuration = c(
    "Strat B: G2M/Glycolysis/MYC interactions, no coupling by construction",
    "Strat B: + coupling correction (double-safety)"
  ),
  R2 = c(R2_hybrid_int_B,R2_hybrid_int_B_corrected)
)
print(results_table)

#Configurazione                                            R²    Cosa correggo
#Naïve (caret, no correzioni)                              0.52  Nulla 
#Coupling-corrected only (caret, FS fuori dal fold)        0.25  Solo coupling — resta FS leakage
#Leakage-free only (loop, FS in fold, no residualization)  0.49  Solo FS — resta coupling
#Fully corrected (loop, FS in fold + residualization)      0.20  Tutto

#Quindi il calo 0.25 → 0.20 di 0.05 punti è esattamente il contributo della feature-selection leakage in questo specifico modello 
#(i 300 geni EGFR-axis top-variance selezionati sui 106 sample invece che sui 70 di training).

coefs <- coef(fit_hybrid_interactions$finalModel, s = fit_hybrid_interactions$bestTune$lambda)
coefs_df <- as.data.frame(as.matrix(coefs))
head(coefs_df)
coefs_df$Feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "Coefficient"
coefs_df <- coefs_df[coefs_df$Coefficient != 0, ]
top30 <- coefs_df[order(abs(coefs_df$Coefficient), decreasing = TRUE), ][1:40,]
top30 <- top30[top30$Feature != "(Intercept)", ]
ggplot(top30, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "TRUE" = "firebrick",
      "FALSE" = "steelblue"
    ),
    guide = "none"
  ) +
  labs(x = "", y = "Beta Coefficient") +
  ggtitle(paste("Top 40 predictors of Hybrid Model (3 pathways with interactions)", sep="")) +
  theme_minimal()
pred <- predict(fit_hybrid_interactions,  x =   X_full)
cor(y_cont, pred)^2 #  0.8

# Multivariate regression: This lets the model:Learn different weights per pathway Capture heterogeneity of response
# A model that learns pathway importance (weights) allows different mechanisms of response (heterogeneity)
# Those are two different problems:
#  | Goal                  | Method                           |
#  | --------------------- | -------------------------------- |
#  | Learn weights         | linear / elastic net             |
#  | Capture heterogeneity | nonlinear / mixture / clustering |


#####################################
#      HIRACHY
#####################################
# I have Patient-Derived Xenografts (PDXs) and Patient/Derived Tumoroids (PDXTs)
# LMX: Liver Metastasis Xenografts,
# LMO: Liver Metastasis Organoids,
# LMH: Liver Metastasis Humanized tumoroids. Questo però va escluso poichè non ci sono campioni trattati con cetuximab che sono LMH
# So I must stratify LMO, LMH and LMX. Right now I mix LMO and LMX. This is wrong since they behave differently. 
# If I fit separate models, I will get: lower SD and a clearer signal
# shared signal + structured variation
# learn different weights per pathway and capture heterogeneity
X <- t(gsva_nt_all)   # ~50 pathways (ideal)
df <- as.data.frame(X)
df$y <- y_cont
df$patient <- factor(substr(rownames(df), 1, 7))
# ATTUALE: vettore di 2 elementi applicato a ~99 righe → riciclo silenzioso
model_type = c("LMO","LMX")
# CORRETTO: estrai dal nome del campione
df$model <- factor(case_when(
  grepl("LMO", rownames(df)) ~ "LMO",
  grepl("LMX", rownames(df)) ~ "LMX",
  TRUE ~ NA_character_
))
library(brms)
# Build Hierarchical Model #allow different weights per model (KEY STEP) OR Even better: varying slopes. This allows pathway effects themselves to vary across patients
  formula <- bf(    
    y ~ model +
      HALLMARK_MTORC1_SIGNALING +
      HALLMARK_E2F_TARGETS +
    (1 + HALLMARK_GLYCOLYSIS + HALLMARK_E2F_TARGETS + HALLMARK_G2M_CHECKPOINT + HALLMARK_MTORC1_SIGNALING | patient) 
    #it means samples from the same patient are correlated. Without this term,the model treats them as independent observations,which inflates noise.
  )
  priors <- c(
    prior(normal(0, 0.2), class = "b"),
    prior(normal(0, 0.5), class = "Intercept"),
    prior(exponential(2), class = "sd"),
    prior(exponential(2), class = "sigma")
  )
  # fit_hier <- brm(formula = formula,data = df,prior = priors,chains = 4,iter = 4000,warmup = 2000,cores = 4,seed = 123,control = list(adapt_delta = 0.99,max_treedepth = 15))
  #saveRDS(fit_hier, "models/fit_hier.rds")
  # poi per ricaricare senza rifare il fitting
  fit_hier <- readRDS("models/fit_hier.rds")  
  #Prior predictive check:
  #fit_prior_check <- brm(formula, data = df, prior = priors,sample_prior = "only", chains = 2)
  #pp_check(fit_prior_check)  # i prior sono sensati?
  bayes_R2(fit_hier)   
 #  Estimate  Est.Error      Q2.5     Q97.5
 #0.7803617 0.06358323 0.6324983 0.8772621
  #Convergenza: controlla sempre sistematicamente, non solo le divergenze:
  rhat_vals <- rhat(fit_hier)
  summary(rhat_vals)  # tutti < 1.01?
  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.9997  1.0004  1.0007  1.0008  1.0010  1.0037
  sum(neff_ratio(fit_hier)>0.1)# 382  # tutti > 0.1? Sono >0.1 382 su 383
  length((neff_ratio(fit_hier)>0.1)) #383
  #Plot posterior distributions This would make an excellent figure:
  mcmc_plot(fit_hier, type = "areas")
  conditional_effects(fit_hier)
  fixef(fit_hier)
#  Estimate  Est.Error        Q2.5      Q97.5
#  Intercept                  0.28290616 0.03480310  0.212818724 0.34996462
#  modelLMX                  -0.03349734 0.02984428 -0.092581152 0.02547553
#  HALLMARK_MTORC1_SIGNALING  0.19418945 0.09651796 -0.001021604 0.37773282
#  HALLMARK_E2F_TARGETS       0.30187959 0.05812151  0.187033286 0.41523262

  ####################
  #If I want different weights per pathway and capture heterogeneity. I need interactions:
  formula <- bf(
    y ~ model * (
        HALLMARK_GLYCOLYSIS +
        HALLMARK_G2M_CHECKPOINT +
        HALLMARK_MTORC1_SIGNALING +
        HALLMARK_E2F_TARGETS +
        HALLMARK_XENOBIOTIC_METABOLISM
    ) + (1 | patient)
  )
  
  # The model expands to: y = β0 + β1pathway + β2model + β3(pathway×model)
  priors <- c(
    prior(normal(0, 0.2), class = "b"),
    prior(normal(0, 0.5), class = "Intercept"),
    prior(exponential(2), class = "sd"),
    prior(exponential(2), class = "sigma")
  )
  # fit_hier2 <- brm(formula = formula,data = df,prior = priors,chains = 4,iter = 4000,warmup = 2000,cores = 4,seed = 123,control = list(adapt_delta = 0.99,max_treedepth = 15))
  #saveRDS(fit_hier2, "models/fit_hier2.rds")
  # poi per ricaricare senza rifare il fitting
  fit_hier2 <- readRDS("models/fit_hier2.rds") 
  # fit_prior_check <- brm(formula, data = df, prior = priors,sample_prior = "only", chains = 2)
  # pp_check(fit_prior_check)  # i prior sono sensati?
  bayes_R2(fit_hier2)     
#Warning messages: There were 17 divergent transitions 
#Estimate Est.Error      Q2.5     Q97.5
#R2 0.7110655 0.07784122 0.5310646 0.8276438
  rhat_vals <- rhat(fit_hier2)
  summary(rhat_vals)  # tutti < 1.01?
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.9997  1.0003  1.0005  1.0007  1.0009  1.0040 
  sum(neff_ratio(fit_hier2)>0.1)  # tutti > 0.1? 87
  length((neff_ratio(fit_hier2)>0.1)) #Sono >0.1 87 su 89
fixef(fit_hier2)
#Estimate Est.Error         Q2.5      Q97.5
#Intercept                                0.30970123 0.05335401  0.20655134 0.41439403
#modelLMX                                -0.06005839 0.06338538 -0.18301870 0.06584747
#HALLMARK_GLYCOLYSIS                     -0.00593477 0.13444783 -0.27255749 0.25498960
#HALLMARK_G2M_CHECKPOINT                  0.11472599 0.11614176 -0.11230671 0.34450629
#HALLMARK_MTORC1_SIGNALING                0.18675051 0.13149232 -0.07163112 0.44667842
#HALLMARK_E2F_TARGETS                     0.29160527 0.11013394  0.07529455 0.50489953
#HALLMARK_XENOBIOTIC_METABOLISM           0.07408419 0.14770166 -0.21437690 0.35926939
#modelLMX:HALLMARK_GLYCOLYSIS            -0.11907934 0.14984475 -0.41073707 0.17664745
#modelLMX:HALLMARK_G2M_CHECKPOINT        -0.19505700 0.13492012 -0.45865886 0.07447491
#modelLMX:HALLMARK_MTORC1_SIGNALING       0.10290101 0.13786079 -0.16837825 0.37526525
#modelLMX:HALLMARK_E2F_TARGETS            0.06133697 0.12400723 -0.18410857 0.30171280
#modelLMX:HALLMARK_XENOBIOTIC_METABOLISM -0.01057795 0.15614558 -0.31725402 0.28955459

#Random slopes performed better. variability is more patient-specific than LMO/LMX-specific.So the dominant heterogeneity axis appears to be patient-level biology rather than merely model class.
#model * pathway assumes only TWO pathway regimes:LMO, LMX. But if heterogeneity is mostly continuous/patient-specific, this is too coarse.
#You can extract pathway heterogeneity directly. This is extremely interesting: data seem to support mostly shared pathway biology + mild patient-specific modulation
# rather than completely distinct patient pathway programs. This is actually plausible for pathway-level GSVA data because:
# pathways are already averaged/compressed signals. much of the extreme gene-level heterogeneity is smoothed away.

#Your model: (1 + glycolysis + mTOR | patient) assumes covariance among random effects. Meaning the model estimates:
#intercept variance
#glycolysis variance
#mTOR variance
#intercept↔glycolysis correlation
#intercept↔mTOR correlation
#glycolysis↔mTOR correlation.
#That is a LOT for 106 samples.You may get a cleaner model using:
formula <- bf(    
  y ~ model + 
  (1 | patient) +
  (0 + HALLMARK_GLYCOLYSIS | patient) +
  (0 + HALLMARK_MTORC1_SIGNALING | patient)
)
priors <- c(
  prior(normal(0, 0.2), class = "b"),
  prior(normal(0, 0.5), class = "Intercept"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)
# fit_hier3 <- brm(formula = formula,data = df,prior = priors,chains = 4,iter = 4000,warmup = 2000,cores = 4,seed = 123,control = list(adapt_delta = 0.99,max_treedepth = 15))
#saveRDS(fit_hier3, "models/fit_hier3.rds")
# poi per ricaricare senza rifare il fitting
fit_hier3 <- readRDS("models/fit_hier3.rds") 
#fit_prior_check <- brm(formula, data = df, prior = priors,sample_prior = "only", chains = 2)
#pp_check(fit_prior_check)  # i prior sono sensati?
bayes_R2(fit_hier3)    
#Estimate Est.Error      Q2.5     Q97.5
#R2 0.4432558 0.1121854 0.2016981 0.6350039
rhat_vals <- rhat(fit_hier3)
summary(rhat_vals)  # tutti < 1.01?
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9997  1.0004  1.0010  1.0012  1.0019  1.0046 
sum(neff_ratio(fit_hier3)>0.1)  # tutti > 0.1? #225
length(neff_ratio(fit_hier3)) #225
fixef(fit_hier3)
#Estimate  Est.Error       Q2.5       Q97.5
#Intercept  0.2672114 0.03338238  0.2023442  0.33267156
#modelLMX  -0.1482165 0.03458087 -0.2159891 -0.07813228

# or even. The || removes covariance estimation. This is often much more stable in omics.
formula <- bf(    
    y ~ model + 
    (1 + HALLMARK_GLYCOLYSIS + HALLMARK_MTORC1_SIGNALING || patient)
)
priors <- c(
  prior(normal(0, 0.2), class = "b"),
  prior(normal(0, 0.5), class = "Intercept"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)
# fit_hier4 <- brm(formula = formula,data = df,prior = priors,chains = 4,iter = 4000,warmup = 2000,cores = 4,seed = 123,control = list(adapt_delta = 0.99,max_treedepth = 15))
#saveRDS(fit_hier4, "models/fit_hier4.rds")
# poi per ricaricare senza rifare il fitting
fit_hier4 <- readRDS("models/fit_hier4.rds") 
#fit_prior_check <- brm(formula, data = df, prior = priors,sample_prior = "only", chains = 2)
#pp_check(fit_prior_check)  # i prior sono sensati?

bayes_R2(fit_hier4)    
#Estimate Est.Error      Q2.5     Q97.5
#R2  0.43658 0.1099966 0.206003 0.6322849
rhat_vals <- rhat(fit_hier4)
summary(rhat_vals)  # tutti < 1.01?
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9997  1.0003  1.0006  1.0007  1.0010  1.0043 
sum(neff_ratio(fit_hier4)>0.1)  # tutti > 0.1? #225
length(neff_ratio(fit_hier4)) #225
fixef(fit_hier4)
#What covariance means biologically: model estimates correlations among: baseline response, glycolysis sensitivity and mTOR sensitivity.
#apparently the data benefit from allowing latent covariance among patient effects.This is common in biological systems because:
#pathways co-vary, metabolic programs co-regulate, tumor states are coordinated.
#the correlated random-slope structure explains real biological variation
# Meaning patient-specific pathway modulation is roughly the same scale as residual noise.
# drug response = shared pathway effects + patient-specific pathway sensitivity + correlated latent tumor states
# yij = αj + β1j x1ij + β2j x2ij + ϵij
# This is substantially more sophisticated than:
# elastic net
# fixed interactions
# subtype-only stratification.
# LMO vs LMX alone was insufficient because tumors are not separating into two clean mechanistic regimes.
# Instead each tumor occupies a continuous latent pathway-response space.
# This is the clearest pathway effect because the credible interval barely excludes zero.
# E2F is consistent with that story. E2F estimate ≈ 0.285.
# Then cluster patients by: glycolysis slope, mTOR slope, intercept.
# You may uncover: metabolic responder states, resistant phenotypes, latent subtypes beyond CMS/LMO/LMX.

# I can inspect which patients are glycolysis-sensitive and which are mTOR-sensitive. For example:
library(tidybayes)
print(n = 200, spread_draws(fit_hier, r_patient[patient,term]))
spread_draws(fit_hier2, r_patient[patient,term])
spread_draws(fit_hier3, r_patient[patient,term])
spread_draws(fit_hier4, r_patient[patient,term])
formula <- bf(    
  y ~ model +
    HALLMARK_G2M_CHECKPOINT +
    HALLMARK_XENOBIOTIC_METABOLISM +
    (1 + HALLMARK_GLYCOLYSIS + HALLMARK_MTORC1_SIGNALING | patient)
)
priors <- c(
  prior(normal(0, 0.2), class = "b"),
  prior(normal(0, 0.5), class = "Intercept"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)
# fit_hier5 <- brm(formula = formula,data = df,prior = priors,chains = 4,iter = 4000,warmup = 2000,cores = 4,seed = 123,control = list(adapt_delta = 0.99,max_treedepth = 15))
#saveRDS(fit_hier5, "models/fit_hier5.rds")
# poi per ricaricare senza rifare il fitting
fit_hier5 <- readRDS("models/fit_hier5.rds") 
#fit_prior_check <- brm(formula, data = df, prior = priors,sample_prior = "only", chains = 2)
#pp_check(fit_prior_check)  # i prior sono sensati?

bayes_R2(fit_hier5)  
#Estimate Est.Error     Q2.5     Q97.5
#R2 0.7479117 0.08422755 0.5429196 0.8699558
rhat_vals <- rhat(fit_hier5)
summary(rhat_vals)  # tutti < 1.01?
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9998  1.0009  1.0014  1.0017  1.0022  1.0096 
sum(neff_ratio(fit_hier5)>0.1)  # tutti > 0.1? 224
length((neff_ratio(fit_hier5)>0.1))#230
fixef(fit_hier5)
#Estimate   Est.Error         Q2.5      Q97.5
#Intercept                       0.26098248 0.05520478  0.15477272 0.36900219
#modelLMX                       -0.02081486 0.03320207 -0.08632227 0.04367614
#HALLMARK_G2M_CHECKPOINT         0.32542005 0.05986721  0.20965605 0.44518043
#HALLMARK_XENOBIOTIC_METABOLISM  0.11160661 0.13255723 -0.14741348 0.37446063
#Plot posterior distributions This would make an excellent figure:
mcmc_plot(fit_hier5, type = "areas")
conditional_effects(fit_hier5)
#  | Pathway    | Old estimate | New estimate |
#  | ---------- | ------------ | ------------ |
#  | G2M        | 0.03         | 0.32         |
#  | Xenobiotic | 0.06         | 0.11         |
# That suggests the removed pathways were mostly adding noise/redundancy. rather than independent biology.
# Why removing E2F and mTOR improved G2M. Because G2M and E2F are highly correlated. And mTOR overlaps partially with proliferative signaling.
# Previously the model had to distribute shared variance among correlated predictors. 
# After simplification G2M absorbs the shared proliferation signal posterior uncertainty shrinks. Classic multicollinearity behavior.

# But the divergent transitions matter! 105 divergent transitions. it is NOT acceptable for inference!
# The posterior geometry is problematic. This usually happens because random-effect covariance matrix is hard to estimate
# data are sparse per patient. centered parameterization struggles.
#I probably do NOT have enough repeated observations per patient to estimate:
#intercept variance, glycolysis variance, mTOR variance, all covariance terms, robustly.
#So the posterior develops funnels.

formula <- bf(
  y ~ model +
    HALLMARK_E2F_TARGETS +
    HALLMARK_XENOBIOTIC_METABOLISM +
    (1 + HALLMARK_GLYCOLYSIS +
       HALLMARK_MTORC1_SIGNALING || patient)
)
# keep patient-specific metabolic modulation. avoid unstable covariance estimation
priors <- c(
  prior(normal(0, 0.2), class = "b"),
  prior(normal(0, 0.5), class = "Intercept"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)
#fit_hier6 <- brm(formula = formula,data = df,prior = priors,chains = 4,iter = 4000,warmup = 2000,cores = 4,seed = 123,control = list(adapt_delta = 0.99,max_treedepth = 15))
#saveRDS(fit_hier6, "models/fit_hier6.rds")
# poi per ricaricare senza rifare il fitting
fit_hier6 <- readRDS("models/fit_hier6.rds") 
#fit_prior_check <- brm(formula, data = df, prior = priors,sample_prior = "only", chains = 2)
#pp_check(fit_prior_check)  # i prior sono sensati?

bayes_R2(fit_hier6)  
#Estimate Est.Error      Q2.5     Q97.5
#R2 0.6387224 0.08984031 0.4392364 0.7844205
rhat_vals <- rhat(fit_hier6)
summary(rhat_vals)  # tutti < 1.01?
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9997  1.0005  1.0009  1.0011  1.0015  1.0045 
sum(neff_ratio(fit_hier6)> 0.1)  # tutti > 0.1? #227
length((neff_ratio(fit_hier6)>0.1))#227

fixef(fit_hier6)
#                                Estimate  Est.Error         Q2.5       Q97.5
#Intercept                       0.30705518 0.04506422  0.21937954 0.39306954
#modelLMX                       -0.03574495 0.03306910 -0.10180438 0.03001595
#HALLMARK_G2M_CHECKPOINT         0.36919961 0.05479692  0.25985991 0.47385014
#HALLMARK_XENOBIOTIC_METABOLISM  0.16004500 0.12650731 -0.08899743 0.40510422

# the covariance structure was not carrying meaningful information. the covariance matrix was weakly identified.
# This is probabl my best model so far: reasonable complexity, good predictive fit, interpretable biology.
# I intentionally dropped glycolysis and mTOR as random slopes only. That means their effects are heterogeneous across patients rather than shared globally.
# Conceptually your model now says: response = shared proliferative biology + shared metabolic adaptation + individualized metabolic sensitivity
# y = β_G2M x_G2M + β_XENO x_XENO + u_gly,j x_gly,j + u_mTOR,j x_mTOR,j + ϵ
# I moved through a very sensible progression
# 1) flat elastic net
# 2) subgroup interactions
# 3) hierarchical intercepts
# 4) varying slopes
# 5) covariance diagnostics
# 6)stabilized hierarchical structure
# Why fit_hier6 is statistically preferable?. The || formulation estimates separate variance terms and no covariance matrix. So instead of estimating:
# Σ = (σ12        ρ12σ1σ2
#      ρ12σ1σ2    σ22    )
# I estimate only independent variances. This dramatically stabilizes hierarchical models in small omics datasets.
#  Σ = (σ12  0
#       0    σ22)

formula <- bf(
  y ~ model +
    HALLMARK_E2F_TARGETS +
    HALLMARK_MTORC1_SIGNALING +
    (1 + HALLMARK_GLYCOLYSIS +     
       HALLMARK_G2M_CHECKPOINT +
       HALLMARK_XENOBIOTIC_METABOLISM || patient)
)
# keep patient-specific metabolic modulation. avoid unstable covariance estimation
priors <- c(
  prior(normal(0, 0.2), class = "b"),
  prior(normal(0, 0.5), class = "Intercept"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)
#fit_hier7 <- brm(formula = formula,data = df,prior = priors,chains = 4,iter = 4000,warmup = 2000,cores = 8,seed = 123,control = list(adapt_delta = 0.99,max_treedepth = 15))
#saveRDS(fit_hier7, "models/fit_hier7.rds")
# poi per ricaricare senza rifare il fitting
fit_hier7 <- readRDS("models/fit_hier7.rds") 
#fit_prior_check <- brm(formula, data = df, prior = priors,sample_prior = "only", chains = 2)
#pp_check(fit_prior_check)  # i prior sono sensati?
fixef(fit_hier7)
#Estimate  Est.Error        Q2.5      Q97.5
#Intercept                  0.30001697 0.02994200  0.24092196 0.35713272
#modelLMX                  -0.04070594 0.02979585 -0.10052192 0.01639862
#HALLMARK_E2F_TARGETS       0.30688912 0.05537473  0.19998759 0.41695528
#HALLMARK_MTORC1_SIGNALING  0.22035536 0.09084431  0.04316455 0.40039816
rhat_vals <- rhat(fit_hier7)
summary(rhat_vals)  # tutti < 1.01?
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9997  1.0004  1.0007  1.0008  1.0012  1.0033  
sum(neff_ratio(fit_hier7)> 0.1)  # tutti > 0.1? #300
length((neff_ratio(fit_hier7)>0.1))#300

#  | Model     | Structure                          | Bayes R² |
#  | --------- | ---------------------------        | -------- |
#  | fit_hier  | correlated random slopes           | 0.78     |
#  | fit_hier2 | fixed model × pathway interactions | 0.71     |     
#  | fit_hier3 | separated uncorrelated slopes      | 0.44     |       
#  | fit_hier4 | no covariance                      | 0.43     |              
#  | fit_hier5 | correlated random slopes           | 0.75     |   
#  | fit_hier6 | independent random slopes          | 0.64     | 
#  | fit_hier7 | independent random slopes          | 0.74     | 

loo_hier <- loo(fit_hier)
loo(fit_hier)
#Pareto k diagnostic values:
#  Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     73    68.9%   233     
#(0.7, 1]   (bad)      32    30.2%   <NA>    
#(1, Inf)   (very bad)  1     0.9%   <NA>    
#Warning message:Found 33 observations with a pareto_k > 0.7 in model 'fit_hier'.
bayes_R2(fit_hier) 
#Estimate  Est.Error      Q2.5     Q97.5
#R2 0.7803617 0.06358323 0.6324983 0.8772621

loo_hier2 <- loo(fit_hier2)
loo(fit_hier2)
#Pareto k diagnostic values:
#  Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     87    82.1%   112     
#(0.7, 1]   (bad)      18    17.0%   <NA>    
#(1, Inf)   (very bad)  1     0.9%   <NA>    
#Warning message:Found 19 observations with a pareto_k > 0.7 in model 'fit_hier2'
bayes_R2(fit_hier2) 
#Estimate  Est.Error      Q2.5     Q97.5
#R2 0.7110655 0.07784122 0.5310646 0.8276438

loo_hier3 <- loo(fit_hier3)
loo(fit_hier3)
#Pareto k diagnostic values:
#  Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     103   97.2%   189     
#(0.7, 1]   (bad)        3    2.8%   <NA>    
#(1, Inf)   (very bad)   0    0.0%   <NA>    
#Warning message:Found 3 observations with a pareto_k > 0.7 in model 'fit_hier3'
bayes_R2(fit_hier3) 
#Estimate Est.Error      Q2.5     Q97.5
#R2 0.4432558 0.1121854 0.2016981 0.6350039

loo_hier4 <- loo(fit_hier4)
loo(fit_hier4)
#Pareto k diagnostic values:
#  Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     102   96.2%   110     
#(0.7, 1]   (bad)        4    3.8%   <NA>    
#(1, Inf)   (very bad)   0    0.0%   <NA>    
#Warning message:Found 4 observations with a pareto_k > 0.7 in model 'fit_hier4'.
bayes_R2(fit_hier4) 
#Estimate Est.Error     Q2.5     Q97.5
#R2  0.43658 0.1099966 0.206003 0.6322849

loo_hier5 <- loo(fit_hier5)
loo(fit_hier5)
#Pareto k diagnostic values:
#  Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     75    70.8%   65      
#(0.7, 1]   (bad)      30    28.3%   <NA>    
#(1, Inf)   (very bad)  1     0.9%   <NA>    
#Warning message:Found 31 observations with a pareto_k > 0.7 in model 'fit_hier5'.
bayes_R2(fit_hier5) 
#Estimate  Est.Error      Q2.5     Q97.5
#R2 0.7479117 0.08422755 0.5429196 0.8699558

loo_hier6 <- loo(fit_hier6)
loo(fit_hier6)
#Pareto k diagnostic values:
#Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     97    91.5%   75      
#(0.7, 1]   (bad)       9     8.5%   <NA>    
#(1, Inf)   (very bad)  0     0.0%   <NA>    
#Warning message:Found 9 observations with a pareto_k > 0.7 in model 'fit_hier6'.
bayes_R2(fit_hier6) 
#Estimate  Est.Error      Q2.5     Q97.5
#R2 0.6387224 0.08984031 0.4392364 0.7844205

loo_hier7 <- loo(fit_hier7)
loo(fit_hier7)
#Pareto k diagnostic values:
#  Count Pct.    Min. ESS
#(-Inf, 0.7]   (good)     89    84.0%   109     
#(0.7, 1]   (bad)      16    15.1%   <NA>    
#(1, Inf)   (very bad)  1     0.9%   <NA>    
#Warning message:Found 17 observations with a pareto_k > 0.7 in model 'fit_hier7'
bayes_R2(fit_hier7)  
#Estimate Est.Error      Q2.5     Q97.5
#R2 0.741212 0.07086819 0.5761253 0.8491388

loo_compare(loo_hier, loo_hier2, loo_hier3, loo_hier4, loo_hier5, loo_hier6, loo_hier7)
#  elpd_diff se_diff
#fit_hier    0.0       0.0  
#fit_hier7  -2.2       2.4  
#fit_hier2  -3.9       4.6  
#fit_hier5  -4.1       2.6  
#fit_hier6 -11.6       4.8  
#fit_hier3 -31.9       9.5  
#fit_hier4 -32.1       9.6  

#In hierarchical Bayesian models, especially with only 106 observations, predictive reliability (LOO diagnostics) is often more important than squeezing out a few extra points of R².
#  | Model     | R²        | Pareto k > 0.7 | Very bad k > 1 | Δelpd |
#  | --------- | --------- | -------------- | -------------- | ----- |
#  | fit_hier  | 0.780     | 33             | 1              | 0.0   |
#  | fit_hier7 | 0.741     | 17             | 1              | -2.2  |
#  | fit_hier2 | 0.711     | 19             | 1              | -3.9  |
#  | fit_hier5 | 0.748     | 31             | 1              | -4.1  |
#  | fit_hier6 | 0.639     | 9              | 0              | -11.6 |
#  | fit_hier3 | 0.443     | 3              | 0              | -31.9 |
#  | fit_hier4 | 0.437     | 4              | 0              | -32.1 |
# Pareto-k: Rules of thumb:
# <5 problematic observations → excellent (fit_hier3 and fit_hier4)
# 5–10 → acceptable (fit_hier6 In my opinion one of the best model)
# 10–20 → caution  (fit_hier2 and  fit_hier7)
# 20 → model may be overfitting (fit_hier to discard)
#  | Model     | Δelpd      |
#  | --------- | ---------- |
#  | fit_hier  | 0          |
#  | fit_hier7 | -2.2 ± 2.4 |
#  | fit_hier2 | -3.9 ± 4.6 |
#  | fit_hier5 | -4.1 ± 2.6 |
# All four are statistically indistinguishable because ∣Δelpd| < se_diff
# I would present fit_hier7 as the primary model.Reason: R2 = 0.74 (close to the maximum 0.78)
# Δelpd = -2.2 ± 2.4 (not meaningfully worse). 17 problematic observations instead of 33
# while fit_hier6 has Excellent diagnostics. No k > 1. Only 9 problematic observations, But R² drops substantially (0.64).A conservative choice, but I lose explanatory power.

# Among the tested hierarchical Bayesian models, Model 7 provided the best balance between explanatory power and predictive robustness. 
# Although the full model achieved a slightly higher Bayesian R² (0.78 vs. 0.74), leave-one-out cross-validation indicated no meaningful 
# difference in predictive performance (Δelpd = -2.2 ± 2.4), while Model 7 reduced the number of influential observations identified 
# by Pareto-k diagnostics from 33 to 17. Consequently, Model 7 was selected as the final model for biological interpretation.   
#y = β0 + β1LMX + β2E2F + β3MTORC1 + u0j + u1jGlycolysis + u2jG2M + u3jXenobiotic, where j indexes patients.
# E2F and MTORC1 have a population-wide effect.
# Glycolysis, G2M checkpoint, and Xenobiotic metabolism are allowed to vary across patients.
# Random-effect covariances are removed (||), making the model much more stable.
# Shared biology: These pathways have a consistent effect across patients. E2F targets, MTORC1 signaling
# Patient-specific biology: The effects of: Glycolysis, G2M checkpoint, Xenobiotic metabolism are heterogeneous.
# In some patients they matter a lot, in others very little. This is exactly what random slopes are intended to capture.
# In fit_hier All pathways were treated as fixed effects. Consequently, the model tried to estimate a single coefficient for Glycolysis, a single coefficient for G2M
# and a single coefficient for Xenobiotic even though these pathways vary between patients. 
# That increases apparent fit R2=0.78 (probably overfitting) but also increases influential observations 33 Pareto-k warnings
# Maximum Rhat: 1.0033. Rule: if <1.01 excellent, if <1.05 acceptable
# From a statistical perspective, fit_hier7 is the model I would defend in a manuscript because it achieves nearly maximal 
# predictive performance while substantially reducing overfitting and improving diagnostic stability.
VarCorr(fit_hier)
VarCorr(fit_hier2)
VarCorr(fit_hier3)
VarCorr(fit_hier4)
VarCorr(fit_hier5)
VarCorr(fit_hier6)
VarCorr(fit_hier7)

# Find a trade-off score. I would define: Score = R2 − 0.005 × (Pareto bad)

#  | Model     |     R² | Pareto bad |      Score |
#  | --------- | -----: | ---------: | ---------: |
#  | fit_hier  | 0.7804 |         33 | **0.6154** |
#  | fit_hier2 | 0.7111 |         19 | **0.6161** |
#  | fit_hier3 | 0.4433 |          3 | **0.4283** |
#  | fit_hier4 | 0.4366 |          4 | **0.4166** |
#  | fit_hier5 | 0.7479 |         31 | **0.5929** |
#  | fit_hier6 | 0.6387 |          9 | **0.5937** |
#  | fit_hier7 | 0.7412 |         17 | **0.6562** |
  
# I  define: Score = R2 − 0.01 × (Pareto bad)
#  | Model     |     Score |
#  | --------- | --------: |
#  | fit_hier7 |     0.571 |
#  | fit_hier2 |     0.521 |
#  | fit_hier6 |     0.549 |
#  | fit_hier  |     0.450 |
#  | fit_hier5 |     0.438 |

#La forza del lavoro non è nell'aver costruito il miglior predittore possibile, ma nell'aver trasformato dati RNA-seq paired pre/post 
#trattamento, in una misura quantitativa e meccanicistica di dipendenza da EGFR. 
#Questo è esattamente il tipo di contributo metodologico che può interessare una rivista di computational oncology o translational bioinformatics.

# =============================================================================
# PTPM — Validazione contro endpoint farmacologici XENTURION (Supplementary Data 6)
#
# Correlazione tra suppression score PTPM e risposta a cetuximab.
# 
# VERSIONE POST-REVISIONE con quattro correzioni:
#   1. Validazione fatta su DUE score: naïve (composite mean) e 
#      coupling-corrected (residualized against baseline pathway score).
#   2. Controllo di specificità contro variazione volume PDX in placebo
#      (deve essere non-significativo).
#   3. FDR (Benjamini-Hochberg) applicato separatamente a Pearson e Spearman.
#   4. Sanity check sul formato ID modelli prima del merge.
#
# Output principali:
#   - Table 6: Pearson correlations (naïve score)
#   - Table 6b: Pearson correlations (coupling-corrected score) [NEW]
#   - Table 7: Spearman correlations (naïve score)
#   - Table 7b: Spearman correlations (coupling-corrected score) [NEW]
#   - Specificity check: placebo control [NEW]
#   - Figure 4: 4-panel scatter (naïve suppression score vs endpoints)
#   - Supplementary Figure S4: 4-panel scatter (coupling-corrected score) [NEW]
# =============================================================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. CARICAMENTO DATI E COSTRUZIONE DEI DUE SCORE
# -----------------------------------------------------------------------------

# Score naïve: composite mean di GSVA(NT) - GSVA(CET) sui 3 pathway EGFR-axis
# Assumo che suppression_score sia già in memoria come vettore named col nome
# del modello (es. "CRC0030LMX0B053XXT")
suppression_score_naive <- suppression_score
stopifnot(length(suppression_score_naive) > 0)
stopifnot(!is.null(names(suppression_score_naive)))

# Score coupling-corrected: residuo di lm(suppression ~ baseline_pathway_score)
# Il baseline_score è la media GSVA(NT) sui 3 pathway EGFR-axis
if (!exists("baseline_score")) {
  # Ricalcola se non è già in memoria
  # baseline_score <- colMeans(gsva_nt) — assumendo gsva_nt disponibile
  stop("baseline_score non disponibile. Ricalcola prima di procedere.")
}

# Allineamento nomi (rimuovi il suffisso _NT se presente)
names(baseline_score) <- sub("_NT$", "", names(baseline_score))
names(suppression_score_naive) <- sub("_NT$", "", names(suppression_score_naive))

# Allinea baseline_score al suppression_score
baseline_aligned <- baseline_score[names(suppression_score_naive)]
stopifnot(sum(is.na(baseline_aligned)) == 0)

# Stima il modello di coupling e calcola i residui
m0 <- lm(suppression_score_naive ~ baseline_aligned)
suppression_score_corrected <- resid(m0)
names(suppression_score_corrected) <- names(suppression_score_naive)

cat(sprintf(
  "Coupling model: R² = %.3f (fraction of suppression explained by baseline)\n",
  summary(m0)$r.squared
))
cat(sprintf(
  "Correlation between naïve and corrected scores: r = %.3f\n\n",
  cor(suppression_score_naive, suppression_score_corrected)
))

# -----------------------------------------------------------------------------
# 2. CARICAMENTO ENDPOINT FARMACOLOGICI
# -----------------------------------------------------------------------------

pharma <- read.table(
  "Cetuximab_response_file_Xenturion.txt",
  row.names        = 1,
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

colnames(pharma) <- c(
  "ATP_5000",         # endpoint primario ex vivo
  "ATP_1250",
  "ATP_20000",
  "Area_5000",
  "Area_1250",
  "Area_20000",
  "PDX_vol_cetux",    # variazione volume tumorale (cetuximab)
  "PDX_vol_placebo",  # variazione volume tumorale (placebo) - CONTROLLO SPECIFICITÀ
  "TGI"               # Tumor Growth Inhibition score
)
pharma$model <- rownames(pharma)

# -----------------------------------------------------------------------------
# 3. SANITY CHECK SUL FORMATO ID MODELLI
# -----------------------------------------------------------------------------

cat("=== Sanity check ID modelli ===\n")
cat("Prime 5 rownames di pharma:\n")
print(head(rownames(pharma), 5))
cat("\nPrime 5 nomi di suppression_score:\n")
print(head(names(suppression_score_naive), 5))

# Estrai ID paziente (primi 7 caratteri) dal suppression score
suppression_patient_ids <- substring(names(suppression_score_naive), 1, 7)

# Verifica overlap
overlap <- intersect(rownames(pharma), unique(suppression_patient_ids))
cat(sprintf("\nOverlap paziente-livello: %d IDs\n", length(overlap)))

if (length(overlap) < 10) {
  stop(sprintf(
    "Overlap insufficiente (%d). Probabile mismatch di formato ID fra pharma (%s) e suppression (%s)",
    length(overlap),
    paste(head(rownames(pharma), 2), collapse = ", "),
    paste(head(names(suppression_score_naive), 2), collapse = ", ")
  ))
}

cat(sprintf("Merge valido: proseguo con %d modelli sovrapposti\n\n", length(overlap)))

# -----------------------------------------------------------------------------
# 4. AGGREGAZIONE A LIVELLO PAZIENTE E MERGE
# -----------------------------------------------------------------------------
# Nota: substring(1, 7) collassa repliche di modello (LMX/LMO/LMH) dello stesso 
# paziente in un singolo score medio. Coerente con XENTURION Supplementary 
# Data 6 che riporta la risposta a livello di paziente.

# Costruzione del data.frame merged con ENTRAMBI gli score

build_merged_df <- function(score_vec, score_name) {
  data.frame(
    model = substring(names(score_vec), 1, 7),
    score = as.numeric(score_vec)
  ) %>%
    group_by(model) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    rename(!!score_name := score) %>%
    inner_join(pharma, by = "model")
}

df_naive     <- build_merged_df(suppression_score_naive,     "suppression_score")
df_corrected <- build_merged_df(suppression_score_corrected, "suppression_score")

cat(sprintf("Modelli nel merge finale (naïve):     %d\n", nrow(df_naive)))
cat(sprintf("Modelli nel merge finale (corrected): %d\n", nrow(df_corrected)))

# Subset con TGI (in vivo)
df_naive_pdx     <- df_naive     %>% filter(!is.na(TGI))
df_corrected_pdx <- df_corrected %>% filter(!is.na(TGI))

cat(sprintf("Modelli con TGI in vivo (naïve):     %d\n", nrow(df_naive_pdx)))
cat(sprintf("Modelli con TGI in vivo (corrected): %d\n\n", nrow(df_corrected_pdx)))

# -----------------------------------------------------------------------------
# 5. CORRELAZIONI — FUNZIONE RIUSABILE
# -----------------------------------------------------------------------------

endpoints <- c("ATP_5000", "ATP_1250", "ATP_20000",
               "Area_5000", "TGI", "PDX_vol_cetux")

compute_correlations <- function(df, endpoints, method = "pearson") {
  results <- lapply(endpoints, function(ep) {
    d <- df %>% filter(!is.na(.data[[ep]]))
    
    ct <- cor.test(d$suppression_score, d[[ep]],
                   method = method,
                   exact  = if (method == "spearman") FALSE else NULL,
                   use    = "complete.obs")
    
    row <- data.frame(
      Endpoint = ep,
      n        = nrow(d),
      estimate = round(as.numeric(ct$estimate), 3),
      p_value  = signif(ct$p.value, 3)
    )
    
    # CI solo per Pearson
    if (method == "pearson" && !is.null(ct$conf.int)) {
      row$CI_low  <- round(ct$conf.int[1], 3)
      row$CI_high <- round(ct$conf.int[2], 3)
    }
    row
  })
  
  out <- do.call(rbind, results)
  
  # FDR Benjamini-Hochberg applicato ai 6 endpoint
  out$p_adj <- signif(p.adjust(out$p_value, method = "BH"), 3)
  out$signif <- with(out, ifelse(p_adj < 0.001, "***",
                                 ifelse(p_adj < 0.01,  "**",
                                        ifelse(p_adj < 0.05,  "*", "ns"))))
  
  # Rinomina estimate in r/rho a seconda del metodo
  colnames(out)[colnames(out) == "estimate"] <- if (method == "pearson") "r" else "rho"
  out
}

# -----------------------------------------------------------------------------
# 6. TABLE 6 e 6b — PEARSON: naïve vs coupling-corrected
# -----------------------------------------------------------------------------

cat("=== Table 6 — Pearson (naïve suppression score) ===\n")
corr_naive_pearson <- compute_correlations(df_naive, endpoints, method = "pearson")
print(corr_naive_pearson)

cat("\n=== Table 6b — Pearson (COUPLING-CORRECTED suppression score) ===\n")
corr_corrected_pearson <- compute_correlations(df_corrected, endpoints, method = "pearson")
print(corr_corrected_pearson)

# -----------------------------------------------------------------------------
# 7. TABLE 7 e 7b — SPEARMAN: naïve vs coupling-corrected
# -----------------------------------------------------------------------------

cat("\n=== Table 7 — Spearman (naïve suppression score) ===\n")
corr_naive_spearman <- compute_correlations(df_naive, endpoints, method = "spearman")
print(corr_naive_spearman)

cat("\n=== Table 7b — Spearman (COUPLING-CORRECTED suppression score) ===\n")
corr_corrected_spearman <- compute_correlations(df_corrected, endpoints, method = "spearman")
print(corr_corrected_spearman)

# -----------------------------------------------------------------------------
# 8. CONTROLLO DI SPECIFICITÀ: placebo
# -----------------------------------------------------------------------------
# La correlazione fra suppression score e volume PDX in braccio placebo deve 
# essere NON significativa. Se lo fosse, indicherebbe che il score predice la 
# crescita tumorale in generale, non la specifica risposta al cetuximab.

cat("\n=== Specificity check: correlazione con volume PDX in PLACEBO ===\n")
cat("(atteso: non significativa)\n\n")

placebo_check <- function(df, score_name, label) {
  d <- df %>% filter(!is.na(PDX_vol_placebo))
  
  ct_p <- cor.test(d$suppression_score, d$PDX_vol_placebo, method = "pearson")
  ct_s <- cor.test(d$suppression_score, d$PDX_vol_placebo, method = "spearman", exact = FALSE)
  
  cat(sprintf("[%s]\n", label))
  cat(sprintf("  n = %d\n", nrow(d)))
  cat(sprintf("  Pearson  r = %+.3f, p = %.3f\n", ct_p$estimate, ct_p$p.value))
  cat(sprintf("  Spearman ρ = %+.3f, p = %.3f\n\n", ct_s$estimate, ct_s$p.value))
  
  invisible(list(pearson = ct_p, spearman = ct_s, n = nrow(d)))
}

specificity_naive     <- placebo_check(df_naive,     "suppression_score", "Naïve score vs placebo volume")
specificity_corrected <- placebo_check(df_corrected, "suppression_score", "Corrected score vs placebo volume")

# Confronto positivo/controllo: quanto la correlazione col cetuximab arm è più 
# forte di quella col placebo arm? Se sono simili, la specificità è compromessa.

cat("=== Confronto direzionale cetuximab vs placebo (positivo/negativo control) ===\n")

d <- df_naive %>% filter(!is.na(PDX_vol_cetux), !is.na(PDX_vol_placebo))
ct_cetux   <- cor.test(d$suppression_score, d$PDX_vol_cetux,   method = "pearson")
ct_placebo <- cor.test(d$suppression_score, d$PDX_vol_placebo, method = "pearson")

cat(sprintf("Cetuximab arm:  r = %+.3f, p = %.3f (atteso: negativo significativo)\n",
            ct_cetux$estimate,   ct_cetux$p.value))
cat(sprintf("Placebo arm:    r = %+.3f, p = %.3f (atteso: non significativo)\n\n",
            ct_placebo$estimate, ct_placebo$p.value))

# -----------------------------------------------------------------------------
# 9. SALVATAGGIO TABELLE
# -----------------------------------------------------------------------------

write.table(corr_naive_pearson,      "Table6_pharma_pearson_naive.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(corr_corrected_pearson,  "Table6b_pharma_pearson_corrected.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(corr_naive_spearman,     "Table7_pharma_spearman_naive.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(corr_corrected_spearman, "Table7b_pharma_spearman_corrected.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------------------------------------------------------
# 10. FIGURA 4 — 4 pannelli con score naïve
# -----------------------------------------------------------------------------

theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.title    = element_text(size = 10),
    axis.text     = element_text(size = 9),
    plot.title    = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

scatter_cor <- function(data, x_var, y_var, x_lab, y_lab, subtitle = "") {
  d <- data %>% filter(!is.na(.data[[y_var]]))
  ggplot(d, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(size = 2, alpha = 0.7, color = "#2166AC") +
    geom_smooth(method  = "lm", se = TRUE, color = "#B2182B",
                fill    = "#FDDBC7", linewidth = 0.8) +
    stat_cor(method     = "pearson",
             aes(label  = paste(..r.label.., ..p.label.., sep = "~`,`~")),
             label.x.npc = "left", label.y.npc = "top",
             size       = 3.2, color = "black") +
    labs(subtitle = subtitle, x = x_lab, y = y_lab) +
    theme_pub
}

make_fig4 <- function(df_in, df_pdx_in, main_title) {
  pA <- scatter_cor(df_in, "suppression_score", "ATP_5000",
                    "PTPM suppression score",
                    "ATP ratio (cetuximab / untreated)\n5,000 cells/well",
                    "Ex vivo PDXT response")
  pB <- scatter_cor(df_pdx_in, "suppression_score", "TGI",
                    "PTPM suppression score",
                    "Tumor Growth Inhibition (TGI)",
                    "In vivo PDX response")
  pC <- scatter_cor(df_in %>% filter(!is.na(PDX_vol_cetux)),
                    "suppression_score", "PDX_vol_cetux",
                    "PTPM suppression score",
                    "Tumor volume change (%)\n3 weeks cetuximab",
                    "In vivo PDX tumor volume")
  pD <- scatter_cor(df_in, "suppression_score", "Area_5000",
                    "PTPM suppression score",
                    "Relative total area (day 7 / day 0)\ncetuximab / untreated",
                    "Ex vivo PDXT longitudinal")
  
  (pA | pB) / (pC | pD) +
    plot_annotation(
      title      = main_title,
      tag_levels = "a",
      tag_prefix = "(", tag_suffix = ")",
      theme = theme(plot.title = element_text(size = 12, face = "bold"),
                    plot.tag   = element_text(size = 11, face = "bold"))
    )
}

fig4_naive <- make_fig4(
  df_naive, df_naive_pdx,
  "PTPM suppression score (naïve) correlates with pharmacological cetuximab response"
)

fig4_corrected <- make_fig4(
  df_corrected, df_corrected_pdx,
  "PTPM suppression score (coupling-corrected) correlates with pharmacological cetuximab response"
)
fig4_corrected
# Salvataggi — formato vettoriale PDF + TIFF alta risoluzione
ggsave("Figure4_PTPM_pharma_validation_naive.pdf",
       plot = fig4_naive, width = 170, height = 140, units = "mm",
       device = cairo_pdf)
ggsave("Figure4_PTPM_pharma_validation_naive.tif",
       plot = fig4_naive, width = 170, height = 140, units = "mm",
       dpi = 600, device = "tiff", compression = "lzw")

ggsave("FigureS4_PTPM_pharma_validation_corrected.pdf",
       plot = fig4_corrected, width = 210, height = 140, units = "mm",
       device = cairo_pdf)
ggsave("FigureS4_PTPM_pharma_validation_corrected.tif",
       plot = fig4_corrected, width = 210, height = 140, units = "mm",
       dpi = 600, device = "tiff", compression = "lzw")

cat("Figure salvate:\n")
cat("  Figure4_PTPM_pharma_validation_naive.pdf/tif\n")
cat("  FigureS4_PTPM_pharma_validation_corrected.pdf/tif\n\n")

# -----------------------------------------------------------------------------
# 11. RIEPILOGO SINTETICO PER IL MANOSCRITTO
# -----------------------------------------------------------------------------

cat("========== RIEPILOGO PER IL MANOSCRITTO ==========\n")
cat("Score utilizzati:\n")
cat(sprintf("  - Naïve (composite mean of 3 EGFR-axis pathways)\n"))
cat(sprintf("  - Coupling-corrected (residualized against baseline pathway score)\n"))
cat(sprintf("    (baseline R² = %.3f, spurious component removed)\n\n",
            summary(m0)$r.squared))

for (score_type in c("naive", "corrected")) {
  corr_p <- if (score_type == "naive") corr_naive_pearson else corr_corrected_pearson
  corr_s <- if (score_type == "naive") corr_naive_spearman else corr_corrected_spearman
  
  cat(sprintf("--- %s SCORE ---\n", toupper(score_type)))
  cat(sprintf("Endpoint significativi (p_adj < 0.05, Pearson):     %d / 6\n",
              sum(corr_p$p_adj < 0.05)))
  cat(sprintf("Endpoint significativi (p_adj < 0.05, Spearman):    %d / 6\n",
              sum(corr_s$p_adj < 0.05)))
  cat(sprintf("Range r Pearson:       [%.2f, %.2f]\n",
              min(corr_p$r), max(corr_p$r)))
  cat(sprintf("Range ρ Spearman:      [%.2f, %.2f]\n\n",
              min(corr_s$rho), max(corr_s$rho)))
}

cat("--- SPECIFICITY CHECK ---\n")
cat(sprintf("Cetuximab arm r = %+.3f (p = %.3f) — atteso significativo\n",
            ct_cetux$estimate, ct_cetux$p.value))
cat(sprintf("Placebo arm   r = %+.3f (p = %.3f) — atteso non significativo\n\n",
            ct_placebo$estimate, ct_placebo$p.value))

cat("========== PIPELINE COMPLETATA ==========\n")

