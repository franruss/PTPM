# =============================================================================
# Shared internal helpers (not exported)
# Identification of paired baseline/treated profiles, used by all run_ptpm*()
# functions (sections 2-7). Two input modes:
#
#   * design mode (any dataset): the user supplies a 'design' data frame that
#     states, for each expression column, the paired unit, the condition and
#     (optionally) the patient used to group the cross-validation;
#   * legacy mode (XENTURION, design = NULL): pairs and treatment are parsed
#     from 'metadata$title' exactly as in PTPM <= 0.1.4, so results on the
#     example dataset are unchanged.
#
# Internally, the averaged profiles are always named "<pair>_NT" (baseline)
# and "<pair>_CET" (treated); these names are internal labels, whatever the
# treatment.
# =============================================================================

# Validate and resolve the gene prefix: "H_" in legacy mode, none otherwise.
.ptpm_resolve_prefix <- function(gene_prefix, legacy) {
  if (is.null(gene_prefix)) return(if (legacy) "H_" else "")
  if (!is.character(gene_prefix) || length(gene_prefix) != 1L || is.na(gene_prefix)) {
    stop("'gene_prefix' must be a single character string (\"\" for no prefix).")
  }
  gene_prefix
}

# Rows whose identifier starts with the prefix (all rows if prefix is "").
.ptpm_has_prefix <- function(x, prefix) {
  if (!nzchar(prefix)) return(rep(TRUE, length(x)))
  startsWith(x, prefix)
}

# Identifier without the prefix (unchanged if prefix is "" or absent).
.ptpm_strip_prefix <- function(x, prefix) {
  if (!nzchar(prefix)) return(x)
  has <- startsWith(x, prefix)
  x[has] <- substring(x[has], nchar(prefix) + 1L)
  x
}

# Returns
#   expr_avg     : genes x profiles matrix, columns "<pair>_NT" / "<pair>_CET"
#                  (replicate samples of a pair and condition are averaged)
#   pair_patient : named character vector, pair -> patient (CV grouping)
#   gene_prefix  : resolved gene prefix
#   mode         : "design" or "legacy"
.ptpm_prepare_pairs <- function(expression, metadata = NULL, design = NULL,
                                baseline_label = "baseline",
                                treated_label = "treated",
                                gene_prefix = NULL, k = 3) {
  legacy <- is.null(design)
  gene_prefix <- .ptpm_resolve_prefix(gene_prefix, legacy)
  if (legacy) {
    out <- .ptpm_pairs_legacy(expression, metadata, k)
  } else {
    out <- .ptpm_pairs_design(expression, design, baseline_label, treated_label, k)
  }
  out$gene_prefix <- gene_prefix
  out$mode <- if (legacy) "legacy" else "design"
  out
}

# -----------------------------------------------------------------------------
# Design mode
# -----------------------------------------------------------------------------
.ptpm_pairs_design <- function(expression, design, baseline_label, treated_label, k) {
  if (!is.data.frame(design)) stop("'design' must be a data frame.")
  for (lab in list(baseline_label, treated_label)) {
    if (!is.character(lab) || length(lab) != 1L || is.na(lab) || !nzchar(lab)) {
      stop("'baseline_label' and 'treated_label' must be single non-empty strings.")
    }
  }
  if (identical(baseline_label, treated_label)) {
    stop("'baseline_label' and 'treated_label' must differ.")
  }
  need <- c("sample", "pair", "condition")
  miss <- setdiff(need, names(design))
  if (length(miss)) {
    stop("'design' must contain the columns 'sample', 'pair' and 'condition' ",
         "(and optionally 'patient'); missing: ", paste(miss, collapse = ", "))
  }
  d <- data.frame(
    sample = as.character(design$sample),
    pair = as.character(design$pair),
    condition = as.character(design$condition),
    patient = if ("patient" %in% names(design)) as.character(design$patient)
              else as.character(design$pair),
    stringsAsFactors = FALSE
  )
  if (anyNA(d$sample) || anyNA(d$pair) || any(!nzchar(d$sample)) || any(!nzchar(d$pair))) {
    stop("'design$sample' and 'design$pair' must not contain missing or empty values.")
  }
  if (anyNA(d$patient) || any(!nzchar(d$patient))) {
    stop("'design$patient' must not contain missing or empty values.")
  }
  if (anyDuplicated(d$sample)) {
    stop("Each expression sample must appear only once in 'design$sample'; duplicated: ",
         paste(unique(d$sample[duplicated(d$sample)]), collapse = ", "))
  }
  # keep only the two conditions of interest
  other <- unique(d$condition[!d$condition %in% c(baseline_label, treated_label)])
  if (length(other)) {
    message("PTPM: design rows with condition ", paste(shQuote(other), collapse = ", "),
            " are ignored (only '", baseline_label, "' and '", treated_label, "' are used).")
  }
  d <- d[d$condition %in% c(baseline_label, treated_label), , drop = FALSE]
  if (!nrow(d)) {
    stop("No design rows with condition '", baseline_label, "' or '", treated_label, "'.")
  }
  missing_expr <- setdiff(d$sample, colnames(expression))
  if (length(missing_expr)) {
    stop("Samples listed in 'design' are missing from the expression columns: ",
         paste(utils::head(missing_expr, 10), collapse = ", "),
         if (length(missing_expr) > 10) ", ..." else "")
  }
  # a pair belongs to exactly one patient
  pp <- unique(d[, c("pair", "patient")])
  if (anyDuplicated(pp$pair)) {
    stop("Each pair must belong to a single patient; pairs with more than one patient: ",
         paste(unique(pp$pair[duplicated(pp$pair)]), collapse = ", "))
  }
  # complete pairs only
  has_b <- unique(d$pair[d$condition == baseline_label])
  has_t <- unique(d$pair[d$condition == treated_label])
  complete <- intersect(has_b, has_t)
  dropped <- setdiff(union(has_b, has_t), complete)
  if (length(dropped)) {
    message("PTPM: ", length(dropped), " pair(s) without both a '", baseline_label,
            "' and a '", treated_label, "' profile were excluded: ",
            paste(utils::head(dropped, 10), collapse = ", "),
            if (length(dropped) > 10) ", ..." else "")
  }
  if (length(complete) < k) {
    stop("Only ", length(complete), " complete pairs were found; at least k = ", k,
         " are needed.")
  }
  d <- d[d$pair %in% complete, , drop = FALSE]
  d$internal <- paste(d$pair, ifelse(d$condition == baseline_label, "NT", "CET"), sep = "_")
  internal_levels <- unique(d$internal)
  if (anyDuplicated(internal_levels)) stop("Pair identifiers produce duplicated profile names.")

  n_rep <- table(d$internal)
  if (any(n_rep > 1L)) {
    message("PTPM: ", sum(n_rep > 1L), " pair/condition combination(s) have replicate ",
            "samples; their expression profiles are averaged.")
  }
  expr_avg <- vapply(internal_levels, function(g) {
    rowMeans(expression[, d$sample[d$internal == g], drop = FALSE])
  }, numeric(nrow(expression)))
  if (is.null(dim(expr_avg))) expr_avg <- matrix(expr_avg, ncol = 1L)
  rownames(expr_avg) <- rownames(expression)
  colnames(expr_avg) <- internal_levels

  pair_patient <- stats::setNames(pp$patient, pp$pair)[complete]
  list(expr_avg = expr_avg, pair_patient = pair_patient)
}

# -----------------------------------------------------------------------------
# Legacy mode (XENTURION, GSE204805): identical to sections 2-7 of PTPM 0.1.4
# -----------------------------------------------------------------------------
.ptpm_pairs_legacy <- function(expression, metadata, k) {
  if (is.null(metadata)) {
    stop("Supply either 'design' (any dataset) or 'metadata' with a 'title' column ",
         "(XENTURION example data).")
  }
  if (!is.data.frame(metadata)) stop("'metadata' must be a data frame.")
  if (!"title" %in% names(metadata)) stop("'metadata' must contain a column named 'title'.")

  # 3. TREATMENT PARSING
  parse_treatment <- function(title) {
    if (grepl("cetux_6w_NT|cetuxi_h_NT|cetux_h_NT", title)) {
      "NT"
    } else if (grepl("cetux|cetuxi", title)) {
      "CET"
    } else if (grepl("BASALE", title)) {
      "BASALE"
    } else {
      NA_character_
    }
  }
  metadata$treatment_geo <- vapply(as.character(metadata$title), parse_treatment,
                                   character(1))

  # 4. PAIR IDENTIFICATION
  metadata$entire_codes <- substr(as.character(metadata$title), 1, 26)
  metadata$pair_id <- metadata$entire_codes
  metadata$pair_id <- gsub("(LMO.*[0-9]{3})[0-9]R", "\\1XR", metadata$pair_id)
  metadata$pair_id <- gsub("(LMX.*[0-9]{3})[0-9]{2}T", "\\1XXT", metadata$pair_id)
  metadata$pair_id <- gsub("UMR[0-9]+$", "", metadata$pair_id)

  pair_table <- table(metadata$pair_id, metadata$treatment_geo)
  if (!all(c("CET", "NT") %in% colnames(pair_table))) {
    stop("Metadata do not contain both CET and NT samples after treatment parsing. ",
         "For datasets other than the XENTURION example, describe the samples with ",
         "the 'design' argument.")
  }

  # 5. PAIRED SAMPLES
  paired_ids <- rownames(pair_table)[pair_table[, "CET"] > 0 & pair_table[, "NT"] > 0]
  if (length(paired_ids) < k) {
    stop("Fewer paired models were identified than the requested number of outer folds.")
  }

  # 6. PAIRED PHENOTYPE
  pheno_paired <- metadata[metadata$pair_id %in% paired_ids &
                             metadata$treatment_geo %in% c("NT", "CET"), , drop = FALSE]
  missing_expr <- setdiff(pheno_paired$entire_codes, colnames(expression))
  if (length(missing_expr) > 0L) {
    stop("Metadata-derived sample identifiers are missing from expression columns: ",
         paste(missing_expr, collapse = ", "))
  }
  pheno_paired$group <- paste(pheno_paired$pair_id, pheno_paired$treatment_geo, sep = "_")
  group_levels <- unique(pheno_paired$group)

  # 7. AVERAGE EXPRESSION PER GROUP
  expr_avg <- sapply(group_levels, function(g) {
    cols <- pheno_paired$entire_codes[pheno_paired$group == g]
    rowMeans(expression[, cols, drop = FALSE])
  })
  if (is.null(dim(expr_avg))) expr_avg <- matrix(expr_avg, ncol = 1L)
  rownames(expr_avg) <- rownames(expression)
  colnames(expr_avg) <- group_levels

  # Patient = first 7 characters of the model code (e.g. "CRC0078")
  pair_patient <- stats::setNames(substr(paired_ids, 1, 7), paired_ids)
  list(expr_avg = expr_avg, pair_patient = pair_patient)
}

# -----------------------------------------------------------------------------
# Balanced patient-grouped outer folds: patients are shuffled and dealt to the
# k folds in turn, so fold sizes (in patients) differ by at most one and no
# fold is empty. caret::groupKFold() assigns patients to folds by sampling
# with replacement, which gives unequal and occasionally empty folds in small
# cohorts; it is kept for the XENTURION example (balanced_folds = FALSE) so
# that the published results are reproduced exactly.
# Returns a list of TRAINING indices, as caret::groupKFold().
# -----------------------------------------------------------------------------
.ptpm_balanced_group_folds <- function(groups, k) {
  u <- unique(groups)
  if (k > length(u)) stop("'k' is larger than the number of patients.")
  f <- sample(rep(seq_len(k), length.out = length(u)))
  names(f) <- u
  gf <- f[groups]
  stats::setNames(lapply(seq_len(k), function(j) which(gf != j)),
                  paste0("Fold", seq_len(k)))
}
