#' Prepare a paired GEO microarray series for PTPM
#'
#' Downloads a GEO series with [GEOquery::getGEO()], converts probe-level
#' intensities to one value per gene symbol and builds the `design` table
#' required by [run_ptpm()] and the other `run_ptpm*()` functions when they
#' are applied to datasets other than the XENTURION example.
#'
#' The function is written for series in which each patient (or model) has a
#' baseline sample and a treated sample, as in window-of-opportunity trials.
#' Sample annotations differ between GEO series, so the fields that identify
#' the patient, the time point and, if present, the untreated control arm are
#' detected automatically from the sample characteristics and titles. If the
#' detection fails or is ambiguous, the function stops and lists the
#' available fields; the fields can then be given explicitly.
#'
#' Processing steps:
#' \enumerate{
#'   \item intensities are log2-transformed when they are not already on a
#'     log scale (`log_transform = "auto"`, same rule as GEO2R);
#'   \item probes are annotated with the first gene symbol of the platform
#'     annotation, probes without a symbol are removed, and for each symbol
#'     the probe with the highest mean intensity across all samples is kept
#'     (an unsupervised step that does not use the phenotype);
#'   \item each sample is classified as baseline or treated with
#'     `baseline_pattern` and `treated_pattern`; samples matching neither (or
#'     both) are left out;
#'   \item patients in the control arm (no treatment between the two
#'     biopsies) are returned separately in `design_control`; they can be
#'     analysed as a negative control in which S contains no drug effect.
#' }
#'
#' @param gse GEO series accession, e.g. `"GSE109756"`.
#' @param patient_field,time_field,arm_field Optional names of sample
#'   annotation columns (as listed in the `samples` element of the output,
#'   e.g. `"patient id:ch1"`) that identify the patient, the time point and
#'   the treatment arm. `NULL` = automatic detection.
#' @param baseline_pattern,treated_pattern Regular expressions (case
#'   insensitive) that recognise baseline and treated samples in the time
#'   field.
#' @param control_pattern Regular expression that recognises the untreated
#'   control arm in the arm field.
#' @param platform Platform accession to use when the series has more than
#'   one; `NULL` = the platform with most samples.
#' @param symbol_field Column of the platform annotation with gene symbols;
#'   `NULL` = automatic detection.
#' @param log_transform `"auto"`, `"yes"` or `"no"`.
#' @param k Minimum number of complete pairs required (the number of outer
#'   folds that will be used).
#' @param cache_dir Directory where GEO files are downloaded.
#' @param patient_regex Optional regular expression with one capture group
#'   that extracts the patient identifier from the patient field, e.g.
#'   `"^(\\d+)[A-Z]"` to obtain "10" from the title "10A;pretreatment;...".
#'   If `patient_field` is `NULL`, the title is used.
#' @param annotation_only Logical; if `TRUE`, return the expression matrix
#'   and the sample annotation table (with a first guess of the time point
#'   from the title) without building the design. Useful to explore an
#'   unfamiliar series and to write the `design` table by hand, e.g. when a
#'   patient has several models or several treated time points.
#' @param eset Optional `ExpressionSet` already loaded (the download is then
#'   skipped). Mainly for testing.
#'
#' @return A list with
#'   \describe{
#'     \item{`expression`}{data frame, first column gene symbol, then one
#'       column per sample (GSM accession);}
#'     \item{`design`}{design table of the treated arm (`sample`, `pair`,
#'       `condition` = "baseline"/"treated", `patient`), to pass to
#'       `run_ptpm(design = ...)`;}
#'     \item{`design_control`}{design table of the control arm, or `NULL`;}
#'     \item{`samples`}{sample annotation table with the classification of
#'       every sample (for checking);}
#'     \item{`fields`}{the annotation fields used.}
#'   }
#' @examples
#' \dontrun{
#' geo <- load_ptpm_geo("GSE109756")
#' geo$samples          # check the classification of every sample
#' res <- run_ptpm(expression = geo$expression, design = geo$design,
#'                 top_k_values = c(3000, 300))
#' }
#' @export
load_ptpm_geo <- function(gse,
                          patient_field = NULL,
                          time_field = NULL,
                          arm_field = NULL,
                          baseline_pattern = paste0(
                            "baseline|pre[- ]?(treat|therap|op|cetux)\\w*|before|",
                            "\\bpre\\b|diagnos\\w*|biopsy ?1\\b|\\bT0\\b|day ?0\\b|week ?0\\b"),
                          treated_pattern = paste0(
                            "post[- ]?(treat|therap|cetux)\\w*|after|on[- ]?treat\\w*|surg\\w*|",
                            "operat\\w*|resect\\w*|\\bpost\\b|biopsy ?2\\b|\\bT1\\b|day ?[1-9][0-9]?\\b|",
                            "week ?[1-9]\\b|weeks?\\b|\\b\\d+ ?(-|to) ?\\d+ ?(days?|months?)\\b|",
                            "\\b\\d+ ?(days?|months?)\\b|\\d+(\\.\\d+)? ?mg\\b"),
                          control_pattern = "control|untreated|no treatment|^none$|placebo",
                          platform = NULL,
                          symbol_field = NULL,
                          log_transform = c("auto", "yes", "no"),
                          k = 3,
                          cache_dir = tools::R_user_dir("PTPM", which = "cache"),
                          patient_regex = NULL,
                          annotation_only = FALSE,
                          eset = NULL) {
  log_transform <- match.arg(log_transform)
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Package 'Biobase' is required.")
  if (is.null(eset)) {
    if (!requireNamespace("GEOquery", quietly = TRUE)) {
      stop("Package 'GEOquery' is required: BiocManager::install(\"GEOquery\").")
    }
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    old <- options(timeout = max(3600, getOption("timeout")))
    on.exit(options(old), add = TRUE)
    message("Downloading ", gse, " from GEO (cached in ", cache_dir, ")...")
    gl <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = TRUE, destdir = cache_dir)
    plats <- vapply(gl, Biobase::annotation, character(1))
    if (!is.null(platform)) {
      if (!platform %in% plats) {
        stop("Platform ", platform, " not found in ", gse, "; available: ",
             paste(plats, collapse = ", "))
      }
      eset <- gl[[match(platform, plats)]]
    } else {
      eset <- gl[[which.max(vapply(gl, ncol, numeric(1)))]]
      if (length(gl) > 1L) {
        message(gse, " has several platforms (", paste(plats, collapse = ", "),
                "); using ", Biobase::annotation(eset), ".")
      }
    }
  }
  .ptpm_geo_to_input(eset, patient_field, time_field, arm_field, baseline_pattern,
                     treated_pattern, control_pattern, symbol_field, log_transform, k,
                     patient_regex = patient_regex, annotation_only = annotation_only)
}

#' @rdname load_ptpm_geo
#' @param ... Further arguments passed to [load_ptpm_geo()].
#' @details `load_ptpm_gse109756()` loads GSE109756 (CHIRON window-of-opportunity
#'   study: head and neck squamous cell carcinoma, biopsies at baseline and at
#'   surgery after about two weeks of cetuximab; Affymetrix HG-U133 Plus 2.0).
#'   Always inspect `$samples` to confirm that patients and time points were
#'   recognised correctly.
#' @export
load_ptpm_gse109756 <- function(...) load_ptpm_geo("GSE109756", ...)


# Conversion of an ExpressionSet into PTPM input (internal; testable offline).
.ptpm_geo_to_input <- function(eset, patient_field, time_field, arm_field,
                               baseline_pattern, treated_pattern, control_pattern,
                               symbol_field, log_transform, k,
                               patient_regex = NULL, annotation_only = FALSE) {
  ex <- Biobase::exprs(eset)
  pd <- Biobase::pData(eset)
  fd <- Biobase::fData(eset)
  if (!nrow(ex) || !ncol(ex)) stop("The expression set is empty.")

  # 1. log2 scale (GEO2R rule) ------------------------------------------------
  qx <- as.numeric(stats::quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
  needs_log <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  do_log <- switch(log_transform, auto = needs_log, yes = TRUE, no = FALSE)
  if (do_log) {
    ex[ex <= 0] <- NA
    ex <- log2(ex)
  }
  message("PTPM: intensities ", if (do_log) "log2-transformed" else "used as provided (log scale)", ".")

  # 2. probes -> gene symbols ----------------------------------------------------
  if (is.null(symbol_field)) {
    cand <- grep("^gene[ ._]?symbol$|^symbol$|^gene_assignment$", names(fd),
                 ignore.case = TRUE, value = TRUE)
    if (!length(cand)) {
      stop("No gene-symbol column found in the platform annotation; columns: ",
           paste(names(fd), collapse = ", "), ". Set 'symbol_field'.")
    }
    symbol_field <- cand[1]
  }
  sym <- as.character(fd[[symbol_field]])
  if (grepl("gene_assignment", symbol_field, ignore.case = TRUE)) {
    sym <- vapply(strsplit(sym, " // ", fixed = TRUE),
                  function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
  }
  sym <- trimws(sub(" ?///.*$", "", sym))
  ok <- !is.na(sym) & nzchar(sym) & sym != "---" & rowSums(is.na(ex)) == 0
  ex <- ex[ok, , drop = FALSE]
  sym <- sym[ok]
  mean_int <- rowMeans(ex)
  ord <- order(sym, -mean_int)
  keep <- ord[!duplicated(sym[ord])]
  ex <- ex[keep, , drop = FALSE]
  rownames(ex) <- sym[keep]
  message("PTPM: ", nrow(ex), " genes (highest-mean probe per symbol).")

  # 3. sample annotation --------------------------------------------------------
  ann_cols <- c(intersect(c("title", "source_name_ch1", "description"), names(pd)),
                grep(":ch1$", names(pd), value = TRUE))
  ann <- data.frame(lapply(pd[, ann_cols, drop = FALSE], as.character),
                    check.names = FALSE, stringsAsFactors = FALSE)
  ann[is.na(ann)] <- ""
  rownames(ann) <- colnames(ex)
  classify <- function(v) {
    b <- grepl(baseline_pattern, v, ignore.case = TRUE, perl = TRUE)
    t <- grepl(treated_pattern, v, ignore.case = TRUE, perl = TRUE)
    out <- rep(NA_character_, length(v))
    out[b & !t] <- "baseline"
    out[t & !b] <- "treated"
    out
  }
  check_field <- function(f, what) {
    if (!is.null(f) && !f %in% names(ann)) {
      stop("'", what, "' = '", f, "' is not a sample annotation field; available: ",
           paste(names(ann), collapse = ", "))
    }
  }
  check_field(time_field, "time_field")
  check_field(patient_field, "patient_field")
  check_field(arm_field, "arm_field")

  if (isTRUE(annotation_only)) {
    guess <- if ("title" %in% names(ann)) classify(ann$title) else NA_character_
    expression <- data.frame(gene_id = rownames(ex), ex, check.names = FALSE,
                             stringsAsFactors = FALSE)
    rownames(expression) <- NULL
    return(list(expression = expression,
                samples = data.frame(sample = colnames(ex), condition_from_title = guess, ann,
                                     check.names = FALSE, stringsAsFactors = FALSE)))
  }

  # time point: field with most samples classified and both classes present
  if (is.null(time_field)) {
    sc <- vapply(names(ann), function(f) {
      cl <- classify(ann[[f]])
      if (all(c("baseline", "treated") %in% cl)) sum(!is.na(cl)) + 0.5 * grepl(":ch1$", f) else 0
    }, numeric(1))
    if (max(sc) == 0) {
      .ptpm_geo_fail(ann, "No annotation field separates baseline from treated samples. ",
                     "Set 'time_field' and, if needed, 'baseline_pattern'/'treated_pattern'.")
    }
    time_field <- names(ann)[which.max(sc)]
  }
  condition <- classify(ann[[time_field]])

  # control arm
  if (is.null(arm_field)) {
    cand <- names(ann)[vapply(names(ann), function(f) {
      m <- grepl(control_pattern, ann[[f]], ignore.case = TRUE, perl = TRUE)
      any(m) && !all(m) && f != time_field
    }, logical(1))]
    arm_field <- if (length(cand)) cand[1] else NA_character_
  }
  is_control <- if (is.na(arm_field)) rep(FALSE, nrow(ann)) else
    grepl(control_pattern, ann[[arm_field]], ignore.case = TRUE, perl = TRUE)

  # patient: explicit field, or the annotation field (or the title without
  # the time-point words) that forms the most complete baseline/treated pairs
  # with at most one sample per patient and time point. Ties are broken in
  # favour of fields named like an identifier.
  title_id <- if ("title" %in% names(ann)) {
    z <- gsub(paste0("(", baseline_pattern, ")|(", treated_pattern, ")"), "",
              ann$title, ignore.case = TRUE, perl = TRUE)
    gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", z))
  } else NULL
  if (is.null(patient_field) && !is.null(patient_regex)) {
    if (!"title" %in% names(ann)) stop("'patient_regex' needs 'patient_field' (no title found).")
    patient_field <- "title"
  }
  if (is.null(patient_field)) {
    cands <- setdiff(names(ann), c(time_field, if (!is.na(arm_field)) arm_field))
    ids <- c(stats::setNames(lapply(cands, function(f) ann[[f]]), cands),
             if (!is.null(title_id)) list(`title (time-point words removed)` = title_id),
             .ptpm_title_id_candidates(if ("title" %in% names(ann)) ann$title else NULL))
    score <- vapply(names(ids), function(f) {
      id <- ids[[f]]
      use <- !is.na(condition) & nzchar(id)
      key <- paste(id[use], condition[use], sep = "\r")
      if (!any(use) || anyDuplicated(key)) return(-1)
      n_pairs <- length(intersect(id[use & condition == "baseline"],
                                  id[use & condition == "treated"]))
      n_pairs + 0.5 * grepl("patient|subject|individual|donor|case|participant|\\bid\\b|_id",
                            f, ignore.case = TRUE)
    }, numeric(1))
    if (max(score) < 1) {
      # a valid identifier that fails only because patients have several
      # treated samples: the series has more than one treated time point
      multi <- vapply(ids, function(id) {
        use <- !is.na(condition) & nzchar(id)
        bl <- id[use & condition == "baseline"]
        trt <- id[use & condition == "treated"]
        !anyDuplicated(bl) && anyDuplicated(trt) > 0 &&
          length(intersect(bl, trt)) >= k
      }, logical(1))
      if (any(multi)) {
        tv <- unique(ann[[time_field]][condition %in% "treated"])
        .ptpm_geo_fail(ann, "Patients have more than one treated sample (several treated time ",
                       "points). Choose one with 'treated_pattern'. Treated values found: ",
                       paste(utils::head(tv, 8), collapse = " | "))
      }
    }
    patient_field <- if (max(score) >= 1) names(ids)[which.max(score)]
                     else "title (time-point words removed)"
    patient <- if (patient_field %in% names(ids)) ids[[patient_field]] else rep("", nrow(ann))
  } else {
    patient <- ann[[patient_field]]
  }
  if (!is.null(patient_regex)) {
    m <- regmatches(patient, regexec(patient_regex, patient, perl = TRUE))
    patient <- vapply(m, function(z) if (length(z) >= 2L) z[2L] else "", character(1))
    patient_field <- paste0(patient_field, " ~ ", patient_regex)
    if (!any(nzchar(patient))) {
      .ptpm_geo_fail(ann, "'patient_regex' did not match any value of the patient field.")
    }
  }

  samples <- data.frame(sample = colnames(ex), patient = patient, condition = condition,
                        arm = ifelse(is_control, "control", "treated"), ann,
                        check.names = FALSE, stringsAsFactors = FALSE)
  rownames(samples) <- NULL

  usable <- !is.na(condition) & nzchar(patient)
  n_treated_types <- length(unique(ann[[time_field]][condition %in% "treated"]))
  dup <- duplicated(paste(patient, condition)[usable]) |
    duplicated(paste(patient, condition)[usable], fromLast = TRUE)
  if (any(dup)) {
    .ptpm_geo_fail(samples[usable, ][dup, c("sample", "patient", "condition", ann_cols)],
                   "Some patients have more than one sample per time point with the detected ",
                   "fields (patient = '", patient_field, "', time = '", time_field, "'). ",
                   "Set 'patient_field' (or 'patient_regex') and 'time_field'; if the series ",
                   "has several treated time points, narrow 'treated_pattern' to one of them.",
                   if (n_treated_types > 1L) paste0(" Treated values found: ",
                     paste(utils::head(unique(ann[[time_field]][condition %in% "treated"]), 6),
                           collapse = " | ")) else "")
  }
  mk_design <- function(sel) {
    s <- samples[sel & usable, , drop = FALSE]
    complete <- intersect(s$patient[s$condition == "baseline"], s$patient[s$condition == "treated"])
    s <- s[s$patient %in% complete, , drop = FALSE]
    if (!nrow(s)) return(NULL)
    data.frame(sample = s$sample, pair = s$patient, condition = s$condition,
               patient = s$patient, stringsAsFactors = FALSE)
  }
  design <- mk_design(!is_control)
  design_control <- mk_design(is_control)
  n_pairs <- if (is.null(design)) 0L else length(unique(design$pair))
  n_ctrl <- if (is.null(design_control)) 0L else length(unique(design_control$pair))
  message("PTPM: fields used - time: '", time_field, "', patient: '", patient_field,
          "', arm: '", if (is.na(arm_field)) "none" else arm_field, "'.")
  message("PTPM: ", n_pairs, " complete baseline/treated pairs in the treated arm",
          if (n_ctrl) paste0(", ", n_ctrl, " in the control arm") else "", "; ",
          sum(!usable), " sample(s) not classified.")
  if (n_pairs < k) {
    .ptpm_geo_fail(samples[, c("sample", "patient", "condition", "arm", ann_cols)],
                   "Only ", n_pairs, " complete pairs were found (k = ", k, " needed). ",
                   "Check the classification above and set the fields explicitly.")
  }

  expression <- data.frame(gene_id = rownames(ex), ex, check.names = FALSE,
                           stringsAsFactors = FALSE)
  rownames(expression) <- NULL
  list(expression = expression, design = design, design_control = design_control,
       samples = samples,
       fields = list(time = time_field, patient = patient_field,
                     arm = if (is.na(arm_field)) NULL else arm_field,
                     symbol = symbol_field, log2_transformed = do_log))
}

.ptpm_geo_fail <- function(tab, ...) {
  op <- options(width = 200)
  on.exit(options(op))
  message("Sample annotation (first rows):")
  print(utils::head(tab, 12), row.names = FALSE)
  stop(..., call. = FALSE)
}

# Patient identifiers commonly encoded in GEO sample titles. Each candidate is
# scored like an annotation field (complete pairs, at most one sample per
# patient and time point), so a wrong guess is simply not selected.
#   "10A;pretreatment;..."  -> first token "10A", without trailing letter "10"
#   "BC_2" / "OC_2"         -> trailing number "2"
#   "P07 baseline"          -> first token "P07"
.ptpm_title_id_candidates <- function(title) {
  if (is.null(title)) return(list())
  first <- sub("^\\s*([^;,|:[:space:]]+).*$", "\\1", title)
  num_of <- function(x, re) ifelse(grepl(re, x, perl = TRUE), sub(re, "\\1", x, perl = TRUE), "")
  list(
    `title: first token` = first,
    `title: first token without trailing letter` = num_of(first, "^([A-Za-z]*[0-9]+)[A-Za-z]$"),
    `title: leading number` = num_of(title, "^\\D*?(\\d+).*$"),
    `title: trailing number` = num_of(title, "^.*?(\\d+)\\D*$")
  )
}
