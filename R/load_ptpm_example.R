#' Load the PTPM example dataset
#'
#' Loads the example data used to demonstrate [run_ptpm()]:
#' the GSE204805 variance-stabilised expression matrix
#' (`GSE204805_vsd.tsv.gz`) and the matching sample metadata
#' (`PTPM_metadata_with_title.csv`, tab-separated).
#'
#' The expression matrix (about 206 MB) is too large to ship with the
#' package, so it is downloaded from the GEO supplementary files of series
#' GSE204805 the first time the function is called and stored in a local
#' cache directory. Later calls read the cached copy without downloading it
#' again. The metadata file is small and is shipped in `inst/extdata`.
#'
#' @param cache_dir Directory where the downloaded expression file is kept.
#'   Defaults to the user cache directory of the package as returned by
#'   [tools::R_user_dir()].
#' @param force_download Logical. If `TRUE`, the expression file is
#'   downloaded again even if a cached copy exists.
#' @param timeout Download timeout in seconds. The file is large, so the
#'   default R timeout of 60 seconds is usually not enough.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`expression`}{Data frame whose first column contains gene
#'       identifiers and whose remaining columns are samples.}
#'     \item{`metadata`}{Data frame of sample metadata including the `title`
#'       column.}
#'   }
#' @examples
#' \dontrun{
#' # The first call downloads about 206 MB from GEO
#' example_data <- load_ptpm_example()
#'
#' result <- run_ptpm(
#'   expression = example_data$expression,
#'   metadata = example_data$metadata
#' )
#' }
#' @export
load_ptpm_example <- function(cache_dir = tools::R_user_dir("PTPM", which = "cache"),
                              force_download = FALSE,
                              timeout = 3600) {
  meta_path <- system.file("extdata", "PTPM_metadata_with_title.csv", package = "PTPM")
  if (!nzchar(meta_path)) {
    stop("Example file 'PTPM_metadata_with_title.csv' not found in the installed PTPM package.")
  }

  # Expression matrix: downloaded from GEO and cached locally.
  expr_path <- .ptpm_get_vsd_file(cache_dir = cache_dir,
                                  force_download = force_download,
                                  timeout = timeout)

  # Metadata: tab-separated as in the original analysis; fall back to commas.
  metadata <- utils::read.table(meta_path, header = TRUE, sep = "\t",
                                quote = "\"", comment.char = "",
                                stringsAsFactors = FALSE, check.names = FALSE)
  if (!"title" %in% names(metadata)) {
    metadata <- utils::read.csv(meta_path, stringsAsFactors = FALSE,
                                check.names = FALSE)
  }
  if (!"title" %in% names(metadata)) {
    stop("The example metadata file does not contain a 'title' column.")
  }

  # Expression: keep the gene identifier as the first column, even when the
  # header has one field fewer than the data rows (as written by write.table).
  con <- gzfile(expr_path, open = "rt")
  first_lines <- readLines(con, n = 2L)
  close(con)
  header <- strsplit(first_lines[1L], "\t", fixed = TRUE)[[1L]]
  n_fields <- length(strsplit(first_lines[2L], "\t", fixed = TRUE)[[1L]])
  header <- gsub('^"|"$', "", header)
  if (length(header) == n_fields - 1L) {
    header <- c("gene_id", header)
  }
  if (length(header) != n_fields) {
    stop("Unexpected layout in 'GSE204805_vsd.tsv.gz': header and data rows ",
         "have different numbers of fields.")
  }

  expression <- utils::read.table(gzfile(expr_path), header = FALSE, skip = 1L,
                                  sep = "\t", quote = "\"", comment.char = "",
                                  col.names = header, check.names = FALSE,
                                  stringsAsFactors = FALSE)

  list(expression = expression, metadata = metadata)
}


# Internal helper: return the path to a local copy of GSE204805_vsd.tsv.gz,
# downloading it from the GEO supplementary files when needed.
.ptpm_get_vsd_file <- function(cache_dir, force_download = FALSE, timeout = 3600) {
  file_name <- "GSE204805_vsd.tsv.gz"
  dest <- file.path(cache_dir, file_name)

  if (file.exists(dest) && !isTRUE(force_download)) {
    return(dest)
  }

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/",
                "?acc=GSE204805&format=file&file=GSE204805%5Fvsd%2Etsv%2Egz")

  # The file is about 206 MB: raise the timeout for this call only.
  old_opts <- options(timeout = max(timeout, getOption("timeout")))
  on.exit(options(old_opts), add = TRUE)

  message("Downloading '", file_name, "' (about 206 MB) from GEO series GSE204805.\n",
          "This is done only once; the file will be cached in:\n  ", cache_dir)

  # Download to a temporary name so an interrupted download never leaves a
  # truncated file that would be mistaken for a valid cached copy.
  tmp <- paste0(dest, ".part")
  status <- tryCatch(
    utils::download.file(url, destfile = tmp, mode = "wb", quiet = FALSE),
    error = function(e) {
      unlink(tmp)
      stop("Could not download '", file_name, "' from GEO: ",
           conditionMessage(e), "\nCheck your internet connection, or download it ",
           "manually from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204805 ",
           "and place it in: ", cache_dir, call. = FALSE)
    }
  )
  if (!identical(status, 0L) || !file.exists(tmp) || file.size(tmp) == 0) {
    unlink(tmp)
    stop("Download of '", file_name, "' from GEO failed.", call. = FALSE)
  }

  if (!file.rename(tmp, dest)) {
    unlink(tmp)
    stop("Could not move the downloaded file to: ", dest, call. = FALSE)
  }

  dest
}
