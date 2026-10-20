% Generated for PTPM 0.1.1
\name{load_ptpm_example}
\alias{load_ptpm_example}
\title{Load the PTPM example dataset}
\usage{
load_ptpm_example()
}
\value{
A list with two elements:
\describe{
\item{\code{expression}}{Data frame whose first column contains gene
identifiers and whose remaining columns are samples.}
\item{\code{metadata}}{Data frame of sample metadata including the \code{title}
column.}
}
}
\description{
Reads the example data shipped with the package in \code{inst/extdata}:
the GSE204805 variance-stabilised expression matrix
(\code{GSE204805_vsd.tsv.gz}) and the matching sample metadata
(\code{PTPM_metadata_with_title.csv}, tab-separated). The returned objects are
already in the format expected by \code{\link[=run_ptpm]{run_ptpm()}}.
}
\examples{
\dontrun{
example_data <- load_ptpm_example()

result <- run_ptpm(
  expression = example_data$expression,
  metadata = example_data$metadata
)
}
}
