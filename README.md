# PTPM: Paired Transcriptomic Perturbation Modeling

[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Version](https://img.shields.io/badge/version-0.1.5-informational)

**PTPM** is an R package implementing a protocol for estimating the
predictability of a paired transcriptomic perturbation phenotype from
baseline profiles.

## Features

- Patient-grouped nested cross-validation
- Fold-internal feature selection
- Fold-internal adjustment for mathematical coupling between baseline
  predictors and the paired outcome
- 2 x 2 decomposition of the performance gap
- Patient-level bootstrap confidence intervals
- Sensitivity analyses on pathway scoring (training-reference GSVA, ssGSEA),
  coupling and pathway circularity

## Installation

PTPM requires R (>= 4.1). Some dependencies (e.g. `GSVA`) are distributed
through Bioconductor, so install them first:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("GSVA")

remotes::install_github("franruss/PTPM")
```

## Dependencies

`caret`, `doParallel`, `foreach`, `glmnet`, `GSVA`, `msigdbr`, `parallel`,
`stats`, `utils`. Optional: `Rcpp`.

## Usage

```r
library(PTPM)

# List the documented functions of the package
help(package = "PTPM")
```

## Citation

If you use PTPM in your work, please cite:

> Russo, F. (2026). *PTPM: Paired Transcriptomic Perturbation Modeling*
> (Version 0.1.5) [R package]. Zenodo.

```bibtex
@software{russo_ptpm_2026,
  author    = {Russo, Francesco},
  title     = {PTPM: Paired Transcriptomic Perturbation Modeling},
  year      = {2026},
  version   = {0.1.5},
  publisher = {Zenodo}
}
```

Citation metadata is also available in [`CITATION.cff`](CITATION.cff).

## Author

**Francesco Russo**
ORCID: [0000-0002-7017-8344](https://orcid.org/0000-0002-7017-8344)
Contact: francesco.russo@cnr.it

## License

This package is released under the GNU General Public License v3.0
(GPL-3). See [`LICENSE`](LICENSE) for the full text.
