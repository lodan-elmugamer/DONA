# Simple package submission guide

## What has already been created

- Standard R package folders.
- `DESCRIPTION`, `NAMESPACE` and MIT licence files.
- User functions for copying and running the pipeline.
- Help comments, tests, a vignette, README and NEWS.
- A pkgdown website configuration.
- Anonymous author details using B288659.

## Check it on the university computer

Install the helper packages:

```r
install.packages(c("devtools", "roxygen2", "testthat", "pkgdown", "rmarkdown", "knitr"))
```

From inside the package folder run:

```r
devtools::document()
devtools::test()
devtools::check()
```

`devtools::document()` creates the final help files in `man/` and refreshes `NAMESPACE`.

## Install it locally

```r
devtools::install()
library(LRIP)
```

## Build the documentation website

```r
pkgdown::build_site()
```

## Create a package file for submission

```r
devtools::build()
```

This creates a file similar to `LRIP_0.1.0.tar.gz`.

## Before a public GitHub release

Replace the anonymous details in:

- `DESCRIPTION`
- `LICENSE`
- `CITATION.cff`

Then create a GitHub repository and push the package. Users could install it with:

```r
remotes::install_github("USERNAME/LRIP")
```

## About CRAN or Bioconductor

The current version is suitable as an installable university software submission and GitHub package. A CRAN or Bioconductor submission would need additional dependency handling, automated tests using small example data and checks on several operating systems.
