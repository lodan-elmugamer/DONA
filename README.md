# DONA

**DONA** is a ligand-receptor interaction pipeline for annotated Seurat objects. It brings together CellChat, CellPhoneDB, LIANA and SingleCellSignalR, then produces tool-specific rankings, plots, expression checks and cross-tool consensus results.

This submission uses the anonymous author identifier **B288659**. The placeholder email in `DESCRIPTION` is deliberately non-functional and should be replaced before public release.

## What it does

- Accepts one or more biological conditions.
- Runs each supported interaction method separately.
- Ranks and plots interactions within each condition.
- Finds interactions supported by multiple tools.
- Performs condition-change analysis when exactly two conditions are present.

## Installation from a local folder

```r
install.packages("devtools")
devtools::install("path/to/DONA_R_package")
```

## Start a new analysis

```r
library(DONA)
dona_copy_pipeline("my_DONA_analysis")
```

Open `my_DONA_analysis/config.R` and change the file paths and metadata column names.

Then run:

```r
dona_run("my_DONA_analysis")
```

Run CellPhoneDB separately after configuring its Python environment:

```r
dona_run_cellphonedb("my_DONA_analysis", env = c(
  VENV_ACTIVATE = "/path/to/venv/bin/activate",
  CPDB_ZIP = "/path/to/cellphonedb.zip"
))
```

Finally run ranking, plotting and consensus analysis:

```r
dona_run_downstream("my_DONA_analysis")
```

## Important dependency note

The package wrapper itself is lightweight, but the analysis scripts require the relevant scientific packages and CellPhoneDB Python environment. These must be installed separately because some are distributed through Bioconductor, GitHub or Python rather than CRAN.

## Anonymous grading

Before public release, replace `B288659` and `redacted@example.invalid` in `DESCRIPTION`, `LICENSE` and `CITATION.cff` with the official author details.
