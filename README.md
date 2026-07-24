# DONA - Decoding Overlapping Networks for ligand–receptor Analysis

**DONA** is a ligand-receptor interaction pipeline for annotated single-cell RNA-seq Seurat objects. It brings together CellChat, CellPhoneDB, LIANA and SingleCellSignalR, then produces tool-specific rankings, plots, expression checks and cross-tool consensus results.

This submission uses the anonymous author identifier **B288659**. The placeholder email in `DESCRIPTION` will be replaced before public release.

## What it does

- Accepts one or more biological conditions.
- Runs each supported interaction method separately.
- Ranks and plots interactions within each condition.
- Finds interactions supported by multiple tools.
- Performs condition-change analysis when exactly two conditions are present.

## Installation from GitHub

```r
install.packages("remotes")
remotes::install_github("B288659-2025/DONA")
library(DONA)
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

The analysis scripts require the relevant packages and CellPhoneDB Python environment. These must be installed separately because some are distributed through Bioconductor, GitHub or Python rather than CRAN.

