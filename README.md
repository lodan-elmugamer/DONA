# DONA - Decoding Overlapping Networks for Ligand–Receptor Analysis

**DONA** is a ligand-receptor interaction pipeline for annotated single-cell RNA-seq Seurat objects. It uses CellChat, CellPhoneDB, LIANA and SingleCellSignalR, then produces tool-specific rankings, plots, expression checks and cross-tool consensus results.

The B number **B288659** will be used as the author identifier for the time being.

## What it does

- Accepts one or more biological conditions.
- Runs each interaction method separately.
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

