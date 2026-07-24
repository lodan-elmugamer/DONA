# Ligand receptor pipeline refactor

This folder is a reusable copy of the original ligand receptor scripts. The scientific calculations, comments and general formatting were retained. Only dataset-specific paths, metadata names, condition values and repeated input preparation were changed where needed.

## Included tools

- CellChat
- CellPhoneDB
- LIANA
- SingleCellSignalR
- Cross-tool comparison

NicheNet and the ligand-receptor database comparison folder are not included.

## Condition support

The inference stage accepts one or more conditions.

- For one condition, set `condition_column <- NULL`. The full object is analysed once and saved under `single_condition`.
- For multiple conditions, set `condition_column` to the metadata column name and leave `condition_names <- NULL` to run all values.
- To run selected values only, set `condition_names <- c("Control", "Treatment")`.
- The comparative plotting, ranking and cross-tool scripts currently require exactly two conditions.

## Before running

Open `config.R` and update:

1. `seurat_file`
2. `output_folder`
3. `celltype_column`
4. `condition_column`
5. Optional `condition_names`
6. Optional cell groups and focused groups

Run every command from this main folder.

## Run the inference stage

```bash
Rscript run_pipeline.R
```

CellPhoneDB requires its Python environment and database. Update `CPDB_ZIP` and `VENV_ACTIVATE` at the top of `cellphonedb/run_cellphonedb.sh` if needed, then run:

```bash
bash cellphonedb/run_cellphonedb.sh
```

For exactly two conditions, run the comparative outputs with:

```bash
Rscript downstream_results.R
```

For one condition, do not run `downstream_results.R`. The individual tool results are already saved in the output folders.
