# ============================================================
# Script: export_cellphonedb_inputs.R
#
# Purpose:
# This script exports CellPhoneDB input files separately for
# diagnosis and relapse samples using the processed sample based
# Seurat object with analysis groups.
#
# Output:
# metadata.txt and sparse matrix files saved to:
# the folder set in config.R/cellphonedb/diagnosis/input
# the folder set in config.R/cellphonedb/relapse/input
# ============================================================

# load pipeline settings
source("config.R")
check_pipeline_config()

# load packages
library(Seurat)
library(Matrix)


# load the configured Seurat object and retain the requested groups
seurat_interest <- load_pipeline_seurat()

# split into one or more conditions
condition_objects <- split_pipeline_conditions(seurat_interest)
condition_info <- get_pipeline_conditions(seurat_interest)

# export each condition
for (i in seq_len(nrow(condition_info))) {
  condition_name <- condition_info$condition[i]
  condition_folder <- condition_info$folder[i]
  condition_object <- condition_objects[[condition_name]]
  input_folder <- file.path(results_folder, "cellphonedb", condition_folder, "input")
  dir.create(input_folder, recursive = TRUE, showWarnings = FALSE)

  # check analysis groups
  table(condition_object@meta.data[[celltype_column]])

  # export metadata
  metadata <- data.frame(Cell = colnames(condition_object), cell_type = condition_object@meta.data[[celltype_column]])
  write.table(metadata, file.path(input_folder, "metadata.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  # get normalized counts
  counts <- GetAssayData(condition_object, assay = "RNA", layer = "data")

  # export counts as sparse matrix
  writeMM(counts, file.path(input_folder, "matrix.mtx"))
  write.table(rownames(counts), file.path(input_folder, "features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(colnames(counts), file.path(input_folder, "barcodes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
