# ============================================================
# Script: run_singlecellsignalr.R
#
# Purpose:
# This script runs SingleCellSignalR separately for diagnosis
# and relapse samples using the processed sample based Seurat
# object with analysis groups.
#
# Reference:
# https://bioconductor.org/packages//release/bioc/vignettes/SingleCellSignalR/inst/doc/SingleCellSignalR-Main.html
#
# Output:
# SingleCellSignalR objects and CSV files saved to:
# the folder set in config.R/singlecellsignalr
# ============================================================

# install packages as needed
# BiocManager::install("SingleCellSignalR")

# load pipeline settings
source("config.R")
check_pipeline_config()

# load the libraries
library(BulkSignalR)
library(SingleCellSignalR)
library(Seurat)


# load the configured Seurat object and retain the requested groups
seurat_interest <- load_pipeline_seurat()

# output folder
singlecellsignalr_results_folder <- file.path(results_folder, "singlecellsignalr")

# split into one or more conditions
condition_objects <- split_pipeline_conditions(seurat_interest)
condition_info <- get_pipeline_conditions(seurat_interest)

# function to combine SingleCellSignalR result lists into one data frame
combine_scsr_results <- function(result_list, id_column_name) {
  
  # return empty data frame if there are no results
  if (length(result_list) == 0) {
    return(data.frame())
  }
  
  # add the list name to each result table
  result_tables <- lapply(names(result_list), function(result_name) {
    result_table <- as.data.frame(result_list[[result_name]])
    result_table[[id_column_name]] <- result_name
    return(result_table)
  })
  
  # combine all result tables
  combined_results <- do.call(rbind, result_tables)
  
  # reset row numbers
  rownames(combined_results) <- NULL
  
  # return combined results
  return(combined_results)
}

# function to run SingleCellSignalR on one Seurat object
run_singlecellsignalr_one_object <- function(state, seurat_object, output_name, output_folder) {
  
  # print progress
  cat("\nRunning SingleCellSignalR for:", output_name, "\n")
  
  # create output folder
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  # use analysis groups instead of numeric Seurat clusters
  Idents(seurat_object) <- celltype_column
  
  # get cell type labels
  clusters <- as.character(Idents(seurat_object))
  
  # remove groups with fewer than 5 cells
  tab <- table(clusters)
  keep_clusters <- names(tab[tab >= 5])
  seurat_scsr <- subset(seurat_object, idents = keep_clusters)
  
  # join layers because the Seurat object has multiple RNA layers
  seurat_scsr <- JoinLayers(seurat_scsr, assay = "RNA")
  
  # get expression matrix
  mat <- GetAssayData(seurat_scsr, assay = "RNA", layer = "data")
  mat <- as.matrix(mat)
  
  # convert from natural log scale to log2 scale
  mat <- mat / log(2)
  
  # get cell type labels after filtering
  clusters <- as.character(Idents(seurat_scsr))
  
  # check matrix and labels match
  table(clusters)
  length(clusters)
  ncol(mat)
  stopifnot(length(clusters) == ncol(mat))
  
  # run SingleCellSignalR LRscore method
  scsr <- SCSRNoNet(mat, normalize = FALSE, method = "log-only", min.count = 1, prop = 0.001, log.transformed = TRUE, populations = clusters)
  
  # perform ligand receptor inference
  scsr <- performInferences(scsr, verbose = TRUE, min.logFC = 1e-10, max.pval = 1, min.LR.score = 0.5)
  
  # extract autocrine and paracrine interactions
  autocrines_results <- autocrines(scsr)
  paracrines_results <- paracrines(scsr)
  
  # combine autocrine and paracrine results into data frames
  autocrines_df <- combine_scsr_results(autocrines_results, "population")
  paracrines_df <- combine_scsr_results(paracrines_results, "population_pair")
  
  # all results, ordered by LRscore
  autocrines_by_score <- autocrines_df[order(-autocrines_df$LR.score),]
  paracrines_by_score <- paracrines_df[order(-paracrines_df$LR.score),]
  
  # significant results only, ordered by LRscore
  autocrines_significant <- autocrines_df[autocrines_df$pval < 0.05,]
  autocrines_significant <- autocrines_significant[order(-autocrines_significant$LR.score),]
  paracrines_significant <- paracrines_df[ paracrines_df$pval < 0.05,]
  paracrines_significant <- paracrines_significant[order(-paracrines_significant$LR.score),]
  
  write.csv(autocrines_df,  file.path(output_folder, paste0(output_name, "_singlecellsignalr_autocrines.csv")),  row.names = FALSE)
  write.csv(paracrines_df, file.path(output_folder, paste0(output_name, "_singlecellsignalr_paracrines.csv")), row.names = FALSE)
  
  
  write.csv(autocrines_by_score,  file.path(output_folder, paste0(output_name, "_singlecellsignalr_autocrines_by_score.csv")),  row.names = FALSE)
  write.csv(paracrines_by_score, file.path(output_folder, paste0(output_name, "_singlecellsignalr_paracrines_by_score.csv")), row.names = FALSE)
  write.csv(autocrines_significant,file.path(output_folder, paste0(output_name, "_singlecellsignalr_autocrines_significant.csv")),  row.names = FALSE)
  write.csv(paracrines_significant,  file.path(output_folder, paste0(output_name, "_singlecellsignalr_paracrines_significant.csv") ), row.names = FALSE)

  saveRDS(scsr, file.path(output_folder, paste0(output_name, "_singlecellsignalr_object.rds")))
  
  # return SingleCellSignalR object
  return(scsr)
}

# run SingleCellSignalR for every condition
singlecellsignalr_results <- list()
for (i in seq_len(nrow(condition_info))) {
  condition_name <- condition_info$condition[i]
  condition_folder <- condition_info$folder[i]
  condition_object <- condition_objects[[condition_name]]
  table(condition_object@meta.data[[celltype_column]])
  output_folder <- file.path(singlecellsignalr_results_folder, condition_folder)
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  singlecellsignalr_results[[condition_name]] <- run_singlecellsignalr_one_object(condition_name, seurat_object = condition_object, output_name = condition_folder, output_folder = output_folder)
}
