# ============================================================
# Script: run_liana.R
#
# Purpose:
# This script runs LIANA using the scording methods sca, cellphonedb, 
# and cellchat
#
# Reference:
# https://saezlab.github.io/liana/articles/liana_tutorial.html
#
# Output:
# LIANA objects and CSV files are saved to:
# the folder set in config.R/liana
# ============================================================

# install if needed
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("SummarizedExperiment")
# remotes::install_github("saezlab/liana")
# BiocManager::install("sparseMatrixStats")
# BiocManager::install("DelayedMatrixStats")

# load pipeline settings
source("config.R")
check_pipeline_config()

# load packages
library(Seurat)
library(liana)
library(SingleCellExperiment)
library(dplyr)

# set seed for reproducibility
set.seed(pipeline_seed)


# load the configured Seurat object and retain the requested groups
seurat_interest <- load_pipeline_seurat()

# output folder
liana_results_folder <- file.path(results_folder, "liana")

# split into one or more conditions
condition_objects <- split_pipeline_conditions(seurat_interest)
condition_info <- get_pipeline_conditions(seurat_interest)

# function to run LIANA on one Seurat object
run_liana_one_object <- function(seurat_object, output_name, output_folder) {  
  cat("\nRunning LIANA for:", output_name, "\n")

  DefaultAssay(seurat_object) <- "RNA"

  seurat_object <- NormalizeData(seurat_object, assay = "RNA",>

  # use analysis groups instead of numeric Seurat clusters
  Idents(seurat_object) <- celltype_column
  
  # check cell identities
  table(Idents(seurat_object))
  
  # get cell type labels
  clusters <- as.character(Idents(seurat_object))
  
  # remove groups with fewer than 5 cells
  tab <- table(clusters)
  keep_clusters <- names(tab[tab >= 5])
  seurat_liana <- subset(seurat_object, idents = keep_clusters)
  
  # check filtered cell identities
  table(Idents(seurat_liana))
  
  # convert Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_liana, assay = "RNA")
  
  # add cell type labels
  sce$celltype <- as.character(Idents(seurat_liana))
  
  # run LIANA 
  liana_results <- liana_wrap(sce, method = c("cellphonedb", "sca"), resource = "Consensus", idents_col = "celltype", expr_prop = 0.10)
  
  # save the raw per-method results, in case individual methods need inspecting
  saveRDS(liana_results, file.path(output_folder, paste0(output_name, "_liana_raw_by_method.rds")))
  
  # aggregate across methods into one table
  liana_aggregated <- liana_aggregate(liana_results)
  
  # save all aggregated LIANA results
  write.csv(liana_aggregated, file.path(output_folder, paste0(output_name, "_liana_results.csv")), row.names = FALSE)
  saveRDS(liana_aggregated, file.path(output_folder, paste0(output_name, "_liana_results.rds")))
  
  # filter to top ranked interactions by using aggregate_rank
  liana_results_sig <- liana_aggregated[liana_aggregated$aggregate_rank < 0.05, ]
  
  # order by aggregate_rank, strongest first
  liana_results_sig <- liana_results_sig[order(liana_results_sig$aggregate_rank), ]
  
  # save significant LIANA results
  write.csv(liana_results_sig, file.path(output_folder, paste0(output_name, "_liana_significant_results.csv")), row.names = FALSE)
  
  # return significant and full results
  return(list(all = liana_aggregated, significant = liana_results_sig))
}

# run LIANA for every condition
liana_results <- list()
for (i in seq_len(nrow(condition_info))) {
  condition_name <- condition_info$condition[i]
  condition_folder <- condition_info$folder[i]
  output_folder <- file.path(liana_results_folder, condition_folder)
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  liana_results[[condition_name]] <- run_liana_one_object(seurat_object = condition_objects[[condition_name]], output_name = condition_folder, output_folder = output_folder)
}
