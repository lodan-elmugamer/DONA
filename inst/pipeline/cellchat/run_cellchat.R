# ============================================================
# Script: run_cellchat.R
#
# Purpose:
# This script runs CellChat separately for each configured
# condition using the processed Seurat object with analysis groups.
#
# Reference:
# https://rpubs.com/HHJ/921311
# https://github.com/jinworks/CellChat/blob/main/R/CellChat_class.R
# https://github.com/jinworks/CellChat/blob/main/tutorial/CellChat_analysis_of_spatial_transcriptomics_data.Rmd
#
# Output:
# CellChat objects, CSV files and PNGs saved to:
# the folder set in config.R/cellchat
# ============================================================

# install the packages as needed
# install.packages("devtools")
# devtools::install_github("jinworks/CellChat")

# load pipeline settings
source("config.R")
check_pipeline_config()

# load the libraries
library(Seurat)
library(CellChat)
library(patchwork)
library(viridis)

options(stringsAsFactors = FALSE)

# set seed for reproducibility
set.seed(pipeline_seed)


# load the configured Seurat object and retain the requested groups
seurat_interest <- load_pipeline_seurat()

# output folder
cellchat_results_folder <- file.path(results_folder, "cellchat")

# split into one or more conditions
condition_objects <- split_pipeline_conditions(seurat_interest)
condition_info <- get_pipeline_conditions(seurat_interest)

# function to run CellChat on one Seurat object
run_cellchat_one_object <- function(state, seurat_obj, output_name, output_folder) {
  
  # print progress
  cat("\nRunning CellChat for:", output_name, "\n")
  
  # set analysis group as identity
  Idents(seurat_obj) <- celltype_column
  
  # get normalized expression matrix
  data.input <- LayerData(seurat_obj, assay = "RNA", layer = "data")

  # get metadata in the same order as the expression matrix
  meta <- seurat_obj@meta.data
  meta <- meta[colnames(data.input), , drop = FALSE]

  labels <- as.character(meta[[celltype_column]])

  keep <- !is.na(labels) & labels != ""

  data.input <- data.input[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  meta[[celltype_column]] <- droplevels(factor(meta[[celltype_column]]))

  rownames(meta) <- colnames(data.input)
  
  # create CellChat object from Seurat object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = celltype_column)
  cellchat@idents <- droplevels(cellchat@idents)  

  # use human CellChat database
  CellChatDB <- CellChatDB.human
  
  # use secreted signaling, ecm receptor and cell-cell contact
  CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
  
  # add database to CellChat object
  cellchat@DB <- CellChatDB.use
  
  # subset expression data to signaling genes
  cellchat <- subsetData(cellchat)
  
  # identify overexpressed genes
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  
  # identify overexpressed ligand receptor interactions
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # calculate communication probability
  cellchat <- computeCommunProb(cellchat)
  
  # filter weak interactions
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # calculate communication at pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  # aggregate communication network
  cellchat <- aggregateNet(cellchat)
  
  # save CellChat object
  saveRDS(cellchat, file.path(output_folder, paste0(output_name, "_cellchat.rds")))
  
  # extract interaction table
  communication_table <- subsetCommunication(cellchat)
  
  # save interaction table
  write.csv(communication_table, file.path(output_folder, paste0(output_name, "_cellchat_interactions.csv")), row.names = FALSE)
  
  # create pathway summary table
  pathway_summary <- data.frame(pathway_name = names(tapply(communication_table$prob, communication_table$pathway_name, sum)), total_prob = as.numeric(tapply(communication_table$prob, communication_table$pathway_name, sum)), mean_prob = as.numeric(tapply(communication_table$prob, communication_table$pathway_name, mean)), n = as.numeric(table(communication_table$pathway_name)))
  
  # order pathway summary by total probability
  pathway_summary <- pathway_summary[order(pathway_summary$total_prob, decreasing = TRUE), ]
  
  # save pathway summary table
  write.csv(pathway_summary, file.path(output_folder, paste0(output_name, "_cellchat_pathway_summary.csv")), row.names = FALSE)
  
  # return CellChat object
  return(cellchat)
}

# run CellChat for every condition
cellchat_results <- list()
for (i in seq_len(nrow(condition_info))) {
  condition_name <- condition_info$condition[i]
  condition_folder <- condition_info$folder[i]
  condition_object <- condition_objects[[condition_name]]
  condition_object <- JoinLayers(condition_object, assay = "RNA")
  table(condition_object@meta.data[[celltype_column]])
  output_folder <- file.path(cellchat_results_folder, condition_folder)
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  cellchat_results[[condition_name]] <- run_cellchat_one_object(condition_name, seurat_obj = condition_object, output_name = condition_folder, output_folder = output_folder)

  # set colors
  group_colours <- viridis(length(levels(cellchat_results[[condition_name]]@idents)))
  names(group_colours) <- levels(cellchat_results[[condition_name]]@idents)

  # plot number of interactions
  png(filename = file.path(output_folder, paste0(condition_folder, "_interactions.png")), width = 2400, height = 1800, res = 300)

  # title would not get saved with the diagram so a panel is needed
  layout(matrix(c(1, 2), nrow = 2), heights = c(0.15, 0.85))

  # title panel
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(x = 0.5, y = 0.5, labels = paste(condition_name, "Communication Network"), cex = 1.6, font = 2)

  # circle plot panel
  netVisual_circle(cellchat_results[[condition_name]]@net$count, vertex.weight = as.numeric(table(cellchat_results[[condition_name]]@idents)), color.use = group_colours, weight.scale = TRUE, label.edge = FALSE, title.name = paste(condition_name, "Communication Network"))
  dev.off()
}
