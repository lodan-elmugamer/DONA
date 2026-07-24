# ============================================================
# Script: plot_cellchat_ligand_receptors.R
#
# Purpose:
# This script plots CellChat ligand receptor interactions.
#
# Reference:
# https://r-graph-gallery.com/320-the-basis-of-bubble-plot.html
#
# Output:
# PNGs and CSV tables saved to:
# the folder set in config.R/cellchat/ligand_receptor_plots
# ============================================================

# load pipeline settings
source("config.R")
check_pipeline_config()
condition_info <- get_pipeline_conditions()

# load libraries
library(CellChat)
library(ggplot2)
library(viridis)

# set base folder
base_folder <- file.path(results_folder, "cellchat")
plot_folder <- file.path(base_folder, "ligand_receptor_plots")
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# load CellChat objects and pathway summaries
cellchat_objects <- setNames(lapply(seq_len(nrow(condition_info)), function(i) readRDS(file.path(base_folder, condition_info$folder[i], paste0(condition_info$folder[i], "_cellchat.rds")))), condition_info$condition)
pathway_summaries <- lapply(seq_len(nrow(condition_info)), function(i) read.csv(file.path(base_folder, condition_info$folder[i], paste0(condition_info$folder[i], "_cellchat_pathway_summary.csv"))))
pathways_to_plot <- unique(unlist(lapply(pathway_summaries, function(x) x$pathway_name)))

clean_label <- function(x) {
  x <- gsub("Normal_NK_fixed", "Normal_NK", x)
  x <- gsub("Normal_", "", x)
  x <- gsub("Malignant_", "", x)
  x <- gsub("_", " ", x)
  return(x)
}

# loop through selected pathways
for (pathway in pathways_to_plot) {
  pathway_tables <- lapply(names(cellchat_objects), function(condition_name) {
    cellchat_object <- cellchat_objects[[condition_name]]
    if (!pathway %in% cellchat_object@netP$pathways) return(data.frame())
    pathway_table <- subsetCommunication(cellchat_object, signaling = pathway)
    pathway_table$condition <- condition_name
    pathway_table
  })
  pathway_table <- do.call(rbind, pathway_tables)
  if (nrow(pathway_table) == 0) next
  pathway_folder <- file.path(plot_folder, pathway)
  dir.create(pathway_folder, recursive = TRUE, showWarnings = FALSE)
  write.csv(pathway_table, file.path(pathway_folder, paste0(pathway, "_combined_ligand_receptor_table.csv")), row.names = FALSE)
  pathway_table$source_label <- clean_label(pathway_table$source)
  pathway_table$target_label <- clean_label(pathway_table$target)
  pathway_table$cell_pair <- paste(pathway_table$source_label, pathway_table$target_label, sep = " -> ")
  pathway_table$ligand_receptor <- paste(pathway_table$ligand, pathway_table$receptor, sep = " - ")
  bubble_plot <- ggplot(pathway_table, aes(x = cell_pair, y = ligand_receptor, colour = prob)) + geom_point(alpha = 0.7, size = 4) + scale_colour_viridis_c(option = "viridis", name = "Interaction score") + facet_wrap(~condition, ncol = 1, scales = "free_x") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 9), plot.title = element_text(hjust = 0.5)) + labs(title = paste0(pathway, " ligand receptor interactions"), x = "Sender -> receiver", y = "Ligand receptor pair")
  ggsave(filename = file.path(pathway_folder, paste0(pathway, "_condition_bubble.png")), plot = bubble_plot, width = 14, height = max(6, 4 * nrow(condition_info)), dpi = 300)
}

# pathway interaction counts
pathway_counts <- do.call(rbind, lapply(seq_len(nrow(condition_info)), function(i) data.frame(pathway_name = pathway_summaries[[i]]$pathway_name, n = pathway_summaries[[i]]$n, condition = condition_info$condition[i])))
count_plot <- ggplot(pathway_counts, aes(x = reorder(pathway_name, n), y = n, fill = condition)) + geom_col(position = "dodge") + coord_flip() + theme_classic() + labs(title = "CellChat pathway interaction counts", x = "Pathway", y = "Number of interactions")
ggsave(filename = file.path(plot_folder, "cellchat_pathway_interaction_counts.png"), plot = count_plot, width = 10, height = 8, dpi = 300)
