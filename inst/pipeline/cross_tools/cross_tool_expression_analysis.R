# ============================================================
# Script: cross_tool_expression_analysis.R
#
# Purpose:
# This script identifies ligand-receptor pairs that are detected
# consistently across multiple tools.
#
# Output:
# Consensus pair table and expression dot plot saved to:
# the folder set in config.R/cross_tool_analysis
# ============================================================

# load pipeline settings
source("config.R")
check_pipeline_config()
condition_info <- get_pipeline_conditions()
conditions <- condition_info$condition

library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)

cross_tool_folder <- file.path(results_folder, "cross_tool_analysis")
dir.create(cross_tool_folder, recursive = TRUE, showWarnings = FALSE)

# load each tool's ranked interaction results
read_ranked <- function(file, tool) read.csv(file, stringsAsFactors = FALSE, check.names = FALSE) %>% select(condition, source, target, ligand, receptor, percentile) %>% mutate(tool = tool)
cellchat_pairs <- read_ranked(file.path(results_folder, "comparison/cellchat/cellchat_ranked_interactions.csv"), "CellChat")
cpdb_pairs <- read_ranked(file.path(results_folder, "comparison/cellphonedb/cellphonedb_ranked_interactions.csv"), "CellPhoneDB")
liana_pairs <- read_ranked(file.path(results_folder, "comparison/liana/liana_ranked_interactions.csv"), "LIANA")
scsr_pairs <- read_ranked(file.path(results_folder, "comparison/singlecellsignalr/singlecellsignalr_ranked_interactions.csv"), "SingleCellSignalR")

# combine all interactions and count how many tools detected each exact interaction
all_pairs <- bind_rows(cellchat_pairs, cpdb_pairs, liana_pairs, scsr_pairs) %>% mutate(ligand = trimws(toupper(ligand)), receptor = trimws(toupper(receptor)))
interaction_tool_counts <- all_pairs %>% distinct(condition, source, target, ligand, receptor, tool, .keep_all = TRUE) %>% group_by(condition, source, target, ligand, receptor) %>% summarise(n_tools = n_distinct(tool), tools = paste(sort(unique(tool)), collapse = ", "), mean_percentile = mean(percentile, na.rm = TRUE), maximum_percentile = max(percentile, na.rm = TRUE), .groups = "drop") %>% arrange(condition, desc(n_tools), desc(mean_percentile))
write.csv(interaction_tool_counts, file.path(cross_tool_folder, "exact_interaction_tool_counts.csv"), row.names = FALSE)

# also count ligand receptor pairs without requiring the same sender and receiver
pair_tool_counts <- all_pairs %>% distinct(condition, ligand, receptor, tool) %>% group_by(condition, ligand, receptor) %>% summarise(n_tools = n_distinct(tool), tools = paste(sort(unique(tool)), collapse = ", "), .groups = "drop") %>% arrange(condition, desc(n_tools))
write.csv(pair_tool_counts, file.path(cross_tool_folder, "ligand_receptor_pair_tool_counts.csv"), row.names = FALSE)

consensus_pairs <- pair_tool_counts %>% filter(n_tools >= 2)
top_consensus_pairs <- consensus_pairs %>% group_by(condition) %>% slice_head(n = 20) %>% ungroup()
write.csv(top_consensus_pairs, file.path(cross_tool_folder, "top_consensus_ligand_receptor_pairs.csv"), row.names = FALSE)

if (nrow(top_consensus_pairs) == 0) {
  cat("No ligand receptor pairs were detected by at least two tools.\n")
} else {
  seurat_interest <- load_pipeline_seurat()
  seurat_interest <- JoinLayers(seurat_interest, assay = "RNA")
  split_complex_genes <- function(gene_string) unlist(strsplit(gene_string, "_"))
  ligand_genes <- unique(unlist(sapply(top_consensus_pairs$ligand, split_complex_genes)))
  receptor_genes <- unique(unlist(sapply(top_consensus_pairs$receptor, split_complex_genes)))
  seurat_gene_names <- rownames(seurat_interest)
  ligands_to_plot <- seurat_gene_names[toupper(seurat_gene_names) %in% toupper(ligand_genes)]
  receptors_to_plot <- seurat_gene_names[toupper(seurat_gene_names) %in% toupper(receptor_genes)]
  Idents(seurat_interest) <- celltype_column
  dotplot_args <- list(object = seurat_interest, group.by = celltype_column, cols = "RdYlBu")
  if (!is.null(condition_column)) dotplot_args$split.by <- condition_column
  if (length(ligands_to_plot) > 0) {
    dotplot_args$features <- ligands_to_plot
    ligand_dotplot <- do.call(DotPlot, dotplot_args) + theme_classic() + ggtitle("Expression of consensus ligand genes") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9), axis.text.y = element_text(size = 9), plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + labs(x = "Ligand gene", y = "Cell group")
    ggsave(file.path(cross_tool_folder, "consensus_ligand_genes_expression_dotplot.png"), plot = ligand_dotplot, width = 14, height = 10, dpi = 300)
  }
  if (length(receptors_to_plot) > 0) {
    dotplot_args$features <- receptors_to_plot
    receptor_dotplot <- do.call(DotPlot, dotplot_args) + theme_classic() + ggtitle("Expression of consensus receptor genes") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9), axis.text.y = element_text(size = 9), plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + labs(x = "Receptor gene", y = "Cell group")
    ggsave(file.path(cross_tool_folder, "consensus_receptor_genes_expression_dotplot.png"), plot = receptor_dotplot, width = 14, height = 10, dpi = 300)
  }
}
