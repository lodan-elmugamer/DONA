# ============================================================
# Script: find_multi_tool_overlaps_and_expression.R
#
# Purpose:
# This script finds exact ligand receptor interactions detected
# by multiple tools, checks their expression in the correct
# source and target cell groups, and ranks them using both
# cross tool support and gene expression.
#
# Output:
# Results are saved to:
# the folder set in config.R/comparison/overlaps
# ============================================================

# load pipeline settings
source("config.R")
check_pipeline_config()
condition_info <- get_pipeline_conditions()
conditions <- condition_info$condition

# load packages
library(Seurat)
library(dplyr)
library(ggplot2)

# set folders
comparison_folder <- file.path(results_folder, "comparison")
output_folder <- file.path(comparison_folder, "overlaps")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# load ranked interactions in long format
all_tools <- read.csv(file.path(comparison_folder, "all_tools_ranked_interactions_long.csv"), stringsAsFactors = FALSE, check.names = FALSE)

# find exact interactions detected by multiple tools within each condition
overlaps <- all_tools %>% group_by(condition, source, target, ligand, receptor) %>% summarise(method_count = n_distinct(method), standalone_count = n_distinct(method[method != "LIANA"]), methods = paste(sort(unique(method)), collapse = "; "), mean_percentile = mean(percentile, na.rm = TRUE), maximum_percentile = max(percentile, na.rm = TRUE), liana_supported = "LIANA" %in% method, .groups = "drop") %>% filter(method_count >= 2) %>% arrange(condition, desc(standalone_count), desc(method_count), desc(mean_percentile))
overlaps$interaction_id <- paste(overlaps$condition, overlaps$source, overlaps$target, overlaps$ligand, overlaps$receptor, sep = " | ")
write.csv(overlaps, file.path(output_folder, "all_multi_tool_exact_overlaps.csv"), row.names = FALSE)

overlap_candidates <- overlaps %>% filter(standalone_count >= 2 | (standalone_count >= 1 & liana_supported))
if (nrow(overlap_candidates) == 0) {
  cat("No eligible multi tool exact overlaps were found.\n")
} else {
  seurat_object <- load_pipeline_seurat()
  seurat_object <- JoinLayers(seurat_object, assay = "RNA")
  DefaultAssay(seurat_object) <- "RNA"
  overlap_candidates$source_group <- ifelse(overlap_candidates$source == "Normal_NK", "Normal_NK_fixed", overlap_candidates$source)
  overlap_candidates$target_group <- ifelse(overlap_candidates$target == "Normal_NK", "Normal_NK_fixed", overlap_candidates$target)

  split_genes <- function(x) {
    x <- x[!is.na(x) & x != ""]
    unique(trimws(unlist(strsplit(x, "_", fixed = TRUE))))
  }

  genes <- unique(c(split_genes(overlap_candidates$ligand), split_genes(overlap_candidates$receptor)))
  missing_genes <- setdiff(genes, rownames(seurat_object))
  if (length(missing_genes) > 0) cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
  genes <- intersect(genes, rownames(seurat_object))
  if (length(genes) == 0) stop("None of the ligand or receptor genes were found in the Seurat object.")

  vars_to_fetch <- c(genes, celltype_column)
  if (!is.null(condition_column)) vars_to_fetch <- c(vars_to_fetch, condition_column)
  expression_data <- FetchData(seurat_object, vars = vars_to_fetch, layer = "data")
  if (is.null(condition_column)) expression_data$.pipeline_condition <- conditions[1] else expression_data$.pipeline_condition <- as.character(expression_data[[condition_column]])

  get_expression <- function(gene, group, condition) {
    cells <- expression_data[expression_data[[celltype_column]] == group & expression_data$.pipeline_condition == condition, , drop = FALSE]
    if (nrow(cells) == 0) return(data.frame(n_cells = 0, mean_expression = NA_real_, percent_expressing = NA_real_))
    data.frame(n_cells = nrow(cells), mean_expression = mean(cells[[gene]], na.rm = TRUE), percent_expressing = mean(cells[[gene]] > 0, na.rm = TRUE) * 100)
  }

  get_interaction_expression <- function(candidate) {
    ligand_genes <- intersect(split_genes(candidate$ligand), genes)
    receptor_genes <- intersect(split_genes(candidate$receptor), genes)
    ligand_rows <- bind_rows(lapply(ligand_genes, function(gene) get_expression(gene, candidate$source_group, candidate$condition) %>% mutate(gene = gene, gene_role = "Ligand", cell_group = candidate$source)))
    receptor_rows <- bind_rows(lapply(receptor_genes, function(gene) get_expression(gene, candidate$target_group, candidate$condition) %>% mutate(gene = gene, gene_role = "Receptor", cell_group = candidate$target)))
    bind_rows(ligand_rows, receptor_rows) %>% mutate(condition = candidate$condition, source = candidate$source, target = candidate$target, ligand = candidate$ligand, receptor = candidate$receptor, interaction_id = candidate$interaction_id, method_count = candidate$method_count, standalone_count = candidate$standalone_count, methods = candidate$methods, mean_percentile = candidate$mean_percentile, maximum_percentile = candidate$maximum_percentile, liana_supported = candidate$liana_supported)
  }

  overlap_expression <- bind_rows(lapply(seq_len(nrow(overlap_candidates)), function(i) get_interaction_expression(overlap_candidates[i, ])))
  clean_group <- function(x) {
    x <- gsub("Normal_", "", x)
    x <- gsub("Malignant_", "", x)
    gsub("_", " ", x)
  }
  overlap_expression <- overlap_expression %>% mutate(mean_expression = round(mean_expression, 3), percent_expressing = round(percent_expressing, 1), interaction = paste0(ligand, " -> ", receptor, " | ", clean_group(source), " -> ", clean_group(target), " | ", standalone_count, " standalone"), gene_label = paste(gene_role, gene, sep = ": "), condition = factor(condition, levels = conditions))
  write.csv(overlap_expression, file.path(output_folder, "all_multi_tool_overlap_expression.csv"), row.names = FALSE)

  expression_summary <- overlap_expression %>% group_by(condition, source, target, ligand, receptor, interaction_id) %>% summarise(minimum_percent_expressing = ifelse(all(is.na(percent_expressing)), NA_real_, min(percent_expressing, na.rm = TRUE)), mean_percent_expressing = ifelse(all(is.na(percent_expressing)), NA_real_, mean(percent_expressing, na.rm = TRUE)), minimum_mean_expression = ifelse(all(is.na(mean_expression)), NA_real_, min(mean_expression, na.rm = TRUE)), mean_expression_support = ifelse(all(is.na(mean_expression)), NA_real_, mean(mean_expression, na.rm = TRUE)), .groups = "drop")
  ranked_overlaps <- overlap_candidates %>% left_join(expression_summary, by = c("condition", "source", "target", "ligand", "receptor", "interaction_id")) %>% arrange(condition, desc(standalone_count), desc(minimum_percent_expressing), desc(mean_percent_expressing), desc(mean_percentile))
  write.csv(ranked_overlaps, file.path(output_folder, "all_multi_tool_overlaps_ranked_with_expression.csv"), row.names = FALSE)
  top_overlaps <- ranked_overlaps %>% group_by(condition) %>% slice_head(n = 15) %>% ungroup()
  write.csv(top_overlaps, file.path(output_folder, "top_multi_tool_exact_overlaps.csv"), row.names = FALSE)
  top_overlap_expression <- overlap_expression %>% filter(interaction_id %in% top_overlaps$interaction_id)
  write.csv(top_overlap_expression, file.path(output_folder, "top_multi_tool_overlap_expression.csv"), row.names = FALSE)
  overlap_plot <- ggplot(top_overlap_expression, aes(x = gene_label, y = interaction, size = percent_expressing, colour = mean_expression)) + geom_point() + facet_wrap(~condition) + theme_classic() + scale_colour_viridis_c(option = "viridis", name = "Mean Expression") + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7), plot.title = element_text(hjust = 0.5)) + labs(title = "Expression support for multi tool interaction overlaps", x = "Ligand or receptor gene", y = "Exact interaction", size = "Percent expressing", colour = "Mean expression")
  ggsave(file.path(output_folder, "multi_tool_overlap_expression_dotplot.png"), plot = overlap_plot, width = 16, height = 11, dpi = 300)
  print(top_overlaps[, c("condition", "source", "target", "ligand", "receptor", "standalone_count", "method_count", "minimum_percent_expressing", "mean_percent_expressing", "methods")])
}
