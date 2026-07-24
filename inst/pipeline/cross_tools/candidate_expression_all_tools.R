# ============================================================
# Script: candidate_expression_all_tools.R
#
# Purpose:
# This script checks ligand and receptor expression for the
# selected candidates from CellChat, CellPhoneDB, LIANA and
# SingleCellSignalR.
#
# Output:
# One combined expression CSV and one dot plot per method saved
# to:
# the folder set in config.R/comparison/candidate_expression
# ============================================================

# load pipeline settings
source("config.R")
check_pipeline_config()
condition_info <- get_pipeline_conditions()
conditions <- condition_info$condition

# load packages
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# set folders
candidates_file <- file.path(results_folder, "comparison/all_tools_classified_interactions.csv")
output_folder <- file.path(results_folder, "comparison/candidate_expression")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# load selected candidates
candidates <- read.csv(candidates_file, stringsAsFactors = FALSE, check.names = FALSE)
candidates <- candidates[candidates$selected_top_two %in% c(TRUE, "TRUE"), ]
candidate_columns <- intersect(c("method", "condition", "change_category", "category_rank", "source", "target", "ligand", "receptor"), colnames(candidates))
candidates <- candidates[, candidate_columns]
if (!"condition" %in% colnames(candidates)) candidates$condition <- NA_character_
candidates <- candidates %>% filter(!is.na(ligand), !is.na(receptor), trimws(ligand) != "", trimws(receptor) != "")

# convert displayed NK name to metadata name
candidates$source_group <- ifelse(candidates$source == "Normal_NK", "Normal_NK_fixed", candidates$source)
candidates$target_group <- ifelse(candidates$target == "Normal_NK", "Normal_NK_fixed", candidates$target)

# load Seurat object
seurat_object <- load_pipeline_seurat()
seurat_object <- JoinLayers(seurat_object, assay = "RNA")
DefaultAssay(seurat_object) <- "RNA"

# split ligand and receptor complexes into individual genes
split_genes <- function(x) {
  x <- x[!is.na(x) & x != ""]
  unique(trimws(unlist(strsplit(x, "_", fixed = TRUE))))
}

candidate_genes <- unique(c(split_genes(candidates$ligand), split_genes(candidates$receptor)))
missing_genes <- setdiff(candidate_genes, rownames(seurat_object))
if (length(missing_genes) > 0) cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
candidate_genes <- intersect(candidate_genes, rownames(seurat_object))
if (length(candidate_genes) == 0) stop("None of the candidate ligand or receptor genes were found in the Seurat object.")

vars_to_fetch <- c(candidate_genes, celltype_column)
if (!is.null(condition_column)) vars_to_fetch <- c(vars_to_fetch, condition_column)
expression_data <- FetchData(seurat_object, vars = vars_to_fetch, layer = "data")
if (is.null(condition_column)) expression_data$.pipeline_condition <- conditions[1] else expression_data$.pipeline_condition <- as.character(expression_data[[condition_column]])

# calculate expression for one gene
summarise_gene <- function(gene_name, cell_group, condition_name) {
  cells <- expression_data[expression_data[[celltype_column]] == cell_group & expression_data$.pipeline_condition == condition_name, , drop = FALSE]
  if (nrow(cells) == 0) return(data.frame(n_cells = 0, mean_expression = NA_real_, percent_expressing = NA_real_))
  data.frame(n_cells = nrow(cells), mean_expression = mean(cells[[gene_name]], na.rm = TRUE), percent_expressing = mean(cells[[gene_name]] > 0, na.rm = TRUE) * 100)
}

# calculate expression for one candidate
summarise_candidate <- function(candidate) {
  ligand_genes <- intersect(split_genes(candidate$ligand), candidate_genes)
  receptor_genes <- intersect(split_genes(candidate$receptor), candidate_genes)
  conditions_to_check <- if (!is.na(candidate$condition) && candidate$condition != "") candidate$condition else conditions
  bind_rows(lapply(conditions_to_check, function(condition_name) {
    ligand_results <- bind_rows(lapply(ligand_genes, function(gene_name) summarise_gene(gene_name, candidate$source_group, condition_name) %>% mutate(gene = gene_name, gene_role = "Ligand", cell_group = candidate$source)))
    receptor_results <- bind_rows(lapply(receptor_genes, function(gene_name) summarise_gene(gene_name, candidate$target_group, condition_name) %>% mutate(gene = gene_name, gene_role = "Receptor", cell_group = candidate$target)))
    bind_rows(ligand_results, receptor_results) %>% mutate(method = candidate$method, change_category = candidate$change_category, category_rank = candidate$category_rank, source = candidate$source, target = candidate$target, ligand = candidate$ligand, receptor = candidate$receptor, condition = condition_name)
  }))
}

expression_results <- bind_rows(lapply(seq_len(nrow(candidates)), function(i) summarise_candidate(candidates[i, ])))
expression_results <- expression_results %>% mutate(mean_expression = round(mean_expression, 3), percent_expressing = round(percent_expressing, 1), condition = factor(condition, levels = conditions))

clean_group <- function(x) {
  x <- gsub("Normal_NK_fixed", "Normal_NK", x)
  x <- gsub("Normal_", "", x)
  x <- gsub("Malignant_", "", x)
  gsub("_", " ", x)
}

expression_results <- expression_results %>% mutate(source_plot = clean_group(source), target_plot = clean_group(target), interaction_plot = paste0(ligand, " -> ", receptor, " | ", source_plot, " -> ", target_plot), gene_plot = paste(gene_role, gene, sep = ": "))
expression_results$interaction_plot <- factor(expression_results$interaction_plot, levels = rev(unique(expression_results$interaction_plot)))
write.csv(expression_results, file.path(output_folder, "all_tools_candidate_expression.csv"), row.names = FALSE)

for (method_name in unique(candidates$method)) {
  method_data <- expression_results[expression_results$method == method_name, ]
  expression_plot <- ggplot(method_data, aes(x = gene_plot, y = interaction_plot, size = percent_expressing, colour = mean_expression)) + geom_point() + scale_colour_viridis_c(option = "viridis", name = "Mean Expression") + facet_wrap(~condition) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7), plot.title = element_text(hjust = 0.5)) + labs(title = paste(method_name, "candidate ligand and receptor expression"), x = "Ligand or receptor gene", y = "Selected interaction", size = "Percent expressing", colour = "Mean expression")
  method_file <- tolower(gsub(" ", "_", method_name))
  ggsave(file.path(output_folder, paste0(method_file, "_candidate_expression_dotplot.png")), plot = expression_plot, width = 15, height = 10, dpi = 300)
}

print(table(candidates$method, candidates$change_category))
cat("Selected candidates:", nrow(candidates), "\n")
cat("Expression rows:", nrow(expression_results), "\n")
