# ============================================================
# Script: rank_singlecellsignalr.R
#
# Purpose:
# This script extracts and ranks focused SingleCellSignalR interactions
# within each condition. When exactly two conditions are present,
# it also compares their presence between conditions.
#
# Output:
# Results are saved to:
# the folder set in config.R/comparison/singlecellsignalr
# ============================================================

# load pipeline settings
source("config.R")
check_pipeline_config()
condition_info <- get_pipeline_conditions()
conditions <- condition_info$condition
n_conditions <- length(conditions)

# load packages
library(dplyr)
library(tidyr)
library(ggplot2)


# set folders
scsr_file <- file.path(results_folder, "singlecellsignalr/csv/singlecellsignalr_combined_results.csv")
comparison_folder <- file.path(results_folder, "comparison/singlecellsignalr")
dir.create(comparison_folder, recursive = TRUE, showWarnings = FALSE)

input_data <- read.csv(scsr_file, stringsAsFactors = FALSE, check.names = FALSE)

# extract and rank one condition
prepare_data <- function(input_table, condition_name) {
  input_table %>% filter(condition == condition_name, pval < 0.05) %>% separate(population_pair, into = c("source", "target"), sep = "_vs_", remove = FALSE) %>% mutate(ligand = L, receptor = R, score = as.numeric(LR.score), pvalue = as.numeric(pval)) %>% filter(keep_focused_pairs(source, target)) %>% group_by(source, target, ligand, receptor) %>% summarise(score = max(score, na.rm = TRUE), pvalue = min(pvalue, na.rm = TRUE), .groups = "drop") %>% mutate(condition = condition_name, interaction_id = paste(source, target, ligand, receptor, sep = " | "), cell_pair = paste(source, target, sep = " -> "), ligand_receptor = paste(ligand, receptor, sep = " -> ")) %>% arrange(desc(score)) %>% mutate(rank = row_number(), percentile = percent_rank(score))
}

# prepare every condition
ranked_results <- bind_rows(lapply(conditions, function(condition_name) prepare_data(input_data, condition_name)))

# save ranked results
write.csv(ranked_results, file.path(comparison_folder, "singlecellsignalr_ranked_interactions.csv"), row.names = FALSE)

# compare conditions only when exactly two are present
if (n_conditions == 2) {
  condition_1 <- conditions[1]
  condition_2 <- conditions[2]
  comparison <- ranked_results %>% select(interaction_id, source, target, ligand, receptor, condition, score, pvalue, rank, percentile) %>% pivot_wider(names_from = condition, values_from = c(score, pvalue, rank, percentile)) %>% mutate(behaviour = case_when(!is.na(.data[[paste0("score_", condition_1)]]) & !is.na(.data[[paste0("score_", condition_2)]]) ~ "Present_in_both", !is.na(.data[[paste0("score_", condition_1)]]) & is.na(.data[[paste0("score_", condition_2)]]) ~ "Disappears", is.na(.data[[paste0("score_", condition_1)]]) & !is.na(.data[[paste0("score_", condition_2)]]) ~ "Appears"))
  write.csv(comparison, file.path(comparison_folder, "singlecellsignalr_condition_comparison.csv"), row.names = FALSE)
}

# select top 20 interactions per condition
top_interactions <- ranked_results %>% group_by(condition) %>% slice_max(score, n = 20, with_ties = FALSE) %>% ungroup() %>% mutate(cell_pair_plot = gsub("Normal_NK_fixed", "Normal_NK", cell_pair), plot_label = paste(ligand_receptor, cell_pair_plot, sep = " | "))

# plot top interactions
interaction_plot <- ggplot(top_interactions, aes(x = reorder(plot_label, score), y = score)) + geom_col() + coord_flip() + facet_wrap(~condition, scales = "free_y") + theme_classic() + theme(axis.text.y = element_text(size = 7), plot.title = element_text(hjust = 0.5)) + labs(title = "Top SingleCellSignalR interactions within each condition", x = "Ligand receptor pair | Sender receiver pair", y = "SingleCellSignalR LRscore")

ggsave(file.path(comparison_folder, "singlecellsignalr_top_ranked_interactions.png"), plot = interaction_plot, width = 14, height = 10, dpi = 300)

# print interaction counts
print(table(ranked_results$condition))
