# ============================================================
# Script: downstream_results.R
#
# Purpose:
# Creates plots, rankings, cross tool comparisons and candidate
# expression outputs after all four tools have finished.
# ============================================================

source("config.R")
check_pipeline_config(require_seurat = TRUE)

seurat_object <- load_pipeline_seurat()
condition_info <- get_pipeline_conditions(seurat_object)
conditions <- condition_info$condition
condition_folders <- condition_info$folder
n_conditions <- length(conditions)

cat("Conditions found:", paste(conditions, collapse = ", "), "\n")

run_r_script <- function(script_file) {
  cat("\n============================================================\n")
  cat("Running:", script_file, "\n")
  cat("============================================================\n")
  source(script_file, local = new.env(parent = globalenv()))
}

run_r_script("cellchat/plot_cellchat_ligand_receptors.R")
run_r_script("cellchat/rank_cellchat.R")
run_r_script("cellphonedb/plot_cellphonedb_results.R")
run_r_script("cellphonedb/rank_cellphonedb.R")
run_r_script("liana/plot_liana_results.R")
run_r_script("liana/rank_liana.R")
run_r_script("singlecellsignalr/plot_singlecellsignalr_results.R")
run_r_script("singlecellsignalr/rank_singlecellsignalr.R")

# combine the ranked outputs into one standard table
library(dplyr)
library(tidyr)

ranked_files <- c(CellChat = file.path(results_folder, "comparison/cellchat/cellchat_ranked_interactions.csv"), CellPhoneDB = file.path(results_folder, "comparison/cellphonedb/cellphonedb_ranked_interactions.csv"), LIANA = file.path(results_folder, "comparison/liana/liana_ranked_interactions.csv"), SingleCellSignalR = file.path(results_folder, "comparison/singlecellsignalr/singlecellsignalr_ranked_interactions.csv"))

all_tools_long <- bind_rows(lapply(names(ranked_files), function(method_name) {
  x <- read.csv(ranked_files[[method_name]], stringsAsFactors = FALSE, check.names = FALSE)
  x$method <- method_name
  x
}))

all_tools_long <- all_tools_long %>% group_by(method, condition) %>% arrange(rank, .by_group = TRUE) %>% mutate(category_rank = row_number(), selected_top_two = category_rank <= 2, change_category = "Top", ranking_value = percentile) %>% ungroup()

if (n_conditions == 2) {
  condition_1 <- conditions[1]
  condition_2 <- conditions[2]
  all_tools_classified <- all_tools_long %>% select(method, interaction_id, source, target, ligand, receptor, condition, score, rank, percentile) %>% pivot_wider(names_from = condition, values_from = c(score, rank, percentile))
  score_1 <- paste0("score_", condition_1)
  score_2 <- paste0("score_", condition_2)
  percentile_1 <- paste0("percentile_", condition_1)
  percentile_2 <- paste0("percentile_", condition_2)
  all_tools_classified <- all_tools_classified %>% mutate(change_category = case_when(!is.na(.data[[score_1]]) & is.na(.data[[score_2]]) ~ "Disappears", is.na(.data[[score_1]]) & !is.na(.data[[score_2]]) ~ "Appears", TRUE ~ "Present_in_both"), score_change = .data[[score_2]] - .data[[score_1]], relative_change = score_change / pmax(abs(.data[[score_1]]), .Machine$double.eps), ranking_value = pmax(.data[[percentile_1]], .data[[percentile_2]], na.rm = TRUE)) %>% group_by(method, change_category) %>% arrange(desc(ranking_value), .by_group = TRUE) %>% mutate(category_rank = row_number(), selected_top_two = category_rank <= 2) %>% ungroup()
} else {
  all_tools_classified <- all_tools_long
}

dir.create(file.path(results_folder, "comparison"), recursive = TRUE, showWarnings = FALSE)
write.csv(all_tools_classified, file.path(results_folder, "comparison/all_tools_classified_interactions.csv"), row.names = FALSE)
write.csv(all_tools_long, file.path(results_folder, "comparison/all_tools_ranked_interactions_long.csv"), row.names = FALSE)

run_r_script("cross_tools/cross_tool_expression_analysis.R")
run_r_script("cross_tools/candidate_expression_all_tools.R")
run_r_script("cross_tools/find_multi_tool_overlaps_and_expression.R")
run_r_script("cross_tools/prioritise_candidates_and_overlaps.R")
