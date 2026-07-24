# ============================================================
# Script: prioritise_candidates_and_overlaps.R
#
# Purpose:
# This script selects the top three expression supported
# candidates from each communication method and condition.
#
# Cross tool overlaps are analysed separately using:
# find_multi_tool_overlaps_and_expression.R
#
# Output:
# Results are saved to:
# the folder set in config.R/comparison/final_prioritisation
# ============================================================

# load pipeline settings
source("config.R")
check_pipeline_config()

# load packages
library(dplyr)

# set folders
comparison_folder <- file.path(results_folder, "comparison")
expression_folder <- file.path(comparison_folder, "candidate_expression")
output_folder <- file.path(comparison_folder, "final_prioritisation")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# set number of candidates
top_per_method <- 3

# load classified interactions and expression results
all_tools <- read.csv(file.path(comparison_folder, "all_tools_classified_interactions.csv"), stringsAsFactors = FALSE, check.names = FALSE)
expression_results <- read.csv(file.path(expression_folder, "all_tools_candidate_expression.csv"), stringsAsFactors = FALSE, check.names = FALSE)
selected_candidates <- all_tools %>% filter(selected_top_two %in% c(TRUE, "TRUE"), !is.na(ligand), !is.na(receptor), trimws(ligand) != "", trimws(receptor) != "")
expression_results <- expression_results %>% filter(!is.na(ligand), !is.na(receptor), trimws(ligand) != "", trimws(receptor) != "", !is.na(gene), trimws(gene) != "")

# summarise expression support for each candidate and condition
candidate_expression_summary <- expression_results %>% group_by(method, condition, change_category, source, target, ligand, receptor) %>% summarise(required_gene_count = n_distinct(paste(gene_role, gene)), minimum_percent_expressing = ifelse(all(is.na(percent_expressing)), NA_real_, min(percent_expressing, na.rm = TRUE)), mean_percent_expressing = ifelse(all(is.na(percent_expressing)), NA_real_, mean(percent_expressing, na.rm = TRUE)), minimum_mean_expression = ifelse(all(is.na(mean_expression)), NA_real_, min(mean_expression, na.rm = TRUE)), mean_expression = ifelse(all(is.na(mean_expression)), NA_real_, mean(mean_expression, na.rm = TRUE)), .groups = "drop")

# use the long ranked table so this works with one or more conditions
candidate_scores <- read.csv(file.path(comparison_folder, "all_tools_ranked_interactions_long.csv"), stringsAsFactors = FALSE, check.names = FALSE) %>% mutate(ranking_value = percentile)

if ("condition" %in% colnames(selected_candidates)) {
  selected_keys <- selected_candidates %>% select(method, condition, change_category, source, target, ligand, receptor) %>% distinct()
  candidate_scores <- candidate_scores %>% inner_join(selected_keys, by = c("method", "condition", "source", "target", "ligand", "receptor"))
} else {
  selected_keys <- selected_candidates %>% select(method, change_category, source, target, ligand, receptor) %>% distinct()
  candidate_scores <- candidate_scores %>% inner_join(selected_keys, by = c("method", "source", "target", "ligand", "receptor"))
}

candidate_priorities <- candidate_scores %>% left_join(candidate_expression_summary, by = c("method", "condition", "change_category", "source", "target", "ligand", "receptor")) %>% filter(!is.na(minimum_percent_expressing), !is.na(mean_percent_expressing), !is.na(mean_expression), is.finite(minimum_percent_expressing), is.finite(mean_percent_expressing), is.finite(mean_expression)) %>% group_by(method, condition) %>% arrange(desc(minimum_percent_expressing), desc(mean_percent_expressing), desc(mean_expression), desc(ranking_value), .by_group = TRUE) %>% mutate(method_priority_rank = row_number(), selected_top_three = method_priority_rank <= top_per_method) %>% ungroup()

write.csv(candidate_priorities, file.path(output_folder, "method_specific_candidate_priorities.csv"), row.names = FALSE)
top_three_per_method <- candidate_priorities %>% filter(selected_top_three)
write.csv(top_three_per_method, file.path(output_folder, "top_three_candidates_per_method.csv"), row.names = FALSE)

cat("\nNumber of complete shortlisted candidates per method and condition:\n")
print(table(candidate_priorities$method, candidate_priorities$condition))
cat("\nTop three candidates per method and condition:\n")
print(top_three_per_method[, c("method", "condition", "source", "target", "ligand", "receptor", "minimum_percent_expressing", "mean_percent_expressing", "mean_expression", "ranking_value")])
