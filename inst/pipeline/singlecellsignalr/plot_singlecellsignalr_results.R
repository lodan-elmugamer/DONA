# ============================================================
# Script: plot_singlecellsignalr_results.R
#
# Purpose:
# This script creates SingleCellSignalR plots.
#
# Reference:
# https://r-graph-gallery.com/320-the-basis-of-bubble-plot.html
#
# Output:
# PNGs and CSV tables saved to:
# the folder set in config.R/singlecellsignalr/figures
# the folder set in config.R/singlecellsignalr/csv
# ============================================================

# load pipeline settings
source("config.R")
check_pipeline_config()
condition_info <- get_pipeline_conditions()

library(ggplot2)
library(viridis)

singlecell_folder <- file.path(results_folder, "singlecellsignalr")
figures_folder <- file.path(singlecell_folder, "figures")
csv_folder <- file.path(singlecell_folder, "csv")
dir.create(figures_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_folder, recursive = TRUE, showWarnings = FALSE)

scr <- do.call(rbind, lapply(seq_len(nrow(condition_info)), function(i) {
  x <- read.csv(file.path(singlecell_folder, condition_info$folder[i], paste0(condition_info$folder[i], "_singlecellsignalr_paracrines.csv")))
  x$condition <- condition_info$condition[i]
  x
}))

# fix the labels
ligand_col <- "L"
receptor_col <- "R"
score_col <- "LR.score"

clean_label <- function(x) {
  x <- gsub("Normal_NK_fixed", "Normal_NK", x)
  x <- gsub("Normal_", "", x)
  x <- gsub("Malignant_", "", x)
  x <- gsub("_", " ", x)
  x <- gsub(" vs ", " -> ", x)
  return(x)
}

scr$ligand_receptor <- paste(scr[[ligand_col]], scr[[receptor_col]], sep = " -> ")
scr$cell_pair <- clean_label(scr$population_pair)
scr$score <- scr[[score_col]]
scr$interaction_label <- paste(scr$ligand_receptor, scr$cell_pair, sep = " | ")
write.csv(scr, file.path(csv_folder, "singlecellsignalr_combined_results.csv"), row.names = FALSE)

scr_significant <- scr[scr$pval < 0.05, ]
exact_summary <- aggregate(score ~ condition + interaction_label, data = scr_significant, FUN = mean)
top_exact <- unique(unlist(lapply(unique(exact_summary$condition), function(condition_name) head(exact_summary[exact_summary$condition == condition_name, ][order(exact_summary[exact_summary$condition == condition_name, ]$score, decreasing = TRUE), "interaction_label"], 10))))
exact_plot_data <- exact_summary[exact_summary$interaction_label %in% top_exact, ]
exact_order <- aggregate(score ~ interaction_label, data = exact_plot_data, FUN = max)
exact_order <- exact_order[order(exact_order$score, decreasing = TRUE), ]
exact_plot_data$interaction_label <- factor(exact_plot_data$interaction_label, levels = rev(exact_order$interaction_label))
write.csv(exact_plot_data, file.path(csv_folder, "singlecellsignalr_top_exact_paracrine_interactions.csv"), row.names = FALSE)

topplot <- ggplot(exact_plot_data, aes(x = condition, y = interaction_label)) + geom_point(aes(colour = score), size = 4, alpha = 0.9) + scale_colour_viridis_c(option = "viridis", name = "Mean LRscore") + theme_classic() + theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 7), plot.title = element_text(hjust = 0.5)) + labs(title = "Top SingleCellSignalR paracrine interactions", x = "Condition", y = "Ligand receptor pair | Sender receiver pair")
ggsave(file.path(figures_folder, "singlecellsignalr_top_exact_paracrine_interactions.png"), plot = topplot, width = 12, height = 10, dpi = 300)

if (nrow(condition_info) == 2) {
  change_summary <- aggregate(score ~ condition + ligand_receptor + cell_pair, data = scr_significant, FUN = mean)
  first <- change_summary[change_summary$condition == condition_info$condition[1], c("ligand_receptor", "cell_pair", "score")]
  second <- change_summary[change_summary$condition == condition_info$condition[2], c("ligand_receptor", "cell_pair", "score")]
  colnames(first)[3] <- "condition_1_score"
  colnames(second)[3] <- "condition_2_score"
  change_data <- merge(first, second, by = c("ligand_receptor", "cell_pair"), all = TRUE)
  change_data[is.na(change_data)] <- 0
  change_data$difference <- change_data$condition_2_score - change_data$condition_1_score
  change_data$abs_difference <- abs(change_data$difference)
  change_data$interaction_label <- paste(change_data$ligand_receptor, change_data$cell_pair, sep = " | ")
  top_change <- head(change_data[order(change_data$abs_difference, decreasing = TRUE), ], 20)
  top_change$interaction_label <- factor(top_change$interaction_label, levels = rev(top_change$interaction_label))
  write.csv(top_change, file.path(csv_folder, "singlecellsignalr_largest_condition_changes.csv"), row.names = FALSE)
  changeplot <- ggplot(top_change, aes(x = interaction_label, y = difference)) + geom_col() + coord_flip() + theme_classic() + theme(axis.text.y = element_text(size = 7), plot.title = element_text(hjust = 0.5)) + labs(title = "Largest SingleCellSignalR changes between conditions", x = "Ligand receptor pair | Sender receiver pair", y = paste(condition_info$condition[2], "LRscore minus", condition_info$condition[1], "LRscore"))
  ggsave(file.path(figures_folder, "singlecellsignalr_largest_condition_changes.png"), plot = changeplot, width = 12, height = 8, dpi = 300)
}
