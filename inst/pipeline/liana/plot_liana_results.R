# ============================================================
# Script: plot_liana_results.R
#
# Purpose:
# This script creates LIANA plots using the aggregate_rank
# consensus score from liana_aggregate() (sca + cellphonedb).
#
# btw: liana_aggregate() column names can varyby
# LIANA version.
#
# Output:
# PNGs and CSV table saved to:
# the folder set in config.R/liana/figures
# the folder set in config.R/liana/csv
# ============================================================

# load pipeline settings
source("config.R")
check_pipeline_config()
condition_info <- get_pipeline_conditions()

library(ggplot2)
library(viridis)

liana_results_folder <- file.path(results_folder, "liana")
liana_figures_folder <- file.path(liana_results_folder, "figures")
liana_csv_folder <- file.path(liana_results_folder, "csv")
dir.create(liana_figures_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(liana_csv_folder, recursive = TRUE, showWarnings = FALSE)

liana <- do.call(rbind, lapply(seq_len(nrow(condition_info)), function(i) {
  x <- read.csv(file.path(liana_results_folder, condition_info$folder[i], paste0(condition_info$folder[i], "_liana_results.csv")))
  x$condition <- condition_info$condition[i]
  x
}))

required_cols <- c("source", "target", "ligand.complex", "receptor.complex", "aggregate_rank")
if (!all(required_cols %in% colnames(liana))) stop("Column names differ from expected. Run colnames(liana) and update this script to match.")
if (all(c("sca.LRscore", "cellphonedb.pvalue") %in% colnames(liana))) liana <- liana[!is.na(liana$sca.LRscore) & !is.na(liana$aggregate_rank) & liana$cellphonedb.pvalue < 0.05, ] else liana <- liana[!is.na(liana$aggregate_rank), ]

clean_label <- function(x) {
  x <- gsub("Normal_NK_fixed", "Normal_NK", x)
  x <- gsub("Normal_", "", x)
  x <- gsub("Malignant_", "", x)
  x <- gsub("_", " ", x)
  return(x)
}

liana$source_label <- clean_label(liana$source)
liana$target_label <- clean_label(liana$target)
liana$ligand_receptor <- paste(liana$ligand.complex, liana$receptor.complex, sep = " -> ")
liana$cell_pair <- paste(liana$source_label, liana$target_label, sep = " -> ")
write.csv(liana, file.path(liana_csv_folder, "liana_combined_filtered_results.csv"), row.names = FALSE)

plot_top_dotplot <- function(data, variable, y_label, title, output_name) {
  summary_data <- aggregate(aggregate_rank ~ condition + data[[variable]], data = data, FUN = mean)
  colnames(summary_data) <- c("condition", "group", "aggregate_rank")
  top_groups <- aggregate(aggregate_rank ~ group, data = summary_data, FUN = min)
  top_groups <- top_groups[order(top_groups$aggregate_rank), ]$group[1:min(10, nrow(top_groups))]
  plot_data <- summary_data[summary_data$group %in% top_groups, ]
  plot_data$group <- factor(plot_data$group, levels = rev(top_groups))
  write.csv(plot_data, file.path(liana_csv_folder, paste0(output_name, ".csv")), row.names = FALSE)
  p <- ggplot(plot_data, aes(x = condition, y = group)) + geom_point(aes(colour = aggregate_rank), size = 4, alpha = 0.9) + scale_colour_viridis_c(option = "viridis", direction = -1, name = "Aggregate rank") + theme_classic() + theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 9), plot.title = element_text(hjust = 0.5)) + labs(title = title, x = "Condition", y = y_label)
  ggsave(file.path(liana_figures_folder, paste0(output_name, ".png")), plot = p, width = 8, height = 7, dpi = 300)
}

plot_top_dotplot(liana, "ligand_receptor", "Ligand receptor pair", "Top LIANA ligand receptor pairs", "liana_top_ligand_receptor_pairs")
plot_top_dotplot(liana, "cell_pair", "Sender receiver pair", "Top LIANA sender receiver pairs", "liana_top_sender_receiver_pairs")

source_target_summary <- aggregate(aggregate_rank ~ condition + ligand_receptor + cell_pair, data = liana, FUN = mean)
top_interactions <- unique(unlist(lapply(unique(source_target_summary$condition), function(condition_name) head(source_target_summary[source_target_summary$condition == condition_name, ][order(source_target_summary[source_target_summary$condition == condition_name, ]$aggregate_rank), "ligand_receptor"], 20))))
plot_data <- source_target_summary[source_target_summary$ligand_receptor %in% top_interactions, ]
p3 <- ggplot(plot_data, aes(x = cell_pair, y = ligand_receptor, colour = aggregate_rank)) + geom_point(size = 4) + facet_wrap(~condition) + scale_colour_viridis_c(option = "viridis", direction = -1) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) + labs(title = "LIANA source target interactions", x = "Sender receiver pair", y = "Ligand receptor pair")
ggsave(file.path(liana_figures_folder, "liana_source_target_dotplot.png"), plot = p3, width = 16, height = 9, dpi = 300)

if (nrow(condition_info) == 2) {
  difference_data <- aggregate(aggregate_rank ~ condition + ligand_receptor + cell_pair, data = liana, FUN = mean)
  first <- difference_data[difference_data$condition == condition_info$condition[1], c("ligand_receptor", "cell_pair", "aggregate_rank")]
  second <- difference_data[difference_data$condition == condition_info$condition[2], c("ligand_receptor", "cell_pair", "aggregate_rank")]
  colnames(first)[3] <- "condition_1_rank"
  colnames(second)[3] <- "condition_2_rank"
  change_data <- merge(first, second, by = c("ligand_receptor", "cell_pair"), all = TRUE)
  change_data$difference <- change_data$condition_2_rank - change_data$condition_1_rank
  change_data$abs_difference <- abs(change_data$difference)
  change_data$interaction_label <- paste(change_data$ligand_receptor, change_data$cell_pair, sep = " | ")
  top_change <- head(change_data[order(change_data$abs_difference, decreasing = TRUE), ], 20)
  top_change$interaction_label <- factor(top_change$interaction_label, levels = rev(top_change$interaction_label))
  write.csv(top_change, file.path(liana_csv_folder, "liana_largest_condition_differences.csv"), row.names = FALSE)
  p4 <- ggplot(top_change, aes(x = interaction_label, y = difference)) + geom_col() + coord_flip() + theme_classic() + theme(axis.text.y = element_text(size = 7), plot.title = element_text(hjust = 0.5)) + labs(title = "Largest LIANA rank changes between conditions", x = "Ligand receptor pair | Sender receiver pair", y = paste(condition_info$condition[2], "rank minus", condition_info$condition[1], "rank"))
  ggsave(file.path(liana_figures_folder, "liana_largest_condition_differences.png"), plot = p4, width = 12, height = 8, dpi = 300)
}
