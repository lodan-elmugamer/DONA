# ============================================================
# Script: plot_cellphonedb_results.R
#
# Purpose:
# This script creates CellPhoneDB plots from CellPhoneDB output files.
#
# Output:
# PNGs and CSV tables saved to:
# the folder set in config.R/cellphonedb/figures
# the folder set in config.R/cellphonedb/csv
# ============================================================

# load pipeline settings
source("config.R")
check_pipeline_config()
condition_info <- get_pipeline_conditions()

library(ggplot2)
library(viridis)

cpdb_folder <- file.path(results_folder, "cellphonedb")
figures_folder <- file.path(cpdb_folder, "figures")
csv_folder <- file.path(cpdb_folder, "csv")
dir.create(figures_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_folder, recursive = TRUE, showWarnings = FALSE)

# SPECIFIC TO MY DATA
clean_label <- function(x) {
  x <- gsub("Normal_NK_fixed", "Normal_NK", x)
  x <- gsub("Normal_", "", x)
  x <- gsub("Malignant_", "", x)
  x <- gsub("_", " ", x)
  return(x)
}

find_latest_file <- function(folder, pattern) {
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop(paste("No file matching", pattern, "was found in", folder))
  files[which.max(file.info(files)$mtime)]
}

make_display_name <- function(gene, partner) {
  gene <- as.character(gene)
  partner <- sub("^simple:|^complex:", "", as.character(partner))
  ifelse(!is.na(gene) & gene != "", gene, partner)
}

make_cpdb_long <- function(condition, output_folder) {
  means <- read.delim(find_latest_file(output_folder, "statistical_analysis_means.*\\.txt$"), check.names = FALSE, stringsAsFactors = FALSE)
  pvalues <- read.delim(find_latest_file(output_folder, "statistical_analysis_pvalues.*\\.txt$"), check.names = FALSE, stringsAsFactors = FALSE)
  cell_pair_cols <- grep("\\|", colnames(means), value = TRUE)
  if (length(cell_pair_cols) == 0) stop("No sender receiver columns containing | were found")
  required_cols <- c("gene_a", "gene_b", "partner_a", "partner_b", "interacting_pair")
  missing_cols <- setdiff(required_cols, colnames(means))
  if (length(missing_cols) > 0) stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  out <- do.call(rbind, lapply(cell_pair_cols, function(cell_pair) data.frame(condition = condition, ligand = means$gene_a, receptor = means$gene_b, partner_a = means$partner_a, partner_b = means$partner_b, interacting_pair = means$interacting_pair, cell_pair_raw = cell_pair, score = means[[cell_pair]], pvalue = pvalues[[cell_pair]], stringsAsFactors = FALSE)))
  out <- out[!is.na(out$score) & !is.na(out$pvalue) & out$pvalue < 0.05 & out$score > 0, ]
  out$ligand_display <- make_display_name(out$ligand, out$partner_a)
  out$receptor_display <- make_display_name(out$receptor, out$partner_b)
  out
}

cpdb <- do.call(rbind, lapply(seq_len(nrow(condition_info)), function(i) make_cpdb_long(condition_info$condition[i], file.path(cpdb_folder, condition_info$folder[i], "output"))))
cpdb$source <- sub("\\|.*", "", cpdb$cell_pair_raw)
cpdb$target <- sub(".*\\|", "", cpdb$cell_pair_raw)
cpdb$source_label <- clean_label(cpdb$source)
cpdb$target_label <- clean_label(cpdb$target)
cpdb$ligand_receptor <- paste(cpdb$ligand_display, cpdb$receptor_display, sep = " -> ")
cpdb$cell_pair <- paste(cpdb$source_label, cpdb$target_label, sep = " -> ")
write.csv(cpdb, file.path(csv_folder, "cellphonedb_combined_significant_results.csv"), row.names = FALSE)

plot_top_dotplot <- function(data, variable, y_label, title, output_name) {
  summary_data <- aggregate(score ~ condition + data[[variable]], data = data, FUN = mean)
  colnames(summary_data)[2] <- "group"
  top_groups <- aggregate(score ~ group, data = summary_data, FUN = max)
  top_groups <- top_groups[order(top_groups$score, decreasing = TRUE), ]$group[1:min(10, nrow(top_groups))]
  plot_data <- summary_data[summary_data$group %in% top_groups, ]
  plot_data$group <- factor(plot_data$group, levels = rev(top_groups))
  write.csv(plot_data, file.path(csv_folder, paste0(output_name, ".csv")), row.names = FALSE)
  p <- ggplot(plot_data, aes(x = condition, y = group)) + geom_point(aes(colour = score), size = 4, alpha = 0.9) + scale_colour_viridis_c(option = "viridis", name = "Interaction score") + theme_classic() + theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 9), plot.title = element_text(hjust = 0.5)) + labs(title = title, x = "Condition", y = y_label, colour = "Mean score")
  ggsave(file.path(figures_folder, paste0(output_name, ".png")), plot = p, width = 8, height = 7, dpi = 300)
}

plot_top_dotplot(cpdb, "ligand_receptor", "Ligand receptor pair", "Top CellPhoneDB ligand receptor pairs", "cellphonedb_top_ligand_receptor_pairs")
plot_top_dotplot(cpdb, "cell_pair", "Sender receiver pair", "Top CellPhoneDB sender receiver pairs", "cellphonedb_top_sender_receiver_pairs")

if (nrow(condition_info) == 2) {
  change_summary <- aggregate(score ~ condition + ligand_receptor + cell_pair, data = cpdb, FUN = mean)
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
  write.csv(top_change, file.path(csv_folder, "cellphonedb_largest_condition_changes.csv"), row.names = FALSE)
  p3 <- ggplot(top_change, aes(x = interaction_label, y = difference)) + geom_col() + coord_flip() + theme_classic() + theme(axis.text.y = element_text(size = 7), plot.title = element_text(hjust = 0.5)) + labs(title = "Largest CellPhoneDB changes between conditions", x = "Ligand receptor pair | Sender receiver pair", y = paste(condition_info$condition[2], "score minus", condition_info$condition[1], "score"))
  ggsave(file.path(figures_folder, "cellphonedb_largest_condition_changes.png"), plot = p3, width = 12, height = 8, dpi = 300)
}
