# ============================================================
# Script: run_pipeline.R
#
# Purpose:
# Runs the four ligand receptor tools and their downstream scripts
# in the same order as the original analysis.
#
# CellPhoneDB itself uses Python and is run separately with:
# bash cellphonedb/run_cellphonedb.sh
# ============================================================

source("config.R")
check_pipeline_config(require_seurat = TRUE)

run_r_script <- function(script_file) {
  cat("\n============================================================\n")
  cat("Running:", script_file, "\n")
  cat("============================================================\n")
  source(script_file, local = new.env(parent = globalenv()))
}

# inference and input export
run_r_script("cellchat/run_cellchat.R")
run_r_script("cellphonedb/export_cellphonedb_inputs.R")
run_r_script("liana/run_liana.R")
run_r_script("singlecellsignalr/run_singlecellsignalr.R")

cat("\nCellChat, LIANA and SingleCellSignalR are complete.\n")
cat("Create virtual environment inside CellPhoneDb folder and run with bash ./run_cellphonedb.sh before running downstream_results.R.\n")
