#!/bin/bash

# ============================================================
# Script: run_cellphonedb.sh
#
# Purpose:
# Runs CellPhoneDB for every condition folder exported by
# export_cellphonedb_inputs.R.
# ============================================================

PROJECT_FOLDER="${PROJECT_FOLDER:-$PWD}"
RESULTS_FOLDER="${RESULTS_FOLDER:-$PROJECT_FOLDER/results/cellphonedb}"
CPDB_ZIP="${CPDB_ZIP:-$PROJECT_FOLDER/resources/cellphonedb.zip}"
VENV_ACTIVATE="${VENV_ACTIVATE:-$PROJECT_FOLDER/cellphonedb/venv/bin/activate}"
THREADS="${THREADS:-32}"

cd "$PROJECT_FOLDER"
source "$VENV_ACTIVATE"

found_input=false
for INPUT_FOLDER in "$RESULTS_FOLDER"/*/input; do
  [ -d "$INPUT_FOLDER" ] || continue
  found_input=true
  CONDITION_FOLDER="$(basename "$(dirname "$INPUT_FOLDER")")"
  OUTPUT_FOLDER="$RESULTS_FOLDER/$CONDITION_FOLDER/output"
  mkdir -p "$OUTPUT_FOLDER"
  echo "Running CellPhoneDB for: $CONDITION_FOLDER"

  python - <<PY
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_statistical_analysis_method.call(
    cpdb_file_path="$CPDB_ZIP",
    meta_file_path="$INPUT_FOLDER/metadata.txt",
    counts_file_path="$INPUT_FOLDER",
    counts_data="hgnc_symbol",
    output_path="$OUTPUT_FOLDER",
    threads=int("$THREADS"),
    debug_seed=123
)
PY
done

if [ "$found_input" = false ]; then
  echo "No CellPhoneDB input folders were found under: $RESULTS_FOLDER"
  exit 1
fi
