# ============================================================
# LR pipeline configuration
#
# Edit the values below before running the scripts.
# Run scripts from the main LR_pipeline_refactored folder.
# ============================================================

# input Seurat object
seurat_file <- "path/to/your_seurat_object.rds"

# main output folder
output_folder <- "results"
results_folder <- file.path(output_folder, "ligand_receptor")

# metadata columns
celltype_column <- "analysis_group_fixed"

# condition settings
# use NULL when the object contains only one condition or has no condition column
# when NULL, the whole object is analysed once and saved under single_condition
condition_column <- NULL

# optional condition values to run
# use NULL to automatically use every value found in condition_column
condition_names <- NULL

# optional cell groups to retain
# use NULL to keep every cell group
groups_of_interest <- NULL

# optional groups used for focused malignant-normal ranking
# use NULL for both to rank all sender-receiver pairs
malignant_focus <- NULL
normal_focus <- NULL

# reproducibility
pipeline_seed <- 123

# make a safe folder name
make_condition_folder <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) x <- "condition"
  return(x)
}

# check configuration before starting
check_pipeline_config <- function(require_seurat = FALSE) {
  if (require_seurat && (!file.exists(seurat_file) || dir.exists(seurat_file))) {
    stop("Update seurat_file in config.R so that it points to a valid Seurat RDS file.")
  }
  if (!nzchar(output_folder)) stop("output_folder cannot be empty.")
  if (!nzchar(celltype_column)) stop("celltype_column cannot be empty.")
}

# load and prepare the input object in one place
load_pipeline_seurat <- function() {
  check_pipeline_config(require_seurat = TRUE)
  seurat_object <- readRDS(seurat_file)
  required_columns <- celltype_column
  if (!is.null(condition_column)) required_columns <- c(required_columns, condition_column)
  missing_columns <- setdiff(required_columns, colnames(seurat_object@meta.data))
  if (length(missing_columns) > 0) stop("Missing metadata column(s): ", paste(missing_columns, collapse = ", "))
  if (!is.null(groups_of_interest)) {
    keep_cells <- rownames(seurat_object@meta.data)[seurat_object@meta.data[[celltype_column]] %in% groups_of_interest]
    seurat_object <- subset(seurat_object, cells = keep_cells)
  }
  return(seurat_object)
}

# split an object into one or more configured conditions
split_pipeline_conditions <- function(seurat_object) {
  if (is.null(condition_column)) {
    return(list(single_condition = seurat_object))
  }
  values_found <- unique(as.character(seurat_object@meta.data[[condition_column]]))
  values_found <- values_found[!is.na(values_found) & nzchar(values_found)]
  conditions_to_run <- condition_names
  if (is.null(conditions_to_run)) conditions_to_run <- values_found
  missing_conditions <- setdiff(conditions_to_run, values_found)
  if (length(missing_conditions) > 0) stop("Condition value(s) not found: ", paste(missing_conditions, collapse = ", "))
  condition_objects <- lapply(conditions_to_run, function(condition_name) {
    keep_cells <- rownames(seurat_object@meta.data)[as.character(seurat_object@meta.data[[condition_column]]) == condition_name]
    subset(seurat_object, cells = keep_cells)
  })
  names(condition_objects) <- conditions_to_run
  return(condition_objects)
}

# return condition names and their output folder names
get_pipeline_conditions <- function(seurat_object = NULL) {
  if (is.null(seurat_object)) seurat_object <- load_pipeline_seurat()
  condition_objects <- split_pipeline_conditions(seurat_object)
  condition_labels <- names(condition_objects)
  condition_folders <- vapply(condition_labels, make_condition_folder, character(1))
  if (anyDuplicated(condition_folders)) stop("Condition names create duplicate folder names. Rename the condition values in condition_names.")
  return(data.frame(condition = condition_labels, folder = condition_folders, stringsAsFactors = FALSE))
}

# focused sender-receiver filter shared by ranking scripts
keep_focused_pairs <- function(source, target) {
  if (is.null(malignant_focus) || is.null(normal_focus)) return(rep(TRUE, length(source)))
  return((source %in% malignant_focus & target %in% normal_focus) | (source %in% normal_focus & target %in% malignant_focus))
}
