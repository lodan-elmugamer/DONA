#' Find the installed DONA pipeline
#'
#' Returns the folder containing the pipeline scripts included with DONA.
#'
#' @return A character string containing the pipeline folder.
#' @export
dona_pipeline_path <- function() {
  path <- system.file("pipeline", package = "DONA")
  if (!nzchar(path)) stop("The DONA pipeline files could not be found.")
  path
}

#' Copy the DONA pipeline to a working folder
#'
#' Copies the editable configuration and analysis scripts to a folder chosen by the user.
#' This is the easiest way to start a new analysis without editing the installed package.
#'
#' @param path Folder where the pipeline should be copied.
#' @param overwrite Whether existing files may be replaced.
#'
#' @return The copied pipeline folder, invisibly.
#' @export
dona_copy_pipeline <- function(path = "DONA_analysis", overwrite = FALSE) {
  source_path <- dona_pipeline_path()
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)) > 0 && !overwrite) stop("The destination is not empty. Use overwrite = TRUE to replace matching files.")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(source_path, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  relative <- substring(files, nchar(source_path) + 2)
  for (i in seq_along(files)) {
    destination <- file.path(path, relative[i])
    if (dir.exists(files[i])) {
      dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    } else {
      dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
      file.copy(files[i], destination, overwrite = overwrite)
    }
  }
  message("DONA pipeline copied to: ", normalizePath(path, mustWork = FALSE))
  invisible(normalizePath(path, mustWork = FALSE))
}

.run_dona_script <- function(project_folder, script) {
  project_folder <- normalizePath(project_folder, mustWork = TRUE)
  script_path <- file.path(project_folder, script)
  if (!file.exists(script_path)) stop("Could not find: ", script_path)
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(project_folder)
  sys.source(script_path, envir = new.env(parent = globalenv()))
  invisible(TRUE)
}

#' Run DONA inference
#'
#' Runs the main R inference workflow using the config.R file inside the selected project folder.
#' CellPhoneDB is run separately because it uses Python.
#'
#' @param project_folder Folder created by [dona_copy_pipeline()].
#'
#' @return TRUE invisibly after completion.
#' @export
dona_run <- function(project_folder = ".") {
  .run_dona_script(project_folder, "run_pipeline.R")
}

#' Run DONA downstream analysis
#'
#' Runs ranking, plotting, expression checks and cross-tool consensus. Condition-change
#' comparisons are only performed when exactly two conditions are available.
#'
#' @param project_folder Folder created by [dona_copy_pipeline()].
#'
#' @return TRUE invisibly after completion.
#' @export
dona_run_downstream <- function(project_folder = ".") {
  .run_dona_script(project_folder, "downstream_results.R")
}

#' Run CellPhoneDB from DONA
#'
#' Starts the included CellPhoneDB shell script. This requires Bash, Python,
#' a CellPhoneDB environment and database paths configured for the user's machine.
#'
#' @param project_folder Folder created by [dona_copy_pipeline()].
#' @param env Optional named character vector of environment variables, such as
#'   RESULTS_FOLDER, VENV_ACTIVATE and CPDB_ZIP.
#'
#' @return The shell command exit status, invisibly.
#' @export
dona_run_cellphonedb <- function(project_folder = ".", env = character()) {
  project_folder <- normalizePath(project_folder, mustWork = TRUE)
  script <- file.path(project_folder, "cellphonedb", "run_cellphonedb.sh")
  if (!file.exists(script)) stop("Could not find: ", script)
  if (.Platform$OS.type == "windows") stop("CellPhoneDB must be run through Bash, WSL or a Linux system.")
  old <- Sys.getenv(names(env), unset = NA_character_)
  on.exit({
  for (name in names(env)) if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, stats::setNames(list(old[[name]]), name))
  }, add = TRUE)
  if (length(env)) for (name in names(env)) do.call(Sys.setenv, stats::setNames(list(env[[name]]), name))
  status <- system2("bash", script, stdout = "", stderr = "")
  if (!identical(status, 0L)) stop("CellPhoneDB ended with exit status ", status, ".")
  invisible(status)
}
