# =============================================================================
# File:    lib/paths.R
# Author:  Masami Fujiwara
#
# Purpose:
#   get_paths(root) builds PATHS at run time from the project folder that the
#   opening block of each script finds (see lib/bootstrap_block.R). PATHS
#   holds the project, analysis, data, log and result folders used by every
#   script.
#
# Used by: every top-level script and both drivers, through the opening block.
#
# Functions:
#   PATHS_LAYOUT_VERSION   version number of the PATHS layout
#   get_paths()            returns PATHS (a named list of 17 entries) and
#                          creates the fixed folders if they are missing
#
# Outputs:
#   Creates, if missing: 01_Analysis, 01_Analysis/Logs, 01_Analysis/Data,
#   02_Results and its subfolders 01_Figures, 02_Tables, 03_Diagnostics and
#   04_Models, and a b1_setup_paths_packages folder inside each of those four.
#
# Method notes:
#   The entries script_file, script_stem, fig_dir, tab_dir, diag_dir and
#   model_dir name the folders of the setup script under the fixed stem
#   b1_setup_paths_packages. They are kept so that PATHS always has the same
#   entries, names and order; the pipeline scripts build their own output
#   folders from their own stems. paths_rdata and clean_data_rdata point to
#   01_Analysis/PATHS.RData and 01_Analysis/Data.RData.
#   Because the fixed folders are created here, a fresh copy of the project
#   needs no setup step.
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

PATHS_LAYOUT_VERSION <- 1L

get_paths <- function(root) {
  project_dir <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(project_dir, ".project_root"))) {
    stop("Not a project folder (no .project_root): ", project_dir)
  }

  analysis_dir <- file.path(project_dir, "01_Analysis")
  logs_dir     <- file.path(analysis_dir, "Logs")
  data_dir     <- file.path(analysis_dir, "Data")

  results_dir  <- file.path(project_dir, "02_Results")
  fig_root     <- file.path(results_dir, "01_Figures")
  tab_root     <- file.path(results_dir, "02_Tables")
  diag_root    <- file.path(results_dir, "03_Diagnostics")
  model_root   <- file.path(results_dir, "04_Models")

  # Fixed entries of the setup script, kept so that PATHS keeps the same entries (see header).
  script_file <- "b1_setup_paths_packages.R"
  script_stem <- "b1_setup_paths_packages"
  fig_dir     <- file.path(fig_root, script_stem)
  tab_dir     <- file.path(tab_root, script_stem)
  diag_dir    <- file.path(diag_root, script_stem)
  model_dir   <- file.path(model_root, script_stem)

  for (d in c(analysis_dir, logs_dir, data_dir, results_dir, fig_root, tab_root,
              diag_root, model_root, fig_dir, tab_dir, diag_dir, model_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  list(
    project_dir      = project_dir,
    analysis_dir     = analysis_dir,
    logs_dir         = logs_dir,
    data_dir         = data_dir,
    results_dir      = results_dir,
    fig_root         = fig_root,
    tab_root         = tab_root,
    diag_root        = diag_root,
    model_root       = model_root,
    script_file      = script_file,
    script_stem      = script_stem,
    fig_dir          = fig_dir,
    tab_dir          = tab_dir,
    diag_dir         = diag_dir,
    model_dir        = model_dir,
    paths_rdata      = file.path(analysis_dir, "PATHS.RData"),
    clean_data_rdata = file.path(analysis_dir, "Data.RData")
  )
}
