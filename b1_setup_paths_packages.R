# =============================================================================
# Script: b1_setup_paths_packages.R
#
# Purpose (project directory convention):
#   - Set working directory to this script's folder (01_Analysis)
#   - Define <Project Directory> as the parent of 01_Analysis
#   - Ensure required fixed folders exist:
#       01_Analysis/{Logs,Data}
#       02_Results/{01_Figures,02_Tables,03_Diagnostics,04_Models}
#   - Create script-specific result subfolders:
#       02_Results/<bucket>/<script_stem>/
#   - Define PATHS (standardized directory list)
#   - Save PATHS to 01_Analysis/PATHS.RData (authoritative copy for all scripts)
#   - Write sessionInfo() to 01_Analysis/Logs/<script_stem>/sessionInfo.txt
#
# =============================================================================
#
# The code for future use (top of every script after b1)
# =============================================================================
## Determine script directory and setwd()
# if (requireNamespace("rstudioapi", quietly = TRUE) &&
#     rstudioapi::isAvailable() &&
#     !is.null(rstudioapi::getActiveDocumentContext()$path) &&
#     nzchar(rstudioapi::getActiveDocumentContext()$path)) {
#   script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
# } else if (!is.null(sys.frames()[[1]]$ofile)) {
#   script_dir <- dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/"))
# } else {
#   stop("Cannot determine script directory. Are you running this interactively?")
# }
# setwd(script_dir)
#
## Find project root and load PATHS
# get_project_root <- function() {
#   here <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
#   normalizePath(file.path(here, ".."), winslash = "/", mustWork = TRUE)
# }
# project_dir <- get_project_root()
# paths_file <- file.path(project_dir, "01_Analysis", "PATHS.RData")
# stopifnot(file.exists(paths_file))
# load(paths_file)   # loads PATHS
# =============================================================================

rm(list = ls())

# ─────────────────────────────────────────────────────────────────────────────
# 0. Set working directory to this script's folder
# ─────────────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  
  # 1) RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) {
      return(dirname(normalizePath(p)))
    }
  }
  
  # 2) Rscript (works with nohup)
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) {
      return(dirname(normalizePath(p)))
    }
  }
  
  # 3) Fallback: current working directory
  # (Safe for nohup if you cd into 01_Analysis first)
  return(normalizePath(getwd()))
}

script_dir <- get_script_dir()
setwd(script_dir)
cat("📁 Working directory:", normalizePath(getwd(), winslash = "/"), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Global options (session-scoped)
# ─────────────────────────────────────────────────────────────────────────────
options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Utilities: script name + project root from current folder
# ─────────────────────────────────────────────────────────────────────────────
get_script_name <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) return(basename(sub("^--file=", "", file_arg)))
  
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(basename(p))
  }
  
  "interactive_session.R"
}

get_project_root <- function() {
  here <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  normalizePath(file.path(here, ".."), winslash = "/", mustWork = TRUE)
}

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

init_project_paths <- function(allow_create_results_subfolders = TRUE) {
  project_dir <- get_project_root()
  
  analysis_dir <- file.path(project_dir, "01_Analysis")
  logs_dir     <- file.path(analysis_dir, "Logs")
  data_dir     <- file.path(analysis_dir, "Data")
  
  results_dir  <- file.path(project_dir, "02_Results")
  fig_root     <- file.path(results_dir, "01_Figures")
  tab_root     <- file.path(results_dir, "02_Tables")
  diag_root    <- file.path(results_dir, "03_Diagnostics")
  model_root   <- file.path(results_dir, "04_Models")
  
  # Create fixed structure only
  dir_create(analysis_dir)
  dir_create(logs_dir)
  dir_create(data_dir)
  
  dir_create(results_dir)
  dir_create(fig_root)
  dir_create(tab_root)
  dir_create(diag_root)
  dir_create(model_root)
  
  script_file <- get_script_name()
  script_stem <- sub("\\.[Rr]$", "", script_file)
  
  # Allowed new folders: only within results buckets by script name
  fig_dir   <- file.path(fig_root, script_stem)
  tab_dir   <- file.path(tab_root, script_stem)
  diag_dir  <- file.path(diag_root, script_stem)
  model_dir <- file.path(model_root, script_stem)
  
  if (isTRUE(allow_create_results_subfolders)) {
    dir_create(fig_dir)
    dir_create(tab_dir)
    dir_create(diag_dir)
    dir_create(model_dir)
  }
  
  list(
    project_dir = project_dir,
    
    analysis_dir = analysis_dir,
    logs_dir     = logs_dir,
    data_dir     = data_dir,
    
    results_dir  = results_dir,
    fig_root     = fig_root,
    tab_root     = tab_root,
    diag_root    = diag_root,
    model_root   = model_root,
    
    script_file  = script_file,
    script_stem  = script_stem,
    
    fig_dir      = fig_dir,
    tab_dir      = tab_dir,
    diag_dir     = diag_dir,
    model_dir    = model_dir,
    
    # canonical locations
    paths_rdata      = file.path(analysis_dir, "PATHS.RData"),
    clean_data_rdata = file.path(analysis_dir, "Data.RData")
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. Initialize PATHS and save authoritative copy to 01_Analysis/
# ─────────────────────────────────────────────────────────────────────────────
PATHS <- init_project_paths()

cat("🧭 Project directory:", PATHS$project_dir, "\n")

# Hard check for b1: it should be run from 01_Analysis after setwd(script_dir)
wd_norm <- normalizePath(getwd(), winslash = "/")
analysis_norm <- normalizePath(PATHS$analysis_dir, winslash = "/")
stopifnot(identical(wd_norm, analysis_norm))

# Save authoritative PATHS for reuse in all subsequent scripts
save(PATHS, file = PATHS$paths_rdata)
cat("✅ Saved PATHS to:", PATHS$paths_rdata, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Diagnostics: sessionInfo (01_Analysis/Logs/<script_stem>/)
# ─────────────────────────────────────────────────────────────────────────────
log_dir <- file.path(PATHS$logs_dir, PATHS$script_stem)
dir_create(log_dir)

session_file <- file.path(log_dir, "sessionInfo.txt")
sink(session_file)
cat(PATHS$script_file, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
print(sessionInfo())
sink()

cat("✅ sessionInfo written to:", session_file, "\n")
cat("✅ ", PATHS$script_file, " completed successfully\n", sep = "")