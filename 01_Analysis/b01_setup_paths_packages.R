# =============================================================================
# Script:  b01_setup_paths_packages.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Setup check for the project. Safe to run any number of times, on any
#   computer, but not needed when switching computers.
#     1. Finds the project folder and writes the marker file .project_root at
#        its top if it is missing. Every script finds the project through it.
#     2. Checks that the shared files in 01_Analysis/lib are present, builds
#        PATHS with lib/paths.R (which also creates the fixed folders) and
#        loads lib/helpers.R.
#     3. Checks that every package the scripts use is installed, and prints one
#        install.packages() line for any that are missing.
#     4. Writes the installed package versions and sessionInfo() to a log.
#
# Run after: none (first script of the pipeline; a00 is run by hand).
#
# Inputs:
#   01_Analysis/*.R and the .R files in the subfolders listed in SCRIPT_SUBDIRS
#     (scanned as text to collect package names; no data files)
#   01_Analysis/lib/paths.R, 01_Analysis/lib/helpers.R (written by hand)
#
# Outputs:
#   .project_root at the top of the project folder (only if missing)
#   The fixed folders under 01_Analysis and 02_Results (created by
#     lib/paths.R if missing)
#   01_Analysis/Logs/b1_setup_paths_packages/sessionInfo.txt (overwritten on
#     each run; the folder keeps the stem b1_setup_paths_packages)
#
# Method notes:
#   Each script builds PATHS itself when it starts (see lib/bootstrap_block.R),
#   so no file of absolute paths is saved and nothing needs to be rerun when
#   the project moves to another computer. Packages are collected from the code
#   (library(), require(), requireNamespace() and pkg::), with comments ignored,
#   so the list follows the scripts.
#
# AI assistance: This script was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

rm(list = ls())

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
# Subfolders of 01_Analysis whose scripts are scanned for packages, besides the
# top level. Archive folders (OLD, ARCHIVE_, BACKUP_) are left out on purpose.
SCRIPT_SUBDIRS <- c("lib", "drivers", "e07_arimax_env_forecasts_abundance")

# Library folder used in the install line printed on Linux. If it does not
# exist on the machine, the first folder of .libPaths() is used instead.
LINUX_LIB <- "/home/masami/R/x86_64-pc-linux-gnu-library/4.6"

# -----------------------------------------------------------------------------
# 1. Find the project folder (and write the marker if it is missing)
# -----------------------------------------------------------------------------
# Same search as the standard opening block. If no .project_root is found, the
# project folder is taken to be the parent of the 01_Analysis folder this
# script sits in (it must also contain 02_Results), and the marker is written
# there.
local({
  a <- commandArgs(FALSE)
  here <- c(unlist(lapply(rev(sys.frames()), function(f) get0("ofile", envir = f, inherits = FALSE))),
            sub("^--file=", "", grep("^--file=", a, value = TRUE)),
            a[which(a == "-f") + 1L],
            if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
              tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL),
            getwd())
  here <- gsub("~+~", " ", here[!is.na(here) & nzchar(here)], fixed = TRUE)  # Rscript on Linux writes spaces as ~+~
  start <- vapply(here, function(h) {
    d <- normalizePath(h, winslash = "/", mustWork = FALSE)
    if (!dir.exists(d)) d <- dirname(d)
    d
  }, character(1), USE.NAMES = FALSE)
  find_up <- function(d, ok) {
    repeat {
      if (ok(d)) return(d)
      if (dirname(d) == d) return(NULL)
      d <- dirname(d)
    }
  }
  for (d in start) {
    r <- find_up(d, function(x) file.exists(file.path(x, ".project_root")))
    if (!is.null(r)) {
      cat("Marker found:   ", file.path(r, ".project_root"), "\n", sep = "")
      assign(".PROJECT_ROOT", r, envir = globalenv())
      return(invisible())
    }
  }
  for (d in start) {
    r <- find_up(d, function(x) basename(x) == "01_Analysis" &&
                   dir.exists(file.path(dirname(x), "02_Results")))
    if (!is.null(r)) {
      root <- dirname(r)
      writeLines(c("Marker file for the project folder. Scripts look upward for this file to",
                   "find the project, so keep it here, at the top of the project folder.",
                   "See 01_Analysis/lib/paths.R."),
                 file.path(root, ".project_root"))
      cat("Marker WRITTEN: ", file.path(root, ".project_root"), "\n", sep = "")
      assign(".PROJECT_ROOT", root, envir = globalenv())
      return(invisible())
    }
  }
  stop("Project folder not found: no .project_root, and no 01_Analysis folder with a ",
       "02_Results folder beside it, above ", paste(unique(here), collapse = " or "))
})

# The shared files every script needs. Checked before sourcing, so a missing
# file gives a clear message rather than a "cannot open file" error.
lib_files <- file.path(.PROJECT_ROOT, "01_Analysis", "lib", c("paths.R", "helpers.R"))
if (!all(file.exists(lib_files))) {
  stop("Missing shared file(s) in 01_Analysis/lib:\n  ",
       paste(lib_files[!file.exists(lib_files)], collapse = "\n  "), call. = FALSE)
}
source(file.path(.PROJECT_ROOT, "01_Analysis", "lib", "paths.R"))
PATHS <- get_paths(.PROJECT_ROOT)
setwd(PATHS$analysis_dir)
source(file.path(PATHS$analysis_dir, "lib", "helpers.R"))
cat("Project folder: ", PATHS$project_dir, "\n", sep = "")
cat("Shared files:   lib/paths.R, lib/helpers.R found\n")

# -----------------------------------------------------------------------------
# 2. Packages used by the scripts
# -----------------------------------------------------------------------------
# Collected from the scripts themselves (library(), require(), requireNamespace()
# and pkg::), so the list follows the code. Comments are ignored. Scripts in a
# subfolder not named in SCRIPT_SUBDIRS are not scanned.
sub_dirs <- file.path(PATHS$analysis_dir, SCRIPT_SUBDIRS)
if (any(!dir.exists(sub_dirs))) {
  cat("NOTE: subfolder(s) in SCRIPT_SUBDIRS not found, skipped: ",
      paste(SCRIPT_SUBDIRS[!dir.exists(sub_dirs)], collapse = ", "), "\n", sep = "")
}
script_files <- c(
  list.files(PATHS$analysis_dir, pattern = "\\.R$", full.names = TRUE),
  unlist(lapply(sub_dirs[dir.exists(sub_dirs)], function(d)
    list.files(d, pattern = "\\.R$", full.names = TRUE)))
)
code <- unlist(lapply(script_files, function(f)
  sub("#.*$", "", readLines(f, warn = FALSE, encoding = "UTF-8"))))
grab <- function(pattern) {
  m <- regmatches(code, gregexpr(pattern, code, perl = TRUE))
  unique(sub(pattern, "\\1", unlist(m), perl = TRUE))
}
pkgs <- unique(c(
  grab("\\b(?:library|require)\\(\\s*[\"']?([A-Za-z][A-Za-z0-9.]*)"),
  grab("\\brequireNamespace\\(\\s*[\"']([A-Za-z][A-Za-z0-9.]*)[\"']"),
  grab("\\b([A-Za-z][A-Za-z0-9.]*)::")
))
base_pkgs <- rownames(installed.packages(priority = "base"))
pkgs <- sort(setdiff(pkgs, base_pkgs))
installed <- vapply(pkgs, function(p) nzchar(system.file(package = p)), logical(1))
versions  <- vapply(pkgs, function(p) if (installed[[p]]) as.character(utils::packageVersion(p)) else "MISSING",
                    character(1))

cat("\nPackages used by ", length(script_files), " scripts: ", length(pkgs), "\n", sep = "")
print(data.frame(package = pkgs, version = versions), row.names = FALSE, right = FALSE)

missing <- pkgs[!installed]
if (length(missing)) {
  vec <- paste0("c(", paste0('"', missing, '"', collapse = ", "), ")")
  linux_lib <- if (dir.exists(LINUX_LIB)) LINUX_LIB else .libPaths()[1]
  cat("\nMISSING: ", paste(missing, collapse = ", "), "\nInstall with:\n", sep = "")
  if (Sys.info()[["sysname"]] == "Linux") {
    cat('  install.packages(', vec, ', lib = "', linux_lib, '", dependencies = TRUE)\n', sep = "")
  } else {
    cat("  install.packages(", vec, ", dependencies = TRUE)\n", sep = "")
  }
} else {
  cat("\nAll packages installed.\n")
}

# -----------------------------------------------------------------------------
# 3. Log
# -----------------------------------------------------------------------------
log_dir <- file.path(PATHS$logs_dir, "b1_setup_paths_packages")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
session_file <- file.path(log_dir, "sessionInfo.txt")
writeLines(c(
  "b01_setup_paths_packages.R",
  paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Machine:  ", Sys.info()[["nodename"]]),
  paste("Project:  ", PATHS$project_dir),
  "",
  "Packages used by the scripts:",
  sprintf("  %-20s %s", pkgs, versions),
  "",
  capture.output(sessionInfo())
), session_file)
cat("\nLog written: ", session_file, "\n", sep = "")
if (length(missing)) cat("Setup INCOMPLETE: install the missing packages above.\n") else cat("Setup complete.\n")
