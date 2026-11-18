# =============================================================================
# File:    lib/bootstrap_block.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Reference copy of the opening block that every pipeline script starts
#   with. It is not sourced: each script carries its own copy, because a
#   script has to find the project folder before it can source anything from
#   lib/. If this block changes, change every script to match.
#
# Used by: every top-level script in 01_Analysis and both drivers (each holds
#   its own copy of the block below).
#
# Functions:
#   None. The block sets .PROJECT_ROOT in the global environment, sources
#   lib/paths.R, builds PATHS with get_paths() and sets the working folder to
#   01_Analysis.
#
# Method notes:
#   The project folder is found by looking upward for the marker file
#   .project_root. Nothing machine-specific is stored, so the project can be
#   copied or synced to any computer and run without a setup step.
#   Tested launch methods: RStudio Source, source() (also nested), Rscript
#   with an absolute or relative path, nohup Rscript, R -f, R CMD BATCH,
#   Rscript -e, and code pasted into the console from inside the project
#   folder.
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

# ---- Bootstrap: find the project folder and build PATHS (same in every script)
# Looks for the marker file .project_root above, in order: the file being
# source()d (RStudio Source, source()), the file given to Rscript / R -f /
# R CMD BATCH (also under nohup), the file open in RStudio, the working folder.
local({
  a <- commandArgs(FALSE)
  here <- c(unlist(lapply(rev(sys.frames()), function(f) get0("ofile", envir = f, inherits = FALSE))),
            sub("^--file=", "", grep("^--file=", a, value = TRUE)),
            a[which(a == "-f") + 1L],
            if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
              tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL),
            getwd())
  here <- gsub("~+~", " ", here[!is.na(here) & nzchar(here)], fixed = TRUE)  # Rscript on Linux writes spaces as ~+~
  for (h in here) {
    d <- normalizePath(h, winslash = "/", mustWork = FALSE)
    if (!dir.exists(d)) d <- dirname(d)
    repeat {
      if (file.exists(file.path(d, ".project_root"))) {
        assign(".PROJECT_ROOT", d, envir = globalenv())
        return(invisible())
      }
      if (dirname(d) == d) break
      d <- dirname(d)
    }
  }
  stop("Project folder not found: no .project_root above ", paste(unique(here), collapse = " or "))
})
source(file.path(.PROJECT_ROOT, "01_Analysis", "lib", "paths.R"))
PATHS <- get_paths(.PROJECT_ROOT)
setwd(PATHS$analysis_dir)
# -----------------------------------------------------------------------------
