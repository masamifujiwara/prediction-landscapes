# =============================================================================
# Script:  check_order.R (drivers/)
# Author:  Masami Fujiwara
#
# Purpose:
#   Checks the run order: every output folder a script reads must be written
#   by a script that runs earlier. The run order is file-name order with the
#   exceptions listed in lib/run_order.R. Exits with status 1 and a list if
#   any input comes from a later script. Run by drivers/run_pipeline.R before
#   the pipeline when CHECK_ORDER is TRUE, or on its own.
#
# Run after: none (a check on the code, not on outputs).
#
# Inputs:
#   01_Analysis/*.R (the code of the top-level scripts)
#   01_Analysis/lib/run_order.R
#
# Outputs:
#   01_Analysis/Logs/check_order/check_order.csv (every read between scripts
#     that it found; overwritten on each run, not archived)
#
# Method notes:
#   This is a check on names, not a run:
#   - A script "writes" the folder named in its fixed script_stem, or in an
#     out_dir / diag_out / ts_dir <- file.path(..., "<folder>") line.
#   - A script "reads" a folder when that folder name appears as a quoted
#     string anywhere else in its code (comments ignored).
#   - Top-level scripts in 01_Analysis are put in the order the pipeline runs
#     them (lib/run_order.R).
#   b01 is skipped as a reader: it names folders only to scan for packages.
#
# AI assistance: This script was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

rm(list = ls())

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

SKIP <- c("b01_setup_paths_packages.R")

an      <- PATHS$analysis_dir
source(file.path(an, "lib", "run_order.R"))
scripts <- run_order(list.files(an, pattern = "\\.R$"))
code <- lapply(scripts, function(s) {
  x <- readLines(file.path(an, s), warn = FALSE, encoding = "UTF-8")
  x <- sub("(^|\\s)#.*$", "", x)          # drop comments (rough, good enough here)
  paste(x, collapse = "\n")
})
names(code) <- scripts

grab <- function(pattern, txt) {
  m <- regmatches(txt, gregexpr(pattern, txt, perl = TRUE))[[1]]
  unique(sub(pattern, "\\1", m, perl = TRUE))
}

# Who writes which folder
writes <- do.call(rbind, lapply(scripts, function(s) {
  f <- c(grab('script_stem\\s*<-\\s*"([^"]+)"', code[[s]]),
         grab('(?:out_dir|diag_out|ts_dir)\\s*<-\\s*file\\.path\\([^)]*"([A-Za-z0-9_]+)"\\s*\\)', code[[s]]))
  if (length(f)) data.frame(folder = unique(f), writer = s, stringsAsFactors = FALSE) else NULL
}))
dup <- writes$folder[duplicated(writes$folder)]
if (length(dup)) cat("Note: folders written by more than one script:", paste(unique(dup), collapse = ", "), "\n")

# Who reads which folder, and is the writer earlier?
reads <- do.call(rbind, lapply(setdiff(scripts, SKIP), function(s) {
  lits <- grab('"([A-Za-z0-9_]+)"', code[[s]])
  own  <- writes$folder[writes$writer == s]
  hit  <- writes[writes$folder %in% setdiff(lits, own) & writes$writer != s, ]
  if (!nrow(hit)) return(NULL)
  data.frame(reader = s, folder = hit$folder, writer = hit$writer, stringsAsFactors = FALSE)
}))
if (is.null(reads)) reads <- data.frame(reader = character(), folder = character(),
                                        writer = character(), ok = logical())
# Position in the run order, the same order the pipeline uses
reads$ok <- match(reads$writer, scripts) < match(reads$reader, scripts)

log_dir <- file.path(PATHS$logs_dir, "check_order")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(reads, file.path(log_dir, "check_order.csv"), row.names = FALSE)

cat(length(scripts), "top-level scripts;", nrow(writes), "output folders;",
    nrow(reads), "reads between scripts.\n")
bad <- reads[!reads$ok, ]
if (nrow(bad)) {
  cat("\nINPUT FROM A LATER SCRIPT:\n")
  print(bad[, c("reader", "folder", "writer")], row.names = FALSE, right = FALSE)
  if (!interactive()) quit(status = 1L)
} else {
  cat("Order OK: every input is written by a script that runs earlier.\n")
}
