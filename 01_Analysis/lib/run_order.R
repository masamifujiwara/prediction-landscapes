# =============================================================================
# File:    lib/run_order.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Run order of the top-level scripts. The run order is file-name order
#   (byte order, so "_" sorts before letters) with the exceptions in
#   RUN_BEFORE: each script named on the left runs immediately before the
#   script named on the right.
#
# Used by: drivers/run_pipeline.R, drivers/check_order.R
#
# Functions:
#   RUN_BEFORE    named vector of the exceptions to file-name order
#   run_order()   sorts a vector of script file names into run order
#
# Method notes:
#   f02_slow_pattern_diagnostics.R runs before f01_forecast_scoring.R,
#   because f01 reads the series diagnostics of f02.
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

RUN_BEFORE <- c("f02_slow_pattern_diagnostics.R" = "f01_forecast_scoring.R")

run_order <- function(scripts) {
  scripts <- sort(scripts, method = "radix")
  for (s in names(RUN_BEFORE)) {
    target <- RUN_BEFORE[[s]]
    if (s %in% scripts && target %in% scripts) {
      scripts <- setdiff(scripts, s)
      scripts <- append(scripts, s, after = match(target, scripts) - 1L)
    }
  }
  scripts
}
