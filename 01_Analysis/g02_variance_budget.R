# =============================================================================
# Script:  g02_variance_budget.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Figure 2, the variance partition: where the month-to-month variation
#   of each series comes from. For each target class (biodiversity,
#   abundance), scale (estuaries, coastwide) and sampling method, the
#   variation is split into four parts computed in f01 for SARIMA (automatic):
#   seasonal cycle, slow change (12 to 24 months ahead), short-lived carry-over
#   (one month ahead) and unpredictable. Bars show the mean share over the
#   series of each group, so the four parts add to 100%.
#
# Run after: f01, lib/figures.R
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_variance_budget.rds (f01)
#
# Outputs:
#   01_Analysis/figures/Figure_2_variance_budget.tiff and .pdf
#   01_Analysis/Data/g2_variance_budget/g2_variance_budget_values.csv
#   01_Analysis/Logs/g2_variance_budget/run_<timestamp>.txt
#   Earlier copies of the figure and the CSV are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   SARIMA (automatic) is the base model for the partition. The bars are means
#   over series, not medians, because the shares of one series add to one, so
#   their means also add to one (medians would not). The script checks that
#   the budget holds the 48 diversity and 128 abundance series of the design.
#   Percentages are printed on segments of 5% or more.
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
source(file.path(PATHS$analysis_dir, "lib", "helpers.R"))
source(file.path(PATHS$analysis_dir, "lib", "figures.R"))
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
script_stem <- "g2_variance_budget"
in_budget <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_variance_budget.rds")
MODEL <- "sarima_auto"; W_MM <- 170; H_MM <- 95

dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(in_budget, "f01_forecast_scoring.R")

b <- readRDS(in_budget) %>% filter(model == MODEL)
if (nrow(b) != 176L) stop(nrow(b), " series in the budget; expected 176", call. = FALSE)
parts <- c(seasonal = "Seasonal cycle", slow = "Slow change (12 to 24 months)", short_lived = "Short-lived carry-over (1 month)",
           unpredictable = "Unpredictable")
if (max(abs(rowSums(b[, names(parts)]) - 1)) > 1e-8) stop("Budget shares do not add to one", call. = FALSE)
v <- b %>% group_by(class, scale, GEAR) %>% summarise(series = n(), across(all_of(names(parts)), mean), .groups = "drop") %>%
  tidyr::pivot_longer(all_of(names(parts)), names_to = "part", values_to = "share") %>%
  mutate(part = factor(parts[part], rev(parts)),
         class = factor(CLASS_LAB[class], CLASS_LAB), group = paste0(SCALE_LAB[scale], ", ", tolower(GEAR_LAB[GEAR])),
         group = factor(group, rev(c("Estuaries, bag seine", "Estuaries, bay trawl", "Coastwide, bag seine", "Coastwide, bay trawl"))))
cols <- setNames(c(PAL[["orange"]], PAL[["blue"]], PAL[["aqua"]], PAL[["grey_light"]]), parts)

p <- ggplot(v, aes(share, group, fill = part)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.4) +
  geom_text(data = function(d) filter(d, share >= 0.05), aes(label = sprintf("%.0f%%", 100 * share)),
            position = position_stack(vjust = 0.5), size = FONT_SIZE / .pt, colour = "grey10") +
  facet_wrap(~ class, ncol = 1) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = cols, breaks = unname(parts), name = NULL) +
  labs(x = "Percentage of the variation of the monthly series", y = NULL) +
  theme_paper() + theme(panel.grid.major.y = element_blank(), plot.margin = margin(3, 10, 2, 2))

out_csv <- file.path(dirs$data, "g2_variance_budget_values.csv")
old <- list.files(dirs$data, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)
readr::write_csv(v %>% mutate(share = round(share, 4)), out_csv)
save_figure(p, "Figure_2_variance_budget", W_MM, H_MM)
cat("Saved:", out_csv, "\n")
print(as.data.frame(v %>% mutate(share = round(share, 3))), row.names = FALSE)
end_log("finished")
cat("g02_variance_budget.R completed\n")
