# =============================================================================
# Script:  g08_quadrants.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Figure 7, accuracy against skill in four management situations.
#   Each series is a point: accuracy (RMSE of the out-of-sample prediction
#   errors divided by the tolerance; below 1 = within the tolerance) against
#   skill against the cumulative average (above 0 = better than the
#   cumulative average), one and 12 months ahead, for SARIMA (automatic).
#   Faint lines at accuracy 1 and skill 0 and corner labels name the four
#   management situations (forecast-led, mean-led, forecast with buffer,
#   precaution).
#
# Run after: b09, f01, lib/figures.R, lib/forecast_rules.R
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds (f01)
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09, tolerance)
#
# Outputs:
#   01_Analysis/figures/Figure_7_quadrants.tiff and .pdf
#   01_Analysis/Data/g8_quadrants/g8_quadrant_points.csv
#   01_Analysis/Logs/g8_quadrants/run_<timestamp>.txt
#   Earlier copies of the figure and the CSV are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   The lines are one reasonable summary of a continuous gradient, and the
#   caption says so. Points are drawn at their values. The quadrant counts in
#   the summary table by group (h02, Online Resource 2, Table S6) use the
#   uncertainty-range rule instead (a series counts on the favorable side only
#   when its bootstrap range clears the line), so a few points near a line
#   are counted there on the cautious side. The counts printed to the log are
#   by value. Biodiversity and abundance differ by shape as well as color
#   (circles, triangles), estuaries filled and coastwide open. The script
#   checks that the 48 diversity and 128 abundance series of the design are
#   present at each lead.
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
source(file.path(PATHS$analysis_dir, "lib", "forecast_rules.R"))
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
script_stem <- "g8_quadrants"
in_lead <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds")
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
MODEL <- "sarima_auto"; LEADS <- c(1L, 12L); W_MM <- 170; H_MM <- 95

dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(c(in_lead, in_settings), c("f01_forecast_scoring.R", "b09_forecast_settings.R"))
settings <- load_settings(in_settings)
tol <- c(diversity = settings$tolerance$diversity_rel, abundance = settings$tolerance$abundance_log)

v <- readRDS(in_lead) %>% filter(model == MODEL, lead %in% LEADS) %>%
  mutate(accuracy = rmse_scaled / tol[class], Class = factor(CLASS_LAB[class], CLASS_LAB), Scale = factor(SCALE_LAB[scale], SCALE_LAB),
         Lead = factor(ifelse(lead == 1L, "1 month ahead", "12 months ahead"), c("1 month ahead", "12 months ahead")),
         Group = factor(paste0(CLASS_LAB[class], ", ", tolower(SCALE_LAB[scale])),
                        c("Biodiversity, estuaries", "Biodiversity, coastwide", "Abundance, estuaries", "Abundance, coastwide")))
if (nrow(v) != 176L * length(LEADS)) stop("Expected 176 series at each lead", call. = FALSE)
pad <- function(r, f = 0.10) r + c(-1, 1) * f * diff(r)   # room for the corner labels outside the points
x_rng <- pad(range(c(v$accuracy, 1))); y_rng <- pad(range(c(v$skill_past, 0)))
corners <- data.frame(x = c(x_rng[1], x_rng[1], x_rng[2], x_rng[2]), y = c(y_rng[2], y_rng[1], y_rng[2], y_rng[1]),
                      hj = c(0, 0, 1, 1), vj = c(1, 0, 1, 0),
                      label = c("Forecast-led", "Mean-led", "Forecast with buffer", "Precaution"))
class_cols <- c(Biodiversity = PAL[["blue"]], Abundance = PAL[["orange"]])

p <- ggplot(v, aes(accuracy, skill_past)) +
  geom_vline(xintercept = 1, colour = "grey60", linewidth = 0.3) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = Class, shape = Group), size = 1.5, stroke = 0.6) +
  geom_label(data = corners, aes(x, y, label = label, hjust = hj, vjust = vj), size = FONT_SIZE / .pt, colour = "grey40",
             fill = "white", label.size = 0, label.padding = unit(1.5, "pt")) +
  facet_wrap(~ Lead) +
  scale_colour_manual(values = class_cols, guide = "none") +
  scale_shape_manual(values = c("Biodiversity, estuaries" = 16, "Biodiversity, coastwide" = 1, "Abundance, estuaries" = 17, "Abundance, coastwide" = 2),
                     name = NULL, guide = guide_legend(override.aes = list(colour = unname(class_cols[c("Biodiversity", "Biodiversity", "Abundance", "Abundance")])))) +
  scale_x_continuous(breaks = scales::pretty_breaks(5), limits = x_rng, expand = c(0, 0)) +
  scale_y_continuous(limits = y_rng, expand = c(0, 0)) +
  labs(x = "Prediction error relative to the tolerance (RMSE / tolerance; below 1 = within the tolerance)",
       y = "Skill against the cumulative average\n(proportion of squared error removed)") + theme_paper()

out_csv <- file.path(dirs$data, "g8_quadrant_points.csv")
old <- list.files(dirs$data, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)
readr::write_csv(v %>% select(class, scale, BAY_ID, GEAR, target, lead, rmse_scaled, accuracy, skill_past, rmse_hi, skill_past_lo) %>%
  mutate(across(where(is.double), ~ round(.x, 4))), out_csv)
save_figure(p, "Figure_7_quadrants", W_MM, H_MM)
cat("Saved:", out_csv, "\n")
cat("Series in each corner by value (Table S6 uses the uncertainty-range rule):\n")
print(as.data.frame(v %>% mutate(q = ifelse(accuracy <= 1, ifelse(skill_past > 0, "forecast-led", "mean-led"),
                                            ifelse(skill_past > 0, "forecast with buffer", "precaution"))) %>%
  count(class, scale, lead, q) %>% tidyr::pivot_wider(names_from = q, values_from = n, values_fill = 0)), row.names = FALSE)
end_log("finished")
cat("g08_quadrants.R completed\n")
