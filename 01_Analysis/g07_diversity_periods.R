# =============================================================================
# Script:  g07_diversity_periods.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Online Resource 3, Figure S7: biodiversity and abundance by period.
#   Skill of SARIMA (automatic) for the estuary series, by the period of the
#   target month (1992 to 2007, 2008 to 2019, 2020 to 2024), one and 12
#   months ahead, against the cumulative average (filled symbols) and the
#   12-month average (open symbols): median over the estuary series of each
#   class and sampling method, with the range of their middle 50%.
#
# Run after: f01, lib/figures.R
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_periods.rds (f01)
#
# Outputs:
#   01_Analysis/figures/Figure_S7_diversity_periods.tiff and .pdf (name set by
#     FIGURE_LABEL)
#   01_Analysis/Data/g7_diversity_periods/g7_diversity_periods_values.csv
#   01_Analysis/Logs/g7_diversity_periods/run_<timestamp>.txt
#   Earlier copies of the figure and the CSV are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   Only estuary series are drawn (42 diversity and 112 abundance series);
#   the coastwide values (3 and 8 series per sampling method) are written to
#   the values file. Hill orders are pooled, not drawn separately. Skill is
#   shown against both references because skill against a long-term average
#   can be inflated by a lasting shift in level, while a reference that
#   follows recent levels is not; splitting by period shows whether skill
#   rests on one part of the record. Period comparisons use all forecast
#   dates (f01). Class is shown by color and shape (circles biodiversity,
#   triangles abundance) and the reference by fill, so color is never the
#   only cue.
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
script_stem <- "g7_diversity_periods"
in_per <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_periods.rds")
MODEL <- "sarima_auto"; LEADS <- c(1L, 12L); W_MM <- 170; H_MM <- 105
FIGURE_LABEL <- "Figure_S7"   # figure number in Online Resource 3

dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(in_per, "f01_forecast_scoring.R")

ref_lab <- c(skill_past = "Cumulative average", skill_m12 = "12-month average")
all_v <- readRDS(in_per) %>% filter(model == MODEL, lead %in% LEADS) %>%
  tidyr::pivot_longer(c(skill_past, skill_m12), names_to = "ref", values_to = "skill") %>%
  group_by(class, scale, GEAR, period, lead, ref) %>%
  summarise(series = n(), months = stats::median(n), median = stats::median(skill, na.rm = TRUE),
            lo = stats::quantile(skill, 0.25, na.rm = TRUE), hi = stats::quantile(skill, 0.75, na.rm = TRUE), .groups = "drop")
if (length(unique(all_v$period)) != 3L) stop("Expected three periods", call. = FALSE)
class_cols <- c(Biodiversity = PAL[["blue"]], Abundance = PAL[["orange"]])
v <- all_v %>% filter(scale == "bay") %>%
  mutate(Class = factor(CLASS_LAB[class], CLASS_LAB), Gear = factor(GEAR_LAB[GEAR], GEAR_LAB),
         Lead = factor(ifelse(lead == 1L, "1 month ahead", "12 months ahead"), c("1 month ahead", "12 months ahead")),
         Period = gsub("-", " to ", period), Ref = factor(ref_lab[ref], ref_lab),
         grp = interaction(Class, Ref, lex.order = TRUE),
         fill_key = paste(Class, ref, sep = "_"))
dodge <- position_dodge(width = 0.7)

# Class by color and shape (circle diversity, triangle abundance); reference by fill
# (filled cumulative average, open 12-month average), so color is never the only cue.
FILL <- c(Biodiversity_skill_past = class_cols[["Biodiversity"]], Biodiversity_skill_m12 = "white",
          Abundance_skill_past = class_cols[["Abundance"]], Abundance_skill_m12 = "white")
p <- ggplot(v, aes(Period, median, colour = Class, shape = Class, fill = fill_key, group = grp)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_linerange(aes(ymin = lo, ymax = hi), position = dodge, linewidth = 0.5, show.legend = FALSE) +
  geom_point(position = dodge, size = 2, stroke = 0.6) +
  facet_grid(Lead ~ Gear) + scale_colour_manual(values = class_cols, name = NULL) +
  scale_shape_manual(values = c(Biodiversity = 21, Abundance = 24), name = NULL) +
  scale_fill_manual(values = FILL, breaks = c("Biodiversity_skill_past", "Biodiversity_skill_m12"), labels = unname(ref_lab),
                    name = "Skill against:", guide = guide_legend(override.aes = list(shape = 21, colour = "grey25",
                                                                                      fill = c("grey25", "white")))) +
  guides(colour = guide_legend(override.aes = list(fill = unname(class_cols)), order = 1),
         shape = guide_legend(override.aes = list(fill = unname(class_cols)), order = 1)) +
  labs(x = "Period of the target month", y = "Skill (proportion of squared error removed)") + theme_paper() +
  theme(legend.box = "vertical")

out_csv <- file.path(dirs$data, "g7_diversity_periods_values.csv")
old <- list.files(dirs$data, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)
readr::write_csv(all_v %>% mutate(across(c(median, lo, hi), ~ round(.x, 4))), out_csv)
save_figure(p, paste0(FIGURE_LABEL, "_diversity_periods"), W_MM, H_MM)
cat("Saved:", out_csv, "\n")
print(as.data.frame(all_v %>% mutate(across(c(median, lo, hi), ~ round(.x, 3)))), row.names = FALSE)
end_log("finished")
cat("g07_diversity_periods.R completed\n")
