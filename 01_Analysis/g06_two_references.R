# =============================================================================
# Script:  g06_two_references.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Figure 4, skill against the two references: skill of SARIMA
#   (automatic) by months ahead against the cumulative average (all data up
#   to the forecast date) and the 12-month average (the last 12 months); both
#   references add the seasonal curve. Lines show the median over the series
#   of each class, scale and sampling method; estuary panels also show the
#   band of the middle 50%.
#
# Run after: f01, lib/figures.R
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds (f01)
#
# Outputs:
#   01_Analysis/figures/Figure_4_two_references.tiff and .pdf
#   01_Analysis/Data/g6_two_references/g6_two_references_values.csv
#   01_Analysis/Logs/g6_two_references/run_<timestamp>.txt
#   Earlier copies of the figure and the CSV are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   The two references are shown side by side, following Wesselkamp et al.
#   (2025), because the choice of reference changes what skill means: the
#   12-month average carries recent levels forward, the cumulative average
#   does not. Coastwide panels have only 3 diversity or 8 abundance series per
#   sampling method, so they show the median line without a band. The two
#   references differ by line type and symbol as well as color (cumulative
#   average solid with circles, 12-month average dashed with triangles, marks
#   at 1, 6, 12, 18 and 24 months).
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
script_stem <- "g6_two_references"
in_lead <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds")
MODEL <- "sarima_auto"; W_MM <- 170; H_MM <- 115

dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(in_lead, "f01_forecast_scoring.R")

ref_lab <- c(skill_past = "Cumulative average (all earlier months)", skill_m12 = "12-month average (last 12 months)")
ref_cols <- setNames(c(PAL[["blue"]], PAL[["orange"]]), ref_lab)
# References differ by line type and symbol as well as color (marks at these leads)
ref_lty <- setNames(c("solid", "dashed"), ref_lab)
ref_shp <- setNames(c(16, 17), ref_lab)
MARK_LEADS <- c(1, 6, 12, 18, 24)
v <- readRDS(in_lead) %>% filter(model == MODEL) %>%
  tidyr::pivot_longer(c(skill_past, skill_m12), names_to = "ref", values_to = "skill") %>%
  group_by(class, scale, GEAR, lead, ref) %>%
  summarise(series = n(), median = stats::median(skill), lo = stats::quantile(skill, 0.25), hi = stats::quantile(skill, 0.75), .groups = "drop") %>%
  mutate(Ref = factor(ref_lab[ref], ref_lab),
         panel = factor(paste0(CLASS_LAB[class], ", ", SCALE_LAB[scale], "\n", GEAR_LAB[GEAR]),
                        levels = as.vector(outer(c("Estuaries\nBag seine", "Estuaries\nBay trawl", "Coastwide\nBag seine", "Coastwide\nBay trawl"),
                                                 c("Biodiversity, ", "Abundance, "), function(a, b) paste0(b, a)))))
if (anyNA(v$panel)) stop("Unexpected class, scale or gear", call. = FALSE)

p <- ggplot(v, aes(lead, median, colour = Ref, fill = Ref, linetype = Ref, shape = Ref)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_ribbon(data = filter(v, scale == "bay"), aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(data = filter(v, lead %in% MARK_LEADS), size = 1.4) +
  facet_wrap(~ panel, nrow = 2) +
  scale_colour_manual(values = ref_cols, name = "Skill against:") + scale_fill_manual(values = ref_cols, name = "Skill against:") +
  scale_linetype_manual(values = ref_lty, name = "Skill against:") + scale_shape_manual(values = ref_shp, name = "Skill against:") +
  scale_x_continuous(breaks = c(1, 6, 12, 18, 24)) +
  labs(x = "Months ahead", y = "Skill (proportion of squared error removed)") + theme_paper()

out_csv <- file.path(dirs$data, "g6_two_references_values.csv")
old <- list.files(dirs$data, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)
readr::write_csv(v %>% select(class, scale, GEAR, lead, ref, series, median, lo, hi) %>% mutate(across(c(median, lo, hi), ~ round(.x, 4))), out_csv)
save_figure(p, "Figure_4_two_references", W_MM, H_MM)
cat("Saved:", out_csv, "\n")
print(as.data.frame(v %>% filter(lead %in% c(1, 12, 24)) %>% select(class, scale, GEAR, lead, ref, median) %>%
  tidyr::pivot_wider(names_from = ref, values_from = median) %>% mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
end_log("finished")
cat("g06_two_references.R completed\n")
