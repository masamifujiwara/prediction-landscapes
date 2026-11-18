# =============================================================================
# Script:  g05_misses_by_lead.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Figure 3, forecast misses and the percentage of forecasts within
#   the tolerance by months ahead (absolute accuracy). Top row: size of the
#   monthly misses (out-of-sample prediction errors) of SARIMA (automatic):
#   median (solid) and middle 50% (dashed) of the monthly misses pooled over
#   the series of each group, with open circles for the median miss of the
#   cumulative average (the usual value for the calendar month, from earlier
#   months only) and a labelled grey line at the tolerance. Bottom row:
#   percentage of monthly forecasts within the tolerance.
#
# Run after: b09, f01, lib/figures.R, lib/forecast_rules.R
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_miss_distribution.rds (f01)
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09, tolerance)
#
# Outputs:
#   01_Analysis/figures/Figure_3_misses_by_lead.tiff and .pdf
#   01_Analysis/Data/g5_misses_by_lead/g5_misses_by_lead_values.csv
#   01_Analysis/Logs/g5_misses_by_lead/run_<timestamp>.txt
#   Earlier copies of the figure and the CSV are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   Diversity and abundance are shown side by side on their own scales
#   because their units differ: diversity misses as a percentage of the
#   series mean, abundance misses in log units (log(CPUE + 1)). The tolerance
#   comes from settings (30% for diversity, 0.5 log units for abundance) and
#   is the same at every lead. The middle 80% of the misses is written to the
#   values file but not drawn. Sampling methods differ by symbol as well as
#   color (squares bag seine, triangles bay trawl, at 1, 6, 12, 18 and 24
#   months), so the figure reads without color.
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
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})
script_stem <- "g5_misses_by_lead"
in_miss <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_miss_distribution.rds")
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
MODEL <- "sarima_auto"; REF <- "ref_past_average"; W_MM <- 170; H_MM <- 125
REF_LEADS <- c(1, 6, 12, 18, 24)   # leads marked with the reference circles

dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(c(in_miss, in_settings), c("f01_forecast_scoring.R", "b09_forecast_settings.R"))
settings <- load_settings(in_settings)
tol <- c(diversity = 100 * settings$tolerance$diversity_rel, abundance = settings$tolerance$abundance_log)

v_all <- readRDS(in_miss) %>% filter(model %in% c(MODEL, REF)) %>%
  mutate(mult = ifelse(class == "diversity", 100, 1), across(c(q10, q25, median, q75, q90), ~ .x * mult),
         Scale = factor(SCALE_LAB[scale], SCALE_LAB), Gear = factor(GEAR_LAB[GEAR], GEAR_LAB))
v <- v_all %>% filter(model == MODEL)
ref <- v_all %>% filter(model == REF, lead %in% REF_LEADS)
if (nrow(v) != 2L * 2L * 2L * 24L) stop("Expected 2 classes x 2 scales x 2 gears x 24 leads", call. = FALSE)
if (nrow(ref) != 2L * 2L * 2L * length(REF_LEADS)) stop("Reference rows missing for ", REF, call. = FALSE)

# Sampling methods differ by symbol as well as color (marks at these leads)
GEAR_SHAPES <- c("Bag seine" = 15, "Bay trawl" = 17)
MARK_LEADS <- c(1, 6, 12, 18, 24)
block <- function(cls, y_lab) {
  d <- v %>% filter(class == cls)
  r <- ref %>% filter(class == cls)
  tol_lab <- data.frame(Scale = factor(levels(d$Scale), levels(d$Scale)), x = 24, y = tol[[cls]])
  top <- ggplot(d, aes(lead, median, colour = Gear)) +
    geom_hline(yintercept = tol[[cls]], colour = "grey55", linewidth = 0.4) +
    geom_text(data = tol_lab, aes(x, y), label = "Tolerance", inherit.aes = FALSE, hjust = 1, vjust = -0.4,
              size = FONT_SIZE / .pt, colour = "grey35") +
    geom_line(aes(y = q25), linetype = "dashed", linewidth = 0.35, show.legend = FALSE) +
    geom_line(aes(y = q75), linetype = "dashed", linewidth = 0.35, show.legend = FALSE) +
    geom_line(linewidth = 0.7) +
    geom_point(data = filter(d, lead %in% MARK_LEADS), aes(shape = Gear), size = 1.5) +
    geom_point(data = r, shape = 1, size = 1.6, stroke = 0.5, show.legend = FALSE) +
    facet_wrap(~ Scale, nrow = 1) + scale_x_continuous(breaks = c(1, 6, 12, 18, 24)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.03))) +
    scale_colour_manual(values = GEAR_COLS, name = NULL) + scale_shape_manual(values = GEAR_SHAPES, name = NULL) +
    labs(title = CLASS_LAB[[cls]], x = NULL, y = y_lab) + theme_paper() +
    theme(plot.title = element_text(size = FONT_SIZE + 1, face = "bold"), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  bot <- ggplot(d, aes(lead, hit_rate, colour = Gear, shape = Gear)) + geom_line(linewidth = 0.7) +
    geom_point(data = filter(d, lead %in% MARK_LEADS), size = 1.5) +
    facet_wrap(~ Scale, nrow = 1) + scale_x_continuous(breaks = c(1, 6, 12, 18, 24)) +
    scale_y_continuous(labels = function(x) round(100 * x), limits = c(min(v$hit_rate) - 0.03, max(v$hit_rate) + 0.03)) +
    scale_colour_manual(values = GEAR_COLS, name = NULL) + scale_shape_manual(values = GEAR_SHAPES, name = NULL) +
    labs(x = "Months ahead", y = "Forecasts within\ntolerance (%)") + theme_paper() + theme(strip.text = element_blank(), strip.background = element_blank())
  top / bot + plot_layout(heights = c(2, 1))
}
p <- (block("diversity", "Monthly miss\n(% of series mean)") | block("abundance", "Monthly miss\n(log units)")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

out_csv <- file.path(dirs$data, "g5_misses_by_lead_values.csv")
old <- list.files(dirs$data, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)
readr::write_csv(v %>% select(class, scale, GEAR, lead, n, q10, q25, median, q75, q90, hit_rate) %>%
  mutate(across(c(q10, q25, median, q75, q90, hit_rate), ~ round(.x, 4))), out_csv)
save_figure(p, "Figure_3_misses_by_lead", W_MM, H_MM)
cat("Saved:", out_csv, "\n")
print(as.data.frame(v %>% filter(lead %in% c(1, 12, 24)) %>% select(class, scale, GEAR, lead, median, q25, q75, hit_rate) %>%
  mutate(across(c(median, q25, q75, hit_rate), ~ round(.x, 2)))), row.names = FALSE)
end_log("finished")
cat("g05_misses_by_lead.R completed\n")
