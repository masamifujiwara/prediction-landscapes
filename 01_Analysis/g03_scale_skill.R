# =============================================================================
# Script:  g03_scale_skill.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Figure 5, estuaries against the coast: skill of SARIMA (automatic)
#   against the cumulative average by months ahead (1 to 24), for the estuary
#   series and the coastwide series, by target class and sampling method.
#   Estuaries: median over the estuary series and the band of their middle
#   50%. Coastwide diversity: each of the three series (q = 0, 1, 2) drawn and
#   labelled. Coastwide abundance: median over the 8 taxa (thick line) and
#   each taxon as a thin line.
#
# Run after: f01, lib/figures.R
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds (f01)
#
# Outputs:
#   01_Analysis/figures/Figure_5_scale_skill.tiff and .pdf
#   01_Analysis/Data/g3_scale_skill/g3_scale_skill_values.csv
#   01_Analysis/Logs/g3_scale_skill/run_<timestamp>.txt
#   Earlier copies of the figure and the CSV are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   Skill = 1 - (squared prediction error of the forecast / that of the
#   cumulative average), out of sample, on the common forecast dates (f01).
#   The three coastwide diversity series are drawn one by one because a median
#   of three would hide one of them. Estuaries and coastwide differ by line
#   type and symbol as well as color (estuaries solid with circles, coastwide
#   dashed with triangles, marks at 1, 12 and 24 months). Coastwide abundance
#   taxa whose skill at 24 months ahead is above 0.3 are labelled. The script
#   checks that the 48 diversity and 128 abundance series of the design are
#   present.
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
script_stem <- "g3_scale_skill"
in_lead <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds")
MODEL <- "sarima_auto"; W_MM <- 170; H_MM <- 120

dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(in_lead, "f01_forecast_scoring.R")

bl <- readRDS(in_lead) %>% filter(model == MODEL)
if (length(unique(paste(bl$class, bl$scale, bl$BAY_ID, bl$GEAR, bl$target))) != 176L) stop("Expected 176 series", call. = FALSE)
labs3 <- function(d) d %>% mutate(Scale = factor(SCALE_LAB[scale], SCALE_LAB), Class = factor(CLASS_LAB[class], CLASS_LAB),
                                   Gear = factor(GEAR_LAB[GEAR], GEAR_LAB))
v <- bl %>% group_by(class, scale, GEAR, lead) %>%
  summarise(series = n(), median = stats::median(skill_past), lo = stats::quantile(skill_past, 0.25),
            hi = stats::quantile(skill_past, 0.75), .groups = "drop") %>% labs3()
# Coastwide: few series (3 or 8 per gear), so each is drawn as a thin line
# instead of a band.
cw <- bl %>% filter(scale == "coastwide") %>% labs3()
y_rng <- range(c(v$lo[v$scale == "bay"], v$hi[v$scale == "bay"], cw$skill_past, 0))
# Coastwide diversity: the three orders as labelled lines; label positions
# spread so that labels at most LABEL_GAP apart do not overlap.
LABEL_GAP <- 0.065 * diff(y_rng)
cw_div <- cw %>% filter(class == "diversity") %>% mutate(q_lab = paste0("q = ", sub("^q", "", target)))
spread_y <- function(y) { o <- order(y); z <- y[o]; for (i in seq_along(z)[-1]) z[i] <- max(z[i], z[i - 1] + LABEL_GAP); y[o] <- z; y }
div_lab <- cw_div %>% filter(lead == max(lead)) %>% group_by(Class, Gear) %>% mutate(y_lab = spread_y(skill_past)) %>% ungroup()
v_line <- v %>% filter(!(scale == "coastwide" & class == "diversity"))
# Scales differ by line type and symbol as well as color: estuaries solid with
# circles, coastwide dashed with triangles (marks at 1, 12 and 24 months).
SCALE_LTY <- c("Estuaries" = "solid", "Coastwide" = "dashed")
SCALE_SHP <- c("Estuaries" = 16, "Coastwide" = 17)
MARK_LEADS <- c(1, 12, 24)

p <- ggplot(v, aes(lead, median, colour = Scale, fill = Scale, linetype = Scale, shape = Scale)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_ribbon(data = function(d) filter(d, scale == "bay"), aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(data = filter(cw, class == "abundance"), aes(lead, skill_past, group = target), colour = PAL[["orange"]], alpha = 0.35, linewidth = 0.25, linetype = "dashed", inherit.aes = FALSE) +
  geom_line(data = cw_div, aes(lead, skill_past, group = target), colour = PAL[["orange"]], linewidth = 0.45, linetype = "dashed", inherit.aes = FALSE) +
  geom_point(data = filter(cw_div, lead %in% MARK_LEADS), aes(lead, skill_past), colour = PAL[["orange"]], shape = 17, size = 1.2, inherit.aes = FALSE) +
  geom_text(data = div_lab, aes(max(lead) + 0.6, y_lab, label = q_lab), inherit.aes = FALSE, hjust = 0, size = FONT_SIZE / .pt, colour = PAL[["orange"]]) +
  geom_line(data = v_line, linewidth = 0.6) + geom_point(data = filter(v_line, lead %in% MARK_LEADS), size = 1.2) +
  geom_text(data = filter(cw, class == "abundance", lead == 24, skill_past > 0.3), aes(lead, skill_past, label = target), inherit.aes = FALSE,
            hjust = 1, vjust = -0.6, size = FONT_SIZE / .pt, colour = "grey25") +
  facet_grid(Class ~ Gear) +
  scale_colour_manual(values = SCALE_COLS, name = NULL) + scale_fill_manual(values = SCALE_COLS, name = NULL) +
  scale_linetype_manual(values = SCALE_LTY, name = NULL) + scale_shape_manual(values = SCALE_SHP, name = NULL) +
  scale_x_continuous(breaks = c(1, 6, 12, 18, 24), expand = expansion(mult = c(0.02, 0.11))) +
  coord_cartesian(ylim = y_rng) +
  labs(x = "Months ahead", y = "Skill against the cumulative average") +
  theme_paper()

out_csv <- file.path(dirs$data, "g3_scale_skill_values.csv")
old <- list.files(dirs$data, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)
readr::write_csv(v %>% select(class, scale, GEAR, lead, series, median, lo, hi) %>% mutate(across(c(median, lo, hi), ~ round(.x, 4))), out_csv)
save_figure(p, "Figure_5_scale_skill", W_MM, H_MM)
cat("Saved:", out_csv, "\n")
print(as.data.frame(v %>% filter(lead %in% c(1, 3, 12, 24)) %>% select(class, scale, GEAR, lead, series, median, lo, hi) %>%
  mutate(across(c(median, lo, hi), ~ round(.x, 3)))), row.names = FALSE)
end_log("finished")
cat("g03_scale_skill.R completed\n")
