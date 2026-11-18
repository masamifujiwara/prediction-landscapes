# =============================================================================
# Script:  g04_shared_surprises.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Figure 6, the months in which many estuaries missed together, with
#   temperature and salinity. Four panels on one time axis: (a) biodiversity
#   and (b) abundance: proportion of estuary series with a large forecast miss
#   one month ahead in each target month, drawn above the line when the
#   observed value was higher than the forecast and below when it was lower;
#   (c) water temperature and (d) salinity: departure from each estuary's and
#   sampling method's usual value for the calendar month, averaged over the
#   estuaries and both sampling methods.
#
# Run after: b08, d01, e09, lib/figures.R
#
# Inputs:
#   01_Analysis/Data/d1_compare_model_metrics_dictionary_v3_scaled/forecasts_joined_diversity.rds (d01)
#   01_Analysis/Data/e9_compare_predictability_abundance/forecasts_joined_abundance.rds (e09)
#   01_Analysis/Data/b8_build_environment_timeseries/env_alpha_index.rds (b08)
#
# Outputs:
#   01_Analysis/figures/Figure_6_shared_surprises.tiff and .pdf
#   01_Analysis/Data/g4_shared_surprises/g4_monthly_shares.csv
#   01_Analysis/Logs/g4_shared_surprises/run_<timestamp>.txt
#   Earlier copies of the figure and the CSV are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   A large miss is an absolute one-month-ahead prediction error of SARIMA
#   (automatic) more than LARGE_SD times the series' RMSE over all scored
#   forecast dates (the f01 definition). Only estuary series are used (42
#   diversity and 112 abundance series), so the proportion is a proportion of
#   estuary series; targets in April and May 2020 (not sampled) are not scored
#   and are shaded grey in every panel. Dotted lines show the average
#   proportion over all months in each direction. Temperature and salinity
#   departures are taken from the monthly means of each estuary and sampling
#   method, using sampled months only (filled months are left out), so the
#   environmental and biological series can be read on one time axis. The
#   month in EVENTS is marked with a dashed line. The log also lists the
#   months with the largest proportions and their rank correlation with the
#   size of the environmental departures.
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
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(tidyr); library(patchwork)})
script_stem <- "g4_shared_surprises"
in_div <- file.path(PATHS$data_dir, "d1_compare_model_metrics_dictionary_v3_scaled", "forecasts_joined_diversity.rds")
in_abd <- file.path(PATHS$data_dir, "e9_compare_predictability_abundance", "forecasts_joined_abundance.rds")
in_env <- file.path(PATHS$data_dir, "b8_build_environment_timeseries", "env_alpha_index.rds")
MODEL <- "sarima_auto"; LARGE_SD <- 2; W_MM <- 170; H_MM <- 175
EVENTS <- data.frame(date = as.Date("2015-07-01"), label = "Jul 2015")
GAP <- data.frame(xmin = as.Date("2020-04-01"), xmax = as.Date("2020-06-01"))   # not sampled (settings$gap_months)
GAP_LABEL <- "After the 2020\nsampling gap"

dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(c(in_div, in_abd, in_env), c("d01_compare_model_metrics_dictionary_v3_scaled.R", "e09_compare_predictability_abundance.R",
                                             "b08_build_environment_timeseries.R"))

read_bay <- function(f, cls, id_col) {
  readRDS(f) %>% filter(model == MODEL, lead == 1L, scored, ok, scale == "bay") %>%
    transmute(class = cls, sid = paste(BAY_ID, GEAR, .data[[id_col]]), target_date, err = y_obs - mean)
}
d <- bind_rows(read_bay(in_div, "diversity", "q"), read_bay(in_abd, "abundance", "SPECIES_CODE")) %>%
  group_by(class, sid) %>% mutate(z = err / sqrt(mean(err^2))) %>% ungroup()
n_ser <- d %>% distinct(class, sid) %>% count(class)
if (!identical(setNames(n_ser$n, n_ser$class)[c("diversity", "abundance")], c(diversity = 42L, abundance = 112L))) {
  stop("Expected 42 bay diversity and 112 bay abundance series", call. = FALSE)
}
m <- d %>% group_by(class, target_date) %>%
  summarise(series = n(), higher = mean(z > LARGE_SD), lower = mean(z < -LARGE_SD), .groups = "drop")
chance <- d %>% group_by(class) %>% summarise(higher = mean(z > LARGE_SD), lower = mean(z < -LARGE_SD), .groups = "drop")
cat("Share of estuary series with a large miss in an average month:\n"); print(as.data.frame(chance), row.names = FALSE)

env <- readRDS(in_env) %>% filter(!filled) %>% group_by(BAY_ID, GEAR, MONTH) %>%
  mutate(temp = TEMP__mean - mean(TEMP__mean, na.rm = TRUE), sal = SAL__mean - mean(SAL__mean, na.rm = TRUE)) %>% ungroup() %>%
  group_by(DATE) %>% summarise(temp = mean(temp, na.rm = TRUE), sal = mean(sal, na.rm = TRUE), .groups = "drop") %>%
  filter(DATE >= min(m$target_date), DATE <= max(m$target_date))

lab_div <- "(a) Biodiversity (42 estuary series)"; lab_abd <- "(b) Abundance (112 estuary series)"
lab_t <- "(c) Water temperature"; lab_s <- "(d) Salinity"
Y_MISS <- "Proportion of estuary series\nwith a large forecast miss"
DIR_LAB <- c(higher = "Observed higher than forecast", lower = "Observed lower than forecast")
pl <- bind_rows(
  m %>% transmute(panel = ifelse(class == "diversity", lab_div, lab_abd), date = target_date,
                  higher = higher, lower = -lower) %>%
    pivot_longer(c(higher, lower), names_to = "dir", values_to = "v") %>% mutate(dir = DIR_LAB[dir]),
  env %>% transmute(panel = lab_t, date = DATE, dir = "env", v = temp),
  env %>% transmute(panel = lab_s, date = DATE, dir = "env", v = sal)) %>%
  mutate(panel = factor(panel, c(lab_div, lab_abd, lab_t, lab_s)))
dir_cols <- setNames(c(PAL[["blue"]], PAL[["orange"]]), unname(DIR_LAB))
# Average proportion over all months in each direction (dotted reference lines).
ref <- chance %>% transmute(panel = ifelse(class == "diversity", lab_div, lab_abd), higher, lower = -lower) %>%
  pivot_longer(c(higher, lower), values_to = "v")

x_sc <- scale_x_date(date_breaks = "4 years", date_labels = "%Y", limits = range(pl$date) + c(-20, 20), expand = expansion(mult = 0.01))
ev_lines <- geom_vline(data = EVENTS, aes(xintercept = date), colour = "grey60", linetype = "dashed", linewidth = 0.3)
panel_plot <- function(pn, bars, y_lab) {
  d <- filter(pl, panel == pn)
  g <- ggplot(d, aes(date, v)) +
    geom_rect(data = GAP, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey80", alpha = 0.6, inherit.aes = FALSE) +
    ev_lines + geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.25) +
    facet_wrap(~ panel) + x_sc + labs(x = NULL, y = y_lab) + theme_paper() + theme(strip.text = element_text(hjust = 0))
  if (bars) {
    g <- g + geom_linerange(aes(ymin = 0, ymax = v, colour = dir), linewidth = 0.45) +
      geom_hline(data = filter(ref, panel == pn), aes(yintercept = v), linetype = "dotted", colour = "grey20", linewidth = 0.35) +
      scale_colour_manual(values = dir_cols, breaks = unname(DIR_LAB), name = "Large forecast miss:") +
      guides(colour = guide_legend(override.aes = list(linewidth = 3))) +
      scale_y_continuous(labels = function(x) format(abs(x), trim = TRUE, drop0trailing = TRUE))  # bars below the line are shares too
  } else {
    g <- g + geom_line(colour = PAL[["grey_dark"]], linewidth = 0.3)
  }
  g
}
pa <- panel_plot(lab_div, TRUE, Y_MISS) +
  geom_text(data = EVENTS, aes(date, Inf, label = label), vjust = 1.4, hjust = 1.08, size = FONT_SIZE / .pt,
            colour = "grey30", inherit.aes = FALSE) +
  annotate("text", x = GAP$xmin, y = Inf, label = GAP_LABEL, vjust = 1.3, hjust = 1.06, size = FONT_SIZE / .pt, colour = "grey30", lineheight = 0.9)
pb <- panel_plot(lab_abd, TRUE, Y_MISS)
pc <- panel_plot(lab_t, FALSE, "Temperature\ndeparture (\u00b0C)"); pd <- panel_plot(lab_s, FALSE, "Salinity\ndeparture")
no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p <- (pa + no_x) / (pb + no_x) / (pc + no_x) / pd + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

out_csv <- file.path(dirs$data, "g4_monthly_shares.csv")
old <- list.files(dirs$data, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)
readr::write_csv(m %>% left_join(env, by = c("target_date" = "DATE")) %>% mutate(across(c(higher, lower, temp, sal), ~ round(.x, 4))), out_csv)
save_figure(p, "Figure_6_shared_surprises", W_MM, H_MM)
cat("Saved:", out_csv, "\n")
cat("\nMonths with the largest share of estuaries missing at once:\n")
print(as.data.frame(m %>% mutate(total = higher + lower) %>% group_by(class) %>% slice_max(total, n = 6) %>%
  mutate(across(c(higher, lower, total), ~ round(.x, 3)))), row.names = FALSE)
cat("\nRank correlation of the monthly share with the size of the environmental departure:\n")
print(as.data.frame(m %>% left_join(env, by = c("target_date" = "DATE")) %>% group_by(class) %>%
  summarise(temp = cor(higher + lower, abs(temp), method = "spearman", use = "complete.obs"),
            salinity = cor(higher + lower, abs(sal), method = "spearman", use = "complete.obs"), .groups = "drop")), row.names = FALSE)
end_log("finished")
cat("g04_shared_surprises.R completed\n")
