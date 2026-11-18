# =============================================================================
# Script:  g12_slow_patterns.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Draws Fig. S3 (Online Resource 3): slow, non-seasonal patterns in the
#   departures from the seasonal curve, from the whole record (no forecasts).
#   Columns: diversity (left) and abundance (right). Rows: (a, b) correlation
#   of the departures with themselves 1 to 36 months apart, median by group;
#   (c, d) proportion of the variance of the departures in a 25-month running
#   mean, against the 95th percentile of the same proportion in AR(1) series
#   with the same lag-1 correlation; (e, f) standardized yearly mean departure,
#   median by group. Also writes the group summaries as CSV files.
#
# Run after: f02
#
# Inputs:
#   01_Analysis/Data/f2_slow_pattern_diagnostics/f2_series_diagnostics.rds,
#     f2_acf_by_lag.rds, f2_yearly_means.rds (f02)
#
# Outputs:
#   01_Analysis/figures/Figure_S3_slow_patterns.tiff and .pdf
#   01_Analysis/Data/g12_slow_patterns/g12_acf_by_group.csv, g12_slow_share.csv,
#     g12_yearly_by_group.csv
#   01_Analysis/Logs/g12_slow_patterns/run_<timestamp>.txt
#   Earlier copies of the figure and CSV files are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   - Description only: departures from the whole-record seasonal curve, with
#     April and May 2020 left out (as in f02).
#   - Yearly means are in standard deviations of each series' departures, so
#     that Hill orders and taxa share one axis.
#   - Coastwide diversity: the three orders are drawn as separate lines in (a)
#     and (e), labeled at the right end, and labeled in (c), because a median
#     of three would hide one series.
#   - Groups are told apart by more than color: gear by color and symbol
#     (circle bag seine, triangle bay trawl), scale by line type and symbol
#     fill (filled estuaries, open coastwide). One legend for all panels.
#   - The counts of series above the AR(1) line are checked against the numbers
#     quoted in the Fig. S3 caption (g20); the script stops if they differ.
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
need <- c("dplyr", "ggplot2", "patchwork", "readr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})

FIGURE_LABEL <- "Figure_S3"   # Online Resource 3 figure number (Fig. S3)
script_stem <- "g12_slow_patterns"
F2 <- file.path(PATHS$data_dir, "f2_slow_pattern_diagnostics")
inputs <- c(diag = file.path(F2, "f2_series_diagnostics.rds"), acf = file.path(F2, "f2_acf_by_lag.rds"),
            yearly = file.path(F2, "f2_yearly_means.rds"))
W_MM <- 170; H_MM <- 190
dirs <- script_dirs(script_stem, c("data", "logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(inputs, rep("f02_slow_pattern_diagnostics.R", 3))
out <- file.path(dirs$data, c("g12_acf_by_group.csv", "g12_slow_share.csv", "g12_yearly_by_group.csv"))
old <- out[file.exists(out)]; if (length(old)) archive_files(old)

keys <- c("class", "scale", "BAY_ID", "GEAR", "target")
dg <- readRDS(inputs[["diag"]]); ac <- readRDS(inputs[["acf"]]); ym <- readRDS(inputs[["yearly"]])
if (nrow(dg) != 176L) stop("Expected 176 series in f02", call. = FALSE)
GROUPS <- c("Bag seine, estuaries", "Bag seine, coastwide", "Bay trawl, estuaries", "Bay trawl, coastwide")
lab <- function(d) d %>% mutate(Class = factor(CLASS_LAB[class], CLASS_LAB), Gear = factor(GEAR_LAB[GEAR], GEAR_LAB),
                                Scale = factor(SCALE_LAB[scale], SCALE_LAB),
                                Group = factor(paste0(GEAR_LAB[GEAR], ", ", tolower(SCALE_LAB[scale])), GROUPS))

# -----------------------------------------------------------------------------
# 1) Values
# -----------------------------------------------------------------------------
acf_g <- ac %>% group_by(class, scale, GEAR, lag) %>%
  summarise(series = n(), med = median(acf), lo = quantile(acf, 0.25), hi = quantile(acf, 0.75), .groups = "drop") %>% lab()
slow <- dg %>% select(all_of(keys), acf1, slow_share, slow_null_q95, slow_p) %>% lab()
yr <- ym %>% inner_join(dg %>% select(all_of(keys), sd_departure), by = keys) %>% mutate(z = mean_dep / sd_departure) %>%
  group_by(class, scale, GEAR, YEAR) %>% summarise(series = n(), med = median(z), .groups = "drop") %>% lab()
# Coastwide diversity, each order separately (a, e) and its label (c).
qlab <- function(t) paste0("q = ", sub("^q", "", t))
acf_cwd <- ac %>% filter(class == "diversity", scale == "coastwide") %>% lab() %>% mutate(qlab = qlab(target))
yr_cwd <- ym %>% filter(class == "diversity", scale == "coastwide") %>%
  inner_join(dg %>% select(all_of(keys), sd_departure), by = keys) %>% mutate(z = mean_dep / sd_departure, qlab = qlab(target)) %>% lab()
if (nrow(distinct(acf_cwd, GEAR, target)) != 6L || nrow(distinct(yr_cwd, GEAR, target)) != 6L) stop("Expected 2 gears x 3 orders coastwide", call. = FALSE)
if (anyNA(yr$med) || anyNA(acf_g$med)) stop("Missing values in the group summaries", call. = FALSE)
readr::write_csv(acf_g, out[1]); readr::write_csv(slow, out[2]); readr::write_csv(yr, out[3])
above <- slow %>% group_by(class, scale, GEAR) %>% summarise(series = n(), above_ar1 = sum(slow_share > slow_null_q95), .groups = "drop")
print(as.data.frame(above), row.names = FALSE)
# Numbers quoted in the Fig. S3 caption (g20): series above the AR(1) line.
CAPTION_ABOVE <- c(diversity_bay = 42L, diversity_coastwide = 6L, abundance_bay = 95L, abundance_coastwide = 13L)
got <- above %>% group_by(class, scale) %>% summarise(a = sum(above_ar1), .groups = "drop")
got <- setNames(got$a, paste(got$class, got$scale, sep = "_"))[names(CAPTION_ABOVE)]
if (!identical(unname(got), unname(CAPTION_ABOVE))) stop("Series above the AR(1) line differ from the Figure S3 caption in g20: ",
  paste(names(got), got, collapse = ", "), call. = FALSE)

# -----------------------------------------------------------------------------
# 2) Panels
# -----------------------------------------------------------------------------
# One legend for gear and scale: color and symbol shape for gear, line type and
# symbol fill (filled bays, open coastwide) for scale.
G_COL <- setNames(GEAR_COLS[c("Bag seine", "Bag seine", "Bay trawl", "Bay trawl")], GROUPS)
G_LT  <- setNames(c("solid", "22", "solid", "22"), GROUPS)
G_SH  <- setNames(c(16, 1, 17, 2), GROUPS)
grp_scales <- list(scale_colour_manual(values = G_COL, name = NULL, drop = FALSE),
                   scale_linetype_manual(values = G_LT, name = NULL, drop = FALSE),
                   scale_shape_manual(values = G_SH, name = NULL, drop = FALSE))
ACF_MARK <- c(1, 6, 12, 18, 24, 30, 36)
# Right-end labels for the coastwide diversity orders (one per line).
end_labels <- function(d, x, y) {
  last <- d %>% group_by(GEAR, target) %>% filter(.data[[x]] == max(.data[[x]])) %>% ungroup()
  ggrepel::geom_text_repel(data = last, aes(x = .data[[x]], y = .data[[y]], label = qlab, colour = Group), inherit.aes = FALSE,
                           size = (FONT_SIZE - 2) / .pt, show.legend = FALSE, direction = "y", hjust = 0, nudge_x = 0.04 * diff(range(d[[x]])),
                           min.segment.length = 0.1, segment.size = 0.2, box.padding = 0.1, seed = 1)
}
pan_acf <- function(cls, tag) {
  d <- acf_g %>% filter(class == cls, !(cls == "diversity" & scale == "coastwide"))
  cw <- if (cls == "diversity") acf_cwd else NULL
  p <- ggplot(d, aes(lag, med, colour = Group, linetype = Group, shape = Group)) +
    geom_hline(yintercept = 0, colour = PAL[["grey_light"]], linewidth = 0.4) +
    geom_ribbon(data = d %>% filter(scale == "bay"), mapping = aes(x = lag, ymin = lo, ymax = hi, fill = Gear),
                inherit.aes = FALSE, colour = NA, alpha = 0.15, show.legend = FALSE) +
    geom_line(linewidth = 0.5) + geom_point(data = d %>% filter(lag %in% ACF_MARK), size = 1.5, stroke = 0.6) +
    grp_scales + scale_fill_manual(values = GEAR_COLS, guide = "none") +
    scale_x_continuous(breaks = ACF_MARK, expand = expansion(mult = c(0.02, if (is.null(cw)) 0.02 else 0.13))) +
    labs(title = CLASS_LAB[[cls]], tag = tag, x = "Months apart", y = "Correlation of departures") + theme_paper()
  if (!is.null(cw)) p <- p +
    geom_line(data = cw, aes(lag, acf, group = interaction(GEAR, target)), linewidth = 0.4, show.legend = FALSE) +
    geom_point(data = cw %>% filter(lag %in% ACF_MARK), aes(lag, acf), size = 1.3, stroke = 0.5, show.legend = FALSE) +
    end_labels(cw, "lag", "acf")
  p
}
pan_slow <- function(cls, tag) {
  d <- slow %>% filter(class == cls); lim <- range(c(d$slow_share, d$slow_null_q95))
  ggplot(d, aes(slow_null_q95, slow_share, colour = Group, linetype = Group, shape = Group)) +
    geom_abline(slope = 1, intercept = 0, colour = PAL[["grey_mid"]], linewidth = 0.4) +
    geom_point(size = 1.5, stroke = 0.6) +
    (if (cls == "diversity") ggrepel::geom_text_repel(
       data = d %>% filter(scale == "coastwide") %>% mutate(qlab = qlab(target)), aes(label = qlab), size = (FONT_SIZE - 2) / .pt,
       show.legend = FALSE, min.segment.length = 0.2, segment.size = 0.2, box.padding = 0.2, seed = 1, nudge_x = 0.06,
       direction = "y", hjust = 0) else NULL) +
    scale_x_continuous(limits = lim) + scale_y_continuous(limits = lim) +
    grp_scales +
    labs(tag = tag, x = "Proportion of variance in the 25-month\nrunning mean, AR(1) series (95th percentile)",
         y = "Proportion of variance in the 25-month\nrunning mean, observed") + theme_paper()
}
pan_year <- function(cls, tag) {
  d <- yr %>% filter(class == cls, !(cls == "diversity" & scale == "coastwide"))
  cw <- if (cls == "diversity") yr_cwd else NULL
  p <- ggplot(d, aes(YEAR, med, colour = Group, linetype = Group, shape = Group)) +
    geom_hline(yintercept = 0, colour = PAL[["grey_light"]], linewidth = 0.4) + geom_line(linewidth = 0.5) +
    geom_point(data = d %>% filter(YEAR %% 5 == 0), size = 1.5, stroke = 0.6) + grp_scales +
    scale_x_continuous(breaks = seq(1980, 2020, 10), expand = expansion(mult = c(0.02, if (is.null(cw)) 0.02 else 0.13))) +
    labs(tag = tag, x = "Year", y = "Standardized yearly\nmean departure") + theme_paper()
  if (!is.null(cw)) p <- p +
    geom_line(data = cw, aes(YEAR, z, group = interaction(GEAR, target)), linewidth = 0.4, show.legend = FALSE) +
    geom_point(data = cw %>% filter(YEAR %% 5 == 0), aes(YEAR, z), size = 1.3, stroke = 0.5, show.legend = FALSE) +
    end_labels(cw, "YEAR", "z")
  p
}
# One legend, taken from panel b only: with ggplot2 on Windows the panels' legends
# did not merge, so the others are switched off.
no_leg <- guides(colour = "none", linetype = "none", shape = "none")
one_leg <- guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
fig <- (pan_acf("diversity", "a") + no_leg | pan_acf("abundance", "b") + one_leg) /
  (pan_slow("diversity", "c") + no_leg | pan_slow("abundance", "d") + no_leg) /
  (pan_year("diversity", "e") + no_leg | pan_year("abundance", "f") + no_leg) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", plot.tag = element_text(size = FONT_SIZE + 1, face = "bold"))
save_figure(fig, paste0(FIGURE_LABEL, "_slow_patterns"), W_MM, H_MM)
for (f in out) cat("Saved:", f, "\n")
end_log("finished")
cat("g12_slow_patterns.R completed\n")
