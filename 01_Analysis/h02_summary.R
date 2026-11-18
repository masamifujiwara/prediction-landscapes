# =============================================================================
# Script:  h02_summary.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds Table S6 (Online Resource 2): forecast accuracy, skill and
#   management situations of SARIMA (automatic), one row per target class,
#   scale and sampling method. Upper table: number of series; median RMSE
#   (with interquartile range) 1, 12 and 24 months ahead; median percentage of
#   monthly forecasts within the tolerance 1 and 24 months ahead; median skill
#   against the cumulative average (1, 12, 24 months) and against the
#   12-month average (1, 12 months); numbers of series in the four management
#   situations 1 and 12 months ahead. Lower table: the division of the
#   month-to-month variation shown in Fig. 2, as means and SD (%) across series.
#
# Run after: f01
#
# Inputs:
#   01_Analysis/Data/f1_forecast_scoring/f1_by_lead.rds, f1_quadrants.rds,
#     f1_variance_budget.rds (f01)
#
# Outputs:
#   01_Analysis/tables/Table_S6_summary.csv (upper table values)
#   01_Analysis/tables/Table_S6_partition.csv (lower table values)
#   01_Analysis/tables/Table_S6_summary.docx
#   01_Analysis/Logs/h2_summary/run_<timestamp>.txt
#   An earlier Table_S6_summary.docx is moved to ARCHIVE_<timestamp> first.
#
# Method notes:
#   SARIMA (automatic) is the base model shown. RMSE is a proportion of the
#   series mean for biodiversity and in log units for abundance. The
#   management situations (forecast-led, mean-led, forecast with buffer,
#   precaution; Fig. 7) are counted by the uncertainty-range rule: a series
#   counts as within the tolerance or better than the cumulative average only
#   when its whole 95% range clears the line. The counts are given here
#   because the figure shows the points only. The lower table uses the same
#   input as g02, and the script checks that the four parts add to 1 and that
#   both tables have the same series counts.
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
source(file.path(PATHS$analysis_dir, "lib", "figures.R"))   # labels
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(officer); library(flextable)})
script_stem <- "h2_summary"
in_lead <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_by_lead.rds")
in_quad <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_quadrants.rds")
in_vb   <- file.path(PATHS$data_dir, "f1_forecast_scoring", "f1_variance_budget.rds")
MODEL <- "sarima_auto"
QUAD <- c("forecast-led", "mean-led", "forecast with buffer", "precaution")
CLASS_ROW <- c(diversity = "Biodiversity (RMSE as a proportion of the series mean)", abundance = "Abundance (RMSE in log units)")
TABLE_LABEL <- "Table S6"   # Online Resource 2 number
CAPTION <- paste(
  "Forecast accuracy, skill, and management situations by group (SARIMA forecasts; medians over the series of",
  "each group, with interquartile ranges in parentheses for accuracy). RMSE: root mean squared out-of-sample prediction error, as a",
  "proportion of the series mean for biodiversity and in log units (log(CPUE + 1)) for abundance, 1, 12, and 24 months ahead. Within",
  "tolerance: percentage of monthly forecasts within 30% of the series mean (biodiversity) or 0.5 log units (abundance), 1 and 24 months",
  "ahead. Skill: 1 minus the ratio of the mean squared prediction error to that of the cumulative average (all data up to the forecast date)",
  "or of the 12-month average (the average of the last 12 months); 0 means no better than that reference. Situations: numbers of",
  "series that are forecast-led (FL), mean-led (ML), forecast with buffer (FB), or precaution (P) (Fig. 7); a series counts as within",
  "the tolerance or better than the cumulative average only when the whole 95% uncertainty range clears the line. Series counts:",
  "biodiversity 21 in the estuaries (7 estuaries \u00d7 q = 0, 1, 2) and 3 coastwide per sampling method; abundance 56 in the",
  "estuaries (7 estuaries \u00d7 8 taxa) and 8 coastwide. Calculated on the 367 forecast dates up to December 2022.",
  "Lower table: the division of the month-to-month variation into the seasonal cycle, slow change, short-lived carry-over,",
  "and the remainder (Fig. 2), with the sum of the first three (predictable), as means \u00b1 SD (%) across the series of each group.")

dirs <- script_dirs(script_stem, c("logs"))
start_log(dirs$logs, script_stem)
stop_if_missing(c(in_lead, in_quad, in_vb), rep("f01_forecast_scoring.R", 3))
tab_dir <- file.path(PATHS$analysis_dir, "tables"); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(tab_dir, "Table_S6_summary.csv"); out_docx <- file.path(tab_dir, "Table_S6_summary.docx")
out_vb_csv <- file.path(tab_dir, "Table_S6_partition.csv")

bl <- readRDS(in_lead) %>% filter(model == MODEL)
qd <- readRDS(in_quad) %>% filter(model == MODEL)
if (!all(QUAD %in% unique(qd$quadrant))) stop("Unexpected quadrant labels: ", paste(unique(qd$quadrant), collapse = ", "), call. = FALSE)
med <- function(col, h) bl %>% filter(lead == h) %>% group_by(class, scale, GEAR) %>%
  summarise(series = n(), !!paste0(col, "_", h) := stats::median(.data[[col]]),
            !!paste0(col, "_", h, "_q1") := unname(stats::quantile(.data[[col]], 0.25)),
            !!paste0(col, "_", h, "_q3") := unname(stats::quantile(.data[[col]], 0.75)), .groups = "drop")
tab <- Reduce(function(a, b) left_join(a, b %>% select(-series), by = c("class", "scale", "GEAR")),
              list(med("rmse_scaled", 1), med("rmse_scaled", 12), med("rmse_scaled", 24), med("hit_rate", 1), med("hit_rate", 24),
                   med("skill_past", 1), med("skill_past", 12), med("skill_past", 24), med("skill_m12", 1), med("skill_m12", 12)))
qc <- qd %>% filter(lead %in% c(1L, 12L)) %>% mutate(quadrant = factor(quadrant, QUAD)) %>%
  count(class, scale, GEAR, lead, quadrant, .drop = FALSE) %>%
  group_by(class, scale, GEAR, lead) %>% summarise(sit = paste(n, collapse = " / "), .groups = "drop") %>%
  pivot_wider(names_from = lead, values_from = sit, names_prefix = "situations_")
tab <- tab %>% left_join(qc, by = c("class", "scale", "GEAR")) %>%
  mutate(class = factor(class, c("diversity", "abundance")), scale = factor(scale, c("bay", "coastwide"))) %>% arrange(class, scale, GEAR)
if (nrow(tab) != 8L || sum(tab$series) != 176L) stop("Table S6 must have 8 rows and 176 series", call. = FALSE)
readr::write_csv(tab %>% mutate(across(where(is.double), ~ round(.x, 3))), out_csv)

# ---- Lower table: division of the variation (Fig. 2), means and SD across series
VB_PARTS <- c(seasonal = "Seasonal cycle", slow = "Slow change", short_lived = "Short-lived carry-over",
              unpredictable = "Remainder", predictable = "Predictable (sum)")
vb <- readRDS(in_vb) %>% filter(model == MODEL) %>% mutate(predictable = seasonal + slow + short_lived)
if (max(abs(vb$predictable + vb$unpredictable - 1)) > 1e-8) stop("The four parts do not add to 1", call. = FALSE)
vb_tab <- vb %>% group_by(class, scale, GEAR) %>%
  summarise(series = n(), across(all_of(names(VB_PARTS)), list(mean = mean, sd = stats::sd)), .groups = "drop") %>%
  mutate(class = factor(class, c("diversity", "abundance")), scale = factor(scale, c("bay", "coastwide"))) %>% arrange(class, scale, GEAR)
if (!identical(vb_tab$series, tab$series)) stop("Series counts differ between the two tables", call. = FALSE)
readr::write_csv(vb_tab %>% mutate(across(where(is.double), ~ round(100 * .x, 2))), out_vb_csv)

f2 <- function(x) sub("^-0[.]00$", "0.00", sprintf("%.2f", x));   # no "-0.00"
 f0 <- function(x) sprintf("%.0f", 100 * x)
iq2 <- function(m, lo, hi) paste0(f2(m), " (", f2(lo), " to ", f2(hi), ")")
iq0 <- function(m, lo, hi) paste0(f0(m), " (", f0(lo), " to ", f0(hi), ")")
disp <- tab %>% transmute(
  Group = paste0(SCALE_LAB[as.character(scale)], ", ", tolower(GEAR_LAB[GEAR])), Series = as.character(series),
  `RMSE, 1 month` = iq2(rmse_scaled_1, rmse_scaled_1_q1, rmse_scaled_1_q3),
  `RMSE, 12 months` = iq2(rmse_scaled_12, rmse_scaled_12_q1, rmse_scaled_12_q3),
  `RMSE, 24 months` = iq2(rmse_scaled_24, rmse_scaled_24_q1, rmse_scaled_24_q3),
  `Within tolerance (%), 1 month` = iq0(hit_rate_1, hit_rate_1_q1, hit_rate_1_q3),
  `Within tolerance (%), 24 months` = iq0(hit_rate_24, hit_rate_24_q1, hit_rate_24_q3),
  `Skill vs cumulative average 1 / 12 / 24` = paste(f2(skill_past_1), f2(skill_past_12), f2(skill_past_24), sep = " / "),
  `Skill vs 12-month average 1 / 12` = paste(f2(skill_m12_1), f2(skill_m12_12), sep = " / "),
  `Situations at 1 month (FL / ML / FB / P)` = situations_1, `Situations at 12 months (FL / ML / FB / P)` = situations_12,
  cls = as.character(class))
body <- bind_rows(lapply(c("diversity", "abundance"), function(k) {
  hdr <- as.data.frame(as.list(setNames(c(CLASS_ROW[[k]], rep("", ncol(disp) - 2)), names(disp)[-ncol(disp)])), check.names = FALSE)
  bind_rows(hdr, disp %>% filter(cls == k) %>% select(-cls))
}))
hdr_rows <- which(body$Group %in% CLASS_ROW)
ft <- flextable(body) %>% font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 8.5, part = "all") %>%
  bold(part = "header") %>% bold(i = hdr_rows, part = "body") %>% align(j = 2:ncol(body), align = "center", part = "all") %>%
  border_remove() %>% hline_top(border = fp_border(width = 1.5), part = "header") %>%
  hline_bottom(border = fp_border(width = 0.75), part = "header") %>% hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  padding(padding.top = 1, padding.bottom = 1, padding.left = 2, padding.right = 2, part = "all")
for (i in hdr_rows) ft <- merge_at(ft, i = i, j = seq_len(ncol(body)), part = "body")
col_w <- c(1.25, 0.45, 1.0, 1.0, 1.0, 0.85, 0.85, 1.0, 0.75, 0.85, 0.85); ft <- width(ft, width = col_w * 9.5 / sum(col_w))

pm <- function(m, s) paste0(sprintf("%.0f", 100 * m), " \u00b1 ", sprintf("%.0f", 100 * s))
VB_CLASS_ROW <- c(diversity = "Biodiversity", abundance = "Abundance")
disp2 <- vb_tab %>% transmute(Group = paste0(SCALE_LAB[as.character(scale)], ", ", tolower(GEAR_LAB[GEAR])), Series = as.character(series),
  `Seasonal cycle (%)` = pm(seasonal_mean, seasonal_sd), `Slow change (%)` = pm(slow_mean, slow_sd),
  `Short-lived carry-over (%)` = pm(short_lived_mean, short_lived_sd), `Remainder (%)` = pm(unpredictable_mean, unpredictable_sd),
  `Predictable, sum (%)` = pm(predictable_mean, predictable_sd), cls = as.character(class))
body2 <- bind_rows(lapply(c("diversity", "abundance"), function(k) {
  hdr <- as.data.frame(as.list(setNames(c(VB_CLASS_ROW[[k]], rep("", ncol(disp2) - 2)), names(disp2)[-ncol(disp2)])), check.names = FALSE)
  bind_rows(hdr, disp2 %>% filter(cls == k) %>% select(-cls))
}))
hdr_rows2 <- which(body2$Group %in% VB_CLASS_ROW)
ft2 <- flextable(body2) %>% font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 8.5, part = "all") %>%
  bold(part = "header") %>% bold(i = hdr_rows2, part = "body") %>% align(j = 2:ncol(body2), align = "center", part = "all") %>%
  border_remove() %>% hline_top(border = fp_border(width = 1.5), part = "header") %>%
  hline_bottom(border = fp_border(width = 0.75), part = "header") %>% hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  padding(padding.top = 1, padding.bottom = 1, padding.left = 2, padding.right = 2, part = "all")
for (i in hdr_rows2) ft2 <- merge_at(ft2, i = i, j = seq_len(ncol(body2)), part = "body")
ft2 <- width(ft2, width = c(1.25, 0.45, rep(1.1, 5)))
if (file.exists(out_docx)) archive_files(c(out_docx))
doc <- read_docx() %>% body_set_default_section(prop_section(page_size = page_size(width = 8.5, height = 11, orient = "landscape"),
    page_margins = page_mar(top = 0.75, bottom = 0.75, left = 0.75, right = 0.75, gutter = 0))) %>%
  body_add_fpar(fpar(ftext(paste0(TABLE_LABEL, " "), fp_text(font.family = "Times New Roman", font.size = 10, bold = TRUE)), ftext(CAPTION, fp_text(font.family = "Times New Roman", font.size = 10)))) %>%
  body_add_flextable(ft) %>%
  body_add_par("") %>%
  body_add_flextable(ft2)
print(doc, target = out_docx)
cat("Saved:", out_csv, "\nSaved:", out_vb_csv, "\nSaved:", out_docx, "\n")
print(as.data.frame(body), row.names = FALSE)
print(as.data.frame(body2), row.names = FALSE)
end_log("finished")
cat("h02_summary.R completed\n")
