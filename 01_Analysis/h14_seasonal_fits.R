# =============================================================================
# Script:  h14_seasonal_fits.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds Table S4 (Online Resource 2): the whole-record seasonal curves
#   (harmonic regression with periods of 12, 6 and 4 months) behind Fig. S1.
#   (a) Coastwide biodiversity; (b) biodiversity in the estuaries, as median
#   [minimum, maximum] over the seven estuaries, with the peak months; (c)
#   abundance of the eight taxa, estuaries as median [minimum, maximum] and
#   coastwide. The note says how far the curve refitted at each forecast date
#   moves.
#
# Run after: c01, c02, e03
#
# Inputs:
#   01_Analysis/Data/c1_seasonality_characterization/seasonality_summary.csv,
#     seasonality_fits.rds (c01)
#   01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/gamma_q012_deseas_scaled.rds (c02)
#   02_Results/03_Diagnostics/c2_deseasonalize_alpha_gamma_scaled/alpha_bay_harmonic_fits_scaled.csv,
#     alpha_residual_seasonality_scaled.csv, alpha_residual_acf_scaled.csv,
#     seasonal_by_origin_checks_scaled.csv (c02)
#   02_Results/03_Diagnostics/e3_deseasonalize_abundance/e3_harmonic_fits.csv,
#     e3_residual_seasonality.csv, e3_seasonal_by_origin_checks.csv (e03)
#
# Outputs:
#   01_Analysis/tables/Table_S4_seasonal_fits.docx
#   01_Analysis/Data/h14_seasonal_fits/seasonal_fits_a_coastwide_diversity.csv,
#     seasonal_fits_b_bay_diversity.csv, seasonal_fits_c_abundance.csv
#   01_Analysis/Logs/h14_seasonal_fits/run_<timestamp>.txt
#   Earlier copies of these CSV and Word files are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   The values describe the whole-record curves, shown for description only.
#   The forecasts use a curve refitted at each forecast date on the months up
#   to that date; the note reports how far that curve moves, from the
#   by-origin checks of c02 and e03.
#   Leftover seasonality is the share of the variance of the departures
#   explained by calendar month (observed months only). Its chance level
#   after a three-harmonic fit is (11 - 2 x 3) / (n - 1). For the coastwide
#   biodiversity series it is computed here from the c02 series; the other
#   values are read from c02 and e03.
#   No Ljung-Box column is shown; the note states that the departures of the
#   coastwide biodiversity series carry month-to-month memory (the memory the
#   forecast models use). The note also names the taxon with the highest
#   leftover seasonality. The script stops if either statement no longer
#   holds (Ljung-Box p not below 0.001; highest leftover not brown shrimp, or
#   another taxon above 4%), so the note cannot go out of date.
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
need <- c("dplyr", "readr", "officer", "flextable")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(dplyr); library(officer); library(flextable)})

TABLE_LABEL <- "Table S4"   # Online Resource 2 number
script_stem <- "h14_seasonal_fits"
C1 <- file.path(PATHS$data_dir, "c1_seasonality_characterization")
C2 <- file.path(PATHS$results_dir, "03_Diagnostics", "c2_deseasonalize_alpha_gamma_scaled")
E3 <- file.path(PATHS$results_dir, "03_Diagnostics", "e3_deseasonalize_abundance")
inputs <- c(c1_sum = file.path(C1, "seasonality_summary.csv"), c1_fits = file.path(C1, "seasonality_fits.rds"),
            c2_gamma = file.path(PATHS$data_dir, "c2_deseasonalize_alpha_gamma_scaled", "gamma_q012_deseas_scaled.rds"),
            c2_fits = file.path(C2, "alpha_bay_harmonic_fits_scaled.csv"), c2_res = file.path(C2, "alpha_residual_seasonality_scaled.csv"),
            c2_acf = file.path(C2, "alpha_residual_acf_scaled.csv"), c2_orig = file.path(C2, "seasonal_by_origin_checks_scaled.csv"),
            e3_fits = file.path(E3, "e3_harmonic_fits.csv"), e3_res = file.path(E3, "e3_residual_seasonality.csv"),
            e3_orig = file.path(E3, "e3_seasonal_by_origin_checks.csv"))
dirs <- script_dirs(script_stem, c("data", "logs"))
tab_dir <- file.path(PATHS$analysis_dir, "tables"); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
start_log(dirs$logs, script_stem)
stop_if_missing(inputs, c(rep("c01_seasonality_characterization.R", 2), rep("c02_deseasonalize_alpha_gamma_scaled.R", 5),
                          rep("e03_deseasonalize_abundance.R", 3)))
out_docx <- file.path(tab_dir, "Table_S4_seasonal_fits.docx")
out_csv <- file.path(dirs$data, paste0("seasonal_fits_", c("a_coastwide_diversity", "b_bay_diversity", "c_abundance"), ".csv"))
old <- c(out_csv, out_docx)[file.exists(c(out_csv, out_docx))]; if (length(old)) archive_files(old)
rd <- function(f) readr::read_csv(f, show_col_types = FALSE, progress = FALSE)

GEARS <- c(BS = "Bag seine", BT = "Bay trawl")
Q_LAB <- c("0" = "q = 0", "1" = "q = 1", "2" = "q = 2")   # plain labels
month_r2 <- function(x, month) {   # as in c02: share of variance explained by calendar-month means
  ok <- is.finite(x); x <- x[ok]; month <- month[ok]
  sum((ave(x, month, FUN = mean) - mean(x))^2) / sum((x - mean(x))^2)
}
mmr <- function(v, d) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), median(v), min(v), max(v))
f1 <- function(v) sprintf("%.1f", v); f2 <- function(v) sprintf("%.2f", v)

# -----------------------------------------------------------------------------
# 1) (a) Coastwide diversity
# -----------------------------------------------------------------------------
c1 <- rd(inputs[["c1_sum"]]); c1_fits <- readRDS(inputs[["c1_fits"]])
gam <- readRDS(inputs[["c2_gamma"]])
left_g <- gam %>% filter(!filled) %>% group_by(GEAR, q) %>%
  summarise(left = month_r2(value_deseas, MONTH), chance = (11 - 2 * 3) / (n() - 1), .groups = "drop")
peak_g <- bind_rows(lapply(names(c1_fits), function(k) {
  p <- strsplit(k, "__q")[[1]]; v <- c1_fits[[k]]$seasonal_curve$value_hat
  data.frame(GEAR = p[1], q = as.integer(p[2]), peak = month.abb[which.max(v)], trough = month.abb[which.min(v)])
}))
a <- c1 %>% mutate(q = as.integer(q)) %>% left_join(left_g, by = c("GEAR", "q")) %>% left_join(peak_g, by = c("GEAR", "q")) %>% arrange(GEAR, q)
if (nrow(a) != 6 || anyNA(a$left) || anyNA(a$peak)) stop("Coastwide biodiversity: expected 6 complete rows", call. = FALSE)
if (any(a$ljungbox_p_lag12 >= 0.001)) stop("A coastwide Ljung-Box p is not below 0.001; the note needs changing", call. = FALSE)
tab_a <- a %>% transmute(Gear = GEARS[GEAR], `Hill number` = Q_LAB[as.character(q)], Months = as.character(n_obs),
                         `Annual amplitude` = f2(seasonal_amplitude), `Peak to trough` = f2(seasonal_peak_to_trough),
                         `Peak, trough` = paste(peak, trough, sep = ", "), `Variance explained (%)` = f1(100 * variance_explained_season),
                         `Leftover seasonality (%)` = f1(100 * left), `Lag-1 correlation` = f2(resid_acf_lag1))

# -----------------------------------------------------------------------------
# 2) (b) Bay diversity, median [min, max] over the seven estuaries
# -----------------------------------------------------------------------------
x <- rd(inputs[["c2_fits"]]) %>%
  left_join(rd(inputs[["c2_res"]]) %>% select(BAY_ID, GEAR, q, chance_level, month_r2_new, month_r2_coastwide_curve), by = c("BAY_ID", "GEAR", "q")) %>%
  left_join(rd(inputs[["c2_acf"]]) %>% select(BAY_ID, GEAR, q, acf1_new), by = c("BAY_ID", "GEAR", "q"))
if (nrow(x) != 42 || anyNA(x$month_r2_new) || anyNA(x$acf1_new)) stop("Bay biodiversity: expected 42 complete rows", call. = FALSE)
b <- x %>% group_by(GEAR, q) %>% arrange(BAY_ID, .by_group = TRUE) %>%
  summarise(n_bays = n(), r2_med = median(r2_seasonal), r2_min = min(r2_seasonal), r2_max = max(r2_seasonal),
            rel_med = median(rel_amp_annual), rel_min = min(rel_amp_annual), rel_max = max(rel_amp_annual),
            peaks = { tb <- table(factor(peak_month, levels = 1:12)); tb <- tb[tb > 0]
                      paste0(month.abb[as.integer(names(tb))], " (", as.integer(tb), ")", collapse = ", ") },
            cw_med = median(month_r2_coastwide_curve), cw_min = min(month_r2_coastwide_curve), cw_max = max(month_r2_coastwide_curve),
            own_med = median(month_r2_new), own_min = min(month_r2_new), own_max = max(month_r2_new), chance = mean(chance_level),
            acf_med = median(acf1_new), acf_min = min(acf1_new), acf_max = max(acf1_new),
            `Variance explained (%)` = mmr(100 * r2_seasonal, 1), `Relative annual amplitude` = mmr(rel_amp_annual, 2),
            `Leftover, coastwide curve (%)` = mmr(100 * month_r2_coastwide_curve, 1), `Leftover, own curve (%)` = mmr(100 * month_r2_new, 1),
            `Lag-1 correlation` = mmr(acf1_new, 2), .groups = "drop") %>% arrange(GEAR, q)
if (any(b$n_bays != 7L)) stop("Bay biodiversity: expected seven estuaries per row", call. = FALSE)
tab_b <- b %>% transmute(Gear = GEARS[GEAR], `Hill number` = Q_LAB[as.character(q)], `Variance explained (%)`,
                         `Relative annual amplitude`, `Peak months` = peaks, `Leftover, coastwide curve (%)`,
                         `Leftover, own curve (%)`, `Lag-1 correlation`)

# -----------------------------------------------------------------------------
# 3) (c) Abundance, bays median [min, max] and coastwide
# -----------------------------------------------------------------------------
ab <- rd(inputs[["e3_fits"]]) %>% mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
  left_join(rd(inputs[["e3_res"]]) %>% mutate(SPECIES_CODE = as.character(SPECIES_CODE)) %>%
              select(scale, BAY_ID, GEAR, SPECIES_CODE, month_r2, chance_level, acf1), by = c("scale", "BAY_ID", "GEAR", "SPECIES_CODE"))
if (nrow(ab) != 128 || anyNA(ab$month_r2) || anyNA(ab$acf1)) stop("Abundance: expected 128 complete rows", call. = FALSE)
cc <- ab %>% group_by(GEAR, SPECIES_CODE, Common_Name) %>%
  summarise(n_bays = sum(scale == "bay"),
            r2_bay = mmr(100 * r2_seasonal[scale == "bay"], 1), r2_cw = f1(100 * r2_seasonal[scale == "coastwide"]),
            left_bay = mmr(100 * month_r2[scale == "bay"], 1), left_cw = f1(100 * month_r2[scale == "coastwide"]),
            acf_bay = mmr(acf1[scale == "bay"], 2), acf_cw = f2(acf1[scale == "coastwide"]),
            r2_bay_med = median(r2_seasonal[scale == "bay"]), .groups = "drop") %>%
  arrange(GEAR, desc(r2_bay_med))
if (nrow(cc) != 16 || any(cc$n_bays != 7L)) stop("Abundance: expected 16 gear-taxon rows with seven estuaries", call. = FALSE)
tab_c <- cc %>% transmute(Gear = GEARS[GEAR], Taxon = Common_Name, `Variance explained, estuaries (%)` = r2_bay,
                          `Variance explained, coastwide (%)` = r2_cw, `Leftover, estuaries (%)` = left_bay, `Leftover, coastwide (%)` = left_cw,
                          `Lag-1 correlation, estuaries` = acf_bay, `Lag-1 correlation, coastwide` = acf_cw)

# -----------------------------------------------------------------------------
# 4) How far the curve refitted at each forecast date moves (note)
# -----------------------------------------------------------------------------
orig <- bind_rows(rd(inputs[["c2_orig"]]) %>% mutate(class = "diversity"), rd(inputs[["e3_orig"]]) %>% mutate(class = "abundance"))
mv <- orig %>% group_by(class) %>% summarise(last_max = max(last_vs_full_sd), step_med = median(median_step_sd), step_max = max(max_step_sd), .groups = "drop")
print(as.data.frame(mv), row.names = FALSE)
chance <- range(c(left_g$chance, b$chance, ab$chance_level))
cat(sprintf("Chance level of leftover seasonality: %.2f%% to %.2f%%\n", 100 * chance[1], 100 * chance[2]))
cat(sprintf("Leftover, own curve: coastwide biodiversity %.1f%% to %.1f%%; bay biodiversity median %.1f%% (max %.1f%%); abundance median %.1f%% (max %.1f%%)\n",
            100 * min(a$left), 100 * max(a$left), 100 * median(x$month_r2_new), 100 * max(x$month_r2_new),
            100 * median(ab$month_r2), 100 * max(ab$month_r2)))

top_left <- ab %>% slice_max(month_r2, n = 1)
if (top_left$Common_Name != "Brown shrimp") stop("Highest leftover seasonality is not brown shrimp; the note needs changing", call. = FALSE)
bsr <- ab %>% filter(Common_Name == "Brown shrimp")
bs <- list(bt_bay = 100 * max(bsr$month_r2[bsr$GEAR == "BT" & bsr$scale == "bay"]), bt_cw = 100 * bsr$month_r2[bsr$GEAR == "BT" & bsr$scale == "coastwide"],
           bs_bay = 100 * max(bsr$month_r2[bsr$GEAR == "BS" & bsr$scale == "bay"]), bs_cw = 100 * bsr$month_r2[bsr$GEAR == "BS" & bsr$scale == "coastwide"])
if (any(ab$month_r2[ab$Common_Name != "Brown shrimp"] > 0.04)) stop("Another taxon has leftover seasonality above 4%; the note needs changing", call. = FALSE)

# -----------------------------------------------------------------------------
# 5) Word table, three panels
# -----------------------------------------------------------------------------
style_ft <- function(df, widths) {
  n <- ncol(df)
  flextable(as.data.frame(df) |> setNames(sub("^Gear$", "Sampling method", names(df)))) %>% font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 8, part = "body") %>% fontsize(size = 8.5, part = "header") %>% bold(part = "header") %>%
    align(j = 1:2, align = "left", part = "all") %>% align(j = 3:n, align = "center", part = "all") %>%
    border_remove() %>% hline_top(border = fp_border(width = 1.5), part = "header") %>%
    hline_bottom(border = fp_border(width = 0.75), part = "header") %>% hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
    merge_v(j = 1) %>% valign(j = 1, valign = "top", part = "body") %>%
    padding(padding.top = 1, padding.bottom = 1, padding.left = 3, padding.right = 3, part = "all") %>%
    width(width = widths * 9.4 / sum(widths))
}
ft_a <- style_ft(tab_a, c(0.8, 1.5, 0.6, 0.8, 0.8, 0.85, 0.95, 0.95, 0.95))
ft_b <- style_ft(tab_b, c(0.8, 1.2, 1.15, 1.15, 1.6, 1.15, 1.1, 1.15))   # wider peak-months column (counts added)
ft_c <- style_ft(tab_c, c(0.8, 1.3, 1.15, 0.9, 1.1, 0.9, 1.1, 0.95))
caption <- paste0(
  "Seasonal curves fitted to the whole record (harmonic regression with 12-, 6-, and 4-month waves), as shown ",
  "in Fig. S1: (a) coastwide biodiversity; (b) biodiversity in the seven estuaries; (c) abundance, log(CPUE + 1), of the eight taxa, ",
  "estuaries and coastwide. In (b) and in the estuary columns of (c), each cell summarizes the seven estuary values as median ",
  "[minimum, maximum]. Peak months (b): the months in which the estuaries' curves peak, with the number of estuaries in parentheses. Every series has ",
  "its own curve. Variance explained: proportion of the series' variance accounted for by the curve. Leftover seasonality: proportion ",
  "of the variance of the departures from the curve still explained by calendar month (expected by chance approximately ",
  sprintf("%.1f", 100 * mean(chance)), "%). Lag-1 correlation: the correlation between the departures of consecutive months. Annual amplitude and peak to trough ",
  "(a) are in effective species; relative annual amplitude (b) is the annual amplitude divided by the series mean.")
note <- paste0(
  "The forecasts do not use these curves: at each forecast date the curve was refitted on the months up to that date. ",
  sprintf("From one forecast date to the next the refitted curve moved by a median of %.2f standard deviations of the departures (at most %.2f); ",
          median(orig$median_step_sd), max(orig$max_step_sd)),
  sprintf("at the last forecast date it differed from the whole-record curve by at most %.2f. ", max(orig$last_vs_full_sd)),
  "In every coastwide biodiversity series the departures were autocorrelated (Ljung-Box test at lag 12, p < 0.001), ",
  "the month-to-month memory that the forecast models use. Leftover with the coastwide curve (b): the same measure ",
  "after subtracting the coastwide curve from each estuary's series instead of the estuary's own curve. April and May 2020 ",
  "(no sampling) are left out of every value. Leftover seasonality was highest for brown shrimp ",
  sprintf("(bay trawl: estuaries up to %.1f%%, coastwide %.1f%%; bag seine: up to %.1f%%, coastwide %.1f%%), ", bs$bt_bay, bs$bt_cw, bs$bs_bay, bs$bs_cw),
  "whose short spring peak the three waves do not follow fully.")
ptxt <- function(s, size = 10, bold = FALSE) fpar(ftext(s, fp_text(font.size = size, bold = bold, font.family = "Times New Roman")))
doc <- read_docx() %>%
  body_set_default_section(prop_section(page_size = page_size(orient = "landscape", width = 8.5, height = 11),
                                        page_margins = page_mar(top = 0.75, bottom = 0.75, left = 0.75, right = 0.75))) %>%
  body_add_fpar(fpar(ftext(paste0(TABLE_LABEL, " "), fp_text(font.size = 10, bold = TRUE, font.family = "Times New Roman")),
                     ftext(caption, fp_text(font.size = 10, font.family = "Times New Roman")))) %>%
  body_add_fpar(ptxt("(a) Coastwide biodiversity", 9, TRUE)) %>% body_add_flextable(ft_a) %>%
  body_add_fpar(ptxt("(b) Biodiversity in the seven estuaries", 9, TRUE)) %>% body_add_flextable(ft_b) %>%
  body_add_break() %>% body_add_fpar(ptxt("(c) Abundance", 9, TRUE)) %>% body_add_flextable(ft_c) %>%
  body_add_fpar(ptxt(note, 9))
print(doc, target = out_docx)
readr::write_csv(a %>% select(GEAR, q, n_obs, seasonal_amplitude, seasonal_peak_to_trough, peak, trough, variance_explained_season, left, chance, resid_acf_lag1, ljungbox_p_lag12), out_csv[1])
readr::write_csv(b %>% select(GEAR, q, n_bays, r2_med:acf_max), out_csv[2])
readr::write_csv(cc, out_csv[3])
for (f in c(out_docx, out_csv)) cat("Saved:", f, "\n")
end_log("finished")
cat("h14_seasonal_fits.R completed\n")
