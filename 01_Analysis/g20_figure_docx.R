# =============================================================================
# Script:  g20_figure_docx.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Writes the Word files for the figures: one file per main figure
#   (Figure_1.docx to Figure_7.docx, image and caption) and Online Resource 3
#   (Online_Resource_Figures.docx) with Figs. S1 to S8, one figure per page,
#   each with its caption. The first page of Online Resource 3 gives its title,
#   then the article title, journal, author, affiliation and e-mail (journal
#   requirement). Every figure caption is kept in this script, in one place.
#   Run last among the g scripts.
#
# Run after: g01 to g08, g11 to g17
#
# Inputs:
#   01_Analysis/figures/Figure_<n>_<name>.tiff (g01 to g08, g11 to g17)
#
# Outputs:
#   01_Analysis/figures/docx/Figure_1.docx to Figure_7.docx
#   01_Analysis/figures/docx/Online_Resource_Figures.docx (Online Resource 3)
#   01_Analysis/Logs/g20_figure_docx/run_<timestamp>.txt
#   Earlier copies of the Word files are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   - Captions follow the journal style: "Fig. N" in bold, no punctuation after
#     the number and none at the end of the caption.
#   - The figure list (FIGS) is the only place that ties a figure number to an
#     image file and a caption; the producing g scripts write the images with
#     the same numbers.
#   - Images are placed at the full text width (6.5 in) unless that would make
#     them taller than 8 in; every canvas is at or below 8.0 / 6.5 (the check
#     in lib/figures.R), so all fill the width. Size comes from each TIFF's own
#     pixel dimensions; the script stops if a TIFF is missing or its aspect
#     differs from the canvas listed here.
#   - In Online Resource 3, figures are at most 7 in tall (OR_MAX_H) so that
#     the caption stays on the same page. A per-figure max_h sets a lower limit
#     where the caption would otherwise run onto a second page (Fig. S1, whose
#     page also carries the article header, and Fig. S3, with a long caption).
#   - TIFF images are inserted as they are (600 dpi); no conversion, so no
#     extra package is needed.
#   - Captions are in plain language; RMSE is always described as based on
#     out-of-sample prediction errors. Captions may mark italics with
#     asterisks (cap_chunks).
#   - Online Resources: 1 supplementary methods and results, 2 tables (h20),
#     3 figures (this script).
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
need <- c("officer")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(library(officer))

OR_FIGURES_TITLE <- "Online Resource 3: Figures"   # Online Resources: 1 Supplementary methods and results, 2 Tables, 3 Figures
ARTICLE_HEADER <- c(   # first page of the Online Resource file (journal requirement: title, journal, author, affiliation, e-mail)
  "Forecasting estuarine biodiversity and population abundance for conservation: spatial scale, forecast horizon, and choice of measure",
  "Biodiversity and Conservation",
  "Masami Fujiwara",
  "Department of Ecology and Conservation Biology, Texas A&M University, College Station, TX 77843-2258, USA",
  "Corresponding author: Masami Fujiwara, fujiwara@tamu.edu")
script_stem <- "g20_figure_docx"
fig_dir <- file.path(PATHS$analysis_dir, "figures"); out_dir <- file.path(fig_dir, "docx")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dirs <- script_dirs(script_stem, c("logs"))
start_log(dirs$logs, script_stem)
TEXT_W <- 6.5; TEXT_H <- 8.0
OR_MAX_H <- 7.0   # Online Resources pages carry a caption below a possibly tall figure (and the title on page 1)

# -----------------------------------------------------------------------------
# 1) Figures: number, image file, canvas (mm) and caption
# -----------------------------------------------------------------------------
cap <- function(...) paste0(...)
FIGS <- list(
  list(n = "1", file = "Figure_1_analysis_overview", w = NA, h = NA, cap = cap(
    "Overview of the analysis of the Texas Parks and Wildlife Department (TPWD) bag seine and bay trawl monitoring data. At each ",
    "monthly forecast date, the seasonal curve was refitted to the months up to that date, the forecast models and reference forecasts ",
    "were applied to the departures from it, and the forecasts were compared with the values observed later. The characteristics of ",
    "the time series data, calculated from the whole record without forecasts, help explain why series differ. The models are listed ",
    "in Table 1")),
  list(n = "2", file = "Figure_2_variance_budget", w = 170, h = 95, cap = cap(
    "Division of the month-to-month variation of the monthly series into four parts, by target class, scale, and sampling method: the ",
    "seasonal cycle, slow change (averaged over 12 to 24 months ahead and set to zero when negative), short-lived carry-over, and the ",
    "remainder. Percentages are means over the series of each group; the four parts add to 100%")),
  list(n = "3", file = "Figure_3_misses_by_lead", w = 170, h = 125, cap = cap(
    "Size of the monthly forecast misses and percentage of forecasts within the tolerance, by months ahead. Top: absolute out-of-sample ",
    "prediction errors of SARIMA forecasts, pooled over the series of each group and the forecast dates (per sampling method, 21 ",
    "biodiversity and 56 abundance series in the estuaries; 3 and 8 coastwide); solid line: median (squares: bag seine; triangles: bay ",
    "trawl); dashed lines: middle 50%. Biodiversity misses are a percentage of the series mean; abundance misses are in log units. Open ",
    "circles: median miss of the cumulative average. Gray line: tolerance. Bottom: percentage of monthly forecasts within the ",
    "tolerance. The cumulative average and the tolerance are defined in the Methods (Forecast models and references; Forecast accuracy ",
    "and skill)")),
  list(n = "4", file = "Figure_4_two_references", w = 170, h = 115, cap = cap(
    "Skill of SARIMA forecasts against the cumulative average and the 12-month average (Methods: Forecast models and references; ",
    "Forecast accuracy and skill), by months ahead. Solid lines and circles: skill against the cumulative average; dashed lines and ",
    "triangles: skill against the 12-month average. Median over the series of each group (numbers as in Fig. 3); estuaries also with ",
    "the band of the middle 50%. Where skill against the cumulative average is above skill against the 12-month average, the 12-month ",
    "average is the better reference and recent levels persist; where it is below, departures from the seasonal curve fade and the ",
    "cumulative average is the better reference")),
  list(n = "5", file = "Figure_5_scale_skill", w = 170, h = 120, cap = cap(
    "Single estuaries against the coast as a whole. Skill of SARIMA forecasts against the cumulative average, by months ahead. ",
    "Estuaries (solid lines and circles): median over the estuary series (numbers as in Fig. 3) and band of the middle 50%. Coastwide ",
    "(dashed lines and triangles): each biodiversity order (labeled); for abundance, the median over the eight taxa (thick line) and ",
    "each taxon (thin lines)")),
  list(n = "6", file = "Figure_6_shared_surprises", w = 170, h = 175, cap = cap(
    "Months in which many estuaries missed together, with water temperature and salinity. (a, b) Proportion of estuary series with a ",
    "large forecast miss (Methods: Forecast accuracy and skill) in each month, pooled over both sampling methods (42 biodiversity and ",
    "112 abundance series), drawn above the line when the observed value was higher than forecast and below when it was lower, both as ",
    "positive proportions. Dotted lines: the average proportion over all months, in each direction. (c) Water temperature and (d) ",
    "salinity: departure from the usual value for the calendar month, averaged over the estuaries and both sampling methods. Gray band: ",
    "April and May 2020, not sampled; no one-month forecasts exist for June 2020, so the rise in bag seine biodiversity between March ",
    "and June 2020 first appears in July")),
  list(n = "7", file = "Figure_7_quadrants", w = 170, h = 95, cap = cap(
    "Accuracy against skill of SARIMA forecasts for every series (48 biodiversity and 128 abundance series), 1 and 12 months ahead. ",
    "Accuracy: RMSE divided by the tolerance (below 1: within the tolerance). Skill: against the cumulative average (above 0: better ",
    "than the cumulative average). The lines and corner labels mark the four management situations (Introduction; placement of the ",
    "series: Methods: Forecast accuracy and skill): forecast-led, forecast with buffer, mean-led, and precaution. Biodiversity: blue ",
    "circles; abundance: orange triangles; estuaries: filled; coastwide: open. Points are drawn at their central values; the counts in ",
    "Table S6 place a series on the favorable side of a line only when its whole 95% uncertainty range clears the line")),
  list(n = "S1", file = "Figure_S1_seasonal_fits", w = 170, h = 200, max_h = 6.0, cap = cap(   # smaller image: page 1 also carries the article header
    "Seasonal curves and the departures from them for eight example series. Each panel: top, observed monthly values (gray) and the ",
    "seasonal curve (a harmonic regression with periods of 12, 6, and 4 months) fitted to the whole record (color); bottom, the departures, which the models forecast. ",
    "Left: bag seine; right: bay trawl; all coastwide (the fits for the single estuaries are summarized in Table S4). (a to c, e to g) ",
    "Biodiversity, q = 0, 1, 2; (d, h) brown shrimp, log(CPUE + 1), an abundant target taxon caught by both sampling methods with a ",
    "strong seasonal cycle, as an example of abundance. The forecasts used a curve refitted at each forecast date on earlier months only; ",
    "the whole-record curve is shown for description. April and May 2020 (no sampling) appear as a gap")),
  list(n = "S2", file = "Figure_S2_horizon_tolerance", w = 170, h = 120, cap = cap(
    "Percentage of series whose SARIMA forecasts stay within the tolerance, for tolerances from 20% to 60% of the series mean ",
    "(biodiversity) and 0.3 to 0.9 log units (abundance): (top) at least one month ahead; (bottom) at every month up to 24 months ahead. ",
    "Within the tolerance: the whole 95% uncertainty range of the RMSE of the out-of-sample prediction errors lies inside the tolerance ",
    "(moving-block bootstrap over forecast dates, blocks of 12 months, 2,000 resamples). Dashed vertical line: the tolerance used in the ",
    "other figures and tables (30%; 0.5 log units). Bag seine: blue with circles; bay trawl: orange with triangles; estuaries: solid ",
    "lines and filled symbols; coastwide: dashed lines and open symbols (estuaries: 21 biodiversity and 56 abundance series per sampling ",
    "method; coastwide: 3 and 8, so the coastwide biodiversity curves move in steps of 33%). Calculated on the 367 forecast dates up to ",
    "December 2022. Because forecast error grew little with lead time (Fig. 3), the bottom curves follow the top ones")),
  list(n = "S3", file = "Figure_S3_slow_patterns", w = 170, h = 190, max_h = 6.0, cap = cap(   # smaller image: the long caption ran onto a second page
    "The departures from the seasonal curve over time, for biodiversity (left: Hill numbers q = 0, 1, 2) and abundance (right: ",
    "log(CPUE + 1) of eight taxa, CPUE being catch per unit effort). Departures are each month's value minus the seasonal curve fitted ",
    "to the whole record; April and May 2020 (not sampled) are left out. Bag seine: blue with circles; bay trawl: orange with triangles; ",
    "estuaries: solid lines and filled symbols; coastwide: dashed lines and open symbols (estuaries: 21 biodiversity and 56 abundance ",
    "series per sampling method; coastwide: 3 and 8). (a, b) Correlation of each series' departures with its own departures 1 to 36 ",
    "months later; lines: median over the series of each group, except coastwide biodiversity, whose three orders are drawn separately ",
    "and labeled; shaded: middle 50% of the estuary series. (c, d) One point per series: the proportion of the variance of the ",
    "departures that lies in their 25-month running mean (slow change), against the 95th percentile of the same proportion in 1,000 ",
    "simulated first-order autoregressive series with the series' own lag-1 correlation. Points above the line carry more slow change ",
    "than their lag-1 correlation alone would give: all 48 biodiversity series and 108 of the 128 abundance series (95 of 112 in the ",
    "estuaries, 13 of 16 coastwide). (e, f) Standardized yearly mean departure: for each series, the mean of its monthly departures in ",
    "each year (at least six sampled months), divided by the standard deviation of its monthly departures so that series in different ",
    "units can be combined; lines: median over the series of each group (coastwide biodiversity: each order)")),
  list(n = "S4", file = "Figure_S4_entropy", w = 170, h = 150, cap = cap(
    "Weighted permutation entropy, expressed as 1 minus the entropy (intrinsic predictability), against the skill of SARIMA forecasts ",
    "against the cumulative average, one point per series. The entropy is calculated from the order of 4 values of the departures from the ",
    "seasonal curve: (top) 4 consecutive months, set against skill 1 month ahead; (bottom) the same calendar month in 4 consecutive ",
    "years, set against skill 12 months ahead. Intrinsic predictability runs from 0 (no order) to 1 (fully ordered); dashed line: the ",
    "median for shuffled copies of the series. Bag seine: blue circles; bay trawl: orange triangles; estuaries: filled; coastwide: open ",
    "(estuaries: 21 biodiversity and 56 abundance series per sampling method; coastwide: 3 and 8). Rank correlations (Spearman) in each ",
    "panel, with the number of series. Skill calculated on the 367 forecast dates up to December 2022; the axes differ between panels. ",
    "Nine in ten series lie below 0.06 (shuffled copies about 0.01), so the departures carry little ordinal structure. Without the 24 ",
    "estuary abundance series with many runs of zero catches, the estuary correlations for abundance are 0.54 (one-month step) and 0.31 ",
    "(twelve-month step). Permutation entropy describes how ordered a series is; it does not separate nonlinear dynamics from random ",
    "variation")),
  list(n = "S5", file = "Figure_S5_horizons", w = 170, h = 175, cap = cap(
    "Forecast horizons of SARIMA forecasts: the number of consecutive months ahead, from one month ahead, for which the forecasts (top) ",
    "stay within the tolerance, (middle) are better than the cumulative average, and (bottom) are better than the 12-month average. ",
    "Tolerance: RMSE of the out-of-sample prediction errors within 30% of the series mean (biodiversity) or 0.5 log units ",
    "(abundance). Both references are the seasonal curve plus the mean departure of all earlier months (cumulative average) or of the last 12 ",
    "months (12-month average). Bars: percentage of the series of each group in each horizon class (0 = fails one month ahead; 24 = passes at every month); ",
    "number of series in brackets. Range: the whole 95% uncertainty range must meet the condition (moving-block bootstrap over forecast ",
    "dates, blocks of 12 months, 2,000 resamples); Value: the central estimate alone. A horizon ends at the first month that fails; skill ",
    "that is lost and later returns is not counted (9 estuary biodiversity and 11 estuary abundance series; 1 and 2 coastwide). ",
    "Calculated on the 367 forecast dates up to December 2022")),
  list(n = "S6", file = "Figure_S6_interval_coverage", w = 170, h = 130, cap = cap(
    "Percentage of observed values inside the 80% (top) and 95% (bottom) prediction intervals of SARIMA forecasts, by months ahead. ",
    "Lines: median over the series of each group (estuaries: 21 biodiversity and 56 abundance series per sampling method; coastwide ",
    "abundance: 8), except coastwide biodiversity, whose three orders are drawn separately and labeled; shaded: middle 50% of the estuary ",
    "series; dashed gray line: the nominal level. Bag seine: blue with circles; bay trawl: orange with triangles; estuaries: solid lines ",
    "and filled symbols; coastwide: dashed lines and open symbols. Calculated on the 367 forecast dates up to December 2022; the vertical ",
    "scales differ between panels. The 80% intervals contained 70% to 75% of the observed values for bag seine biodiversity in the ",
    "estuaries, and the coverage of coastwide bay trawl biodiversity declined from 83% one month ahead to 70% to 72% at 16 to 24 months; ",
    "the other groups had 78% to 84%")),
  list(n = "S7", file = "Figure_S7_diversity_periods", w = 170, h = 105, cap = cap(
    "Skill of SARIMA forecasts for the estuary series by period of the target month, one and 12 months ahead, against the cumulative average ",
    "(filled) and against the 12-month average (open). Biodiversity: blue circles; abundance: orange triangles. Median over ",
    "the estuary series of each class and sampling method (21 biodiversity and 56 abundance series per sampling method), with the range ",
    "of the middle 50%. Target months per period: 188, 144, and 57 one month ahead (177, 144, and 56 twelve months ahead); April and May ",
    "2020 were not sampled. Where skill is high against the cumulative average but near zero against the 12-month average, as for bag seine ",
    "biodiversity in 2020 to 2024, the forecasts follow a level that had already shifted rather than anticipate the shift")),
  list(n = "S8", file = "Figure_S8_taxon_skill", w = 170, h = 120, cap = cap(
    "Skill of SARIMA abundance forecasts against the cumulative average for each of the eight taxa, by estuary (north to south) and coastwide, ",
    "for each sampling method, 1 month (circles) and 12 months (triangles) ahead. Skill: proportion of the cumulative average's squared error ",
    "that the forecast removes (0 = no better); cumulative average: the seasonal curve plus the mean departure of all earlier months. Filled: ",
    "the 95% uncertainty range lies above zero (moving-block bootstrap over forecast dates); open: it includes zero. The line joins the ",
    "two values of a series; taxa are ordered by their median one-month skill. Scientific names are in Table S3; menhaden is the genus ",
    "*Brevoortia*. Calculated on the 367 forecast dates up to December 2022. One month ahead, bay trawl forecasts had higher skill than ",
    "bag seine forecasts for the same taxon and estuary in 49 of 56 comparisons. By 12 months, the whole 95% uncertainty range of the ",
    "skill was above zero in only 5 of the 56 bag seine series and 12 of the 56 bay trawl series; 8 of these 17 were blue crab"))
)

# -----------------------------------------------------------------------------
# 2) Image size from the TIFF header (width and height in pixels)
# -----------------------------------------------------------------------------
tiff_dims <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  bo <- rawToChar(readBin(con, "raw", 2)); endian <- if (bo == "II") "little" else if (bo == "MM") "big" else stop("Not a TIFF: ", path, call. = FALSE)
  rd <- function(n, size) readBin(con, "integer", n = n, size = size, signed = FALSE, endian = endian)
  rd(1, 2); off <- readBin(con, "integer", size = 4, endian = endian)
  seek(con, off); nt <- rd(1, 2); w <- h <- NA
  for (i in seq_len(nt)) {
    tag <- rd(1, 2); type <- rd(1, 2); readBin(con, "integer", size = 4, endian = endian)
    v <- if (type == 3) { x <- rd(1, 2); rd(1, 2); x } else readBin(con, "integer", size = 4, endian = endian)
    if (tag == 256) w <- v; if (tag == 257) h <- v
  }
  c(w = w, h = h)
}
place <- function(f, max_h = TEXT_H) {
  path <- file.path(fig_dir, paste0(f$file, ".tiff"))
  if (!file.exists(path)) stop("Missing image for Figure ", f$n, ": ", path, " (run its g script first)", call. = FALSE)
  d <- tiff_dims(path); asp <- d[["h"]] / d[["w"]]
  if (!is.na(f$w) && abs(asp - f$h / f$w) > 0.01) stop("Figure ", f$n, ": image aspect ", round(asp, 3), " differs from the listed canvas", call. = FALSE)
  w <- min(TEXT_W, max_h / asp)
  list(path = path, width = w, height = w * asp)
}

# -----------------------------------------------------------------------------
# 3) Word files
# -----------------------------------------------------------------------------
tnr <- function(size, bold = FALSE) fp_text(font.size = size, bold = bold, font.family = "Times New Roman")
# Captions may mark italic text with asterisks (*Brevoortia*): odd segments are set in italics.
cap_chunks <- function(txt) {
  seg <- strsplit(txt, "*", fixed = TRUE)[[1]]
  if (length(seg) %% 2 == 0) stop("Unpaired asterisk in a caption: ", substr(txt, 1, 60), call. = FALSE)
  Map(function(s, i) ftext(s, fp_text(font.size = 11, italic = i %% 2 == 0, font.family = "Times New Roman")), seg, seq_along(seg)) |> unname()
}
add_figure <- function(doc, f, max_h = TEXT_H) {
  im <- place(f, max_h)
  doc |> body_add_fpar(fpar(external_img(im$path, width = im$width, height = im$height)), style = "Normal") |>
    body_add_fpar(do.call(fpar, c(list(ftext(paste0("Fig. ", f$n, " "), tnr(11, TRUE))), cap_chunks(f$cap))))
}
sec <- prop_section(page_size = page_size(orient = "portrait", width = 8.5, height = 11),
                    page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1))
main <- Filter(function(f) !grepl("^S", f$n), FIGS); supp <- Filter(function(f) grepl("^S", f$n), FIGS)
written <- character(0)
for (f in main) {
  out <- file.path(out_dir, paste0("Figure_", f$n, ".docx"))
  if (file.exists(out)) archive_files(out)
  print(add_figure(read_docx() |> body_set_default_section(sec), f), target = out)
  written <- c(written, out)
}
out_or <- file.path(out_dir, "Online_Resource_Figures.docx")
if (file.exists(out_or)) archive_files(out_or)
doc <- read_docx() |> body_set_default_section(sec) |> body_add_fpar(fpar(ftext(OR_FIGURES_TITLE, tnr(14, TRUE))))
for (h in ARTICLE_HEADER) doc <- body_add_fpar(doc, fpar(ftext(h, tnr(11))))
for (i in seq_along(supp)) {
  doc <- add_figure(doc, supp[[i]], max_h = if (is.null(supp[[i]]$max_h)) OR_MAX_H else supp[[i]]$max_h)
  if (i < length(supp)) doc <- body_add_break(doc)
}
print(doc, target = out_or); written <- c(written, out_or)
for (f in written) cat("Saved:", f, "\n")
end_log("finished")
cat("g20_figure_docx.R completed\n")
