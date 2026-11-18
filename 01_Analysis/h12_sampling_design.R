# =============================================================================
# Script:  h12_sampling_design.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds Table S1 (Online Resource 2): sampling design and community
#   biodiversity by estuary and sampling method. Columns: years sampled,
#   sampling months, station visits, mean stations per month; median and
#   interquartile range of the monthly Hill numbers (q = 0, 1, 2); Cohen's d
#   of each estuary against the other estuaries of the same sampling method,
#   with a bootstrap 95% range.
#
# Run after: b03, b06, b09
#
# Inputs:
#   01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds (b03)
#   01_Analysis/Data/b6_compute_hill_alpha_q012/alpha_q012.rds (b06)
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09: estuaries, seed)
#
# Outputs:
#   01_Analysis/tables/Table_S1_sampling_design.docx
#   01_Analysis/Data/h12_sampling_design/sampling_design_by_bay_gear.csv
#   01_Analysis/Data/h12_sampling_design/diversity_median_cohens_d_by_bay_gear.csv
#   01_Analysis/Logs/h12_sampling_design/run_<timestamp>.txt
#   Earlier CSV files and the earlier Word file are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   The seven estuaries come from the settings (Sabine Lake is dropped in a00,
#   so it is in neither the estuary nor the coastwide series).
#   Biodiversity values are those used in the analysis: the b06 Hill numbers
#   for each estuary and month at the common sample coverage of 0.90, leaving
#   out cells with no catch or failed coverage.
#   Cohen's d uses the pooled SD; its 95% range comes from 2,000 bootstrap
#   resamples with the project seed (settings$seed). The bootstrap treats
#   months as independent, so the range understates the uncertainty; the
#   caption says so.
#   Years are written in plain text ("1982 to 2024"), without dashes.
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
source(file.path(PATHS$analysis_dir, "lib", "forecast_rules.R"))
need <- c("dplyr", "tidyr", "readr", "officer", "flextable")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(officer); library(flextable)})

# -----------------------------------------------------------------------------
# 1) Settings, folders, log
# -----------------------------------------------------------------------------
TABLE_LABEL <- "Table S1"   # Online Resource 2 number
B_BOOT <- 2000L
script_stem <- "h12_sampling_design"
in_station  <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds")
in_alpha    <- file.path(PATHS$data_dir, "b6_compute_hill_alpha_q012", "alpha_q012.rds")
in_settings <- file.path(PATHS$data_dir, "b9_forecast_settings", "settings.rds")
dirs <- script_dirs(script_stem, c("data", "logs"))
tab_dir <- file.path(PATHS$analysis_dir, "tables"); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
start_log(dirs$logs, script_stem)
stop_if_missing(c(in_station, in_alpha, in_settings),
                c("b03_harmonize_station_metadata.R", "b06_compute_hill_alpha_q012.R", "b09_forecast_settings.R"))
settings <- load_settings(in_settings)
BAYS <- as.character(settings$bays)
BAY_NAMES <- c("2" = "Galveston Bay", "3" = "Matagorda Bay", "4" = "San Antonio Bay", "5" = "Aransas Bay",
               "6" = "Corpus Christi Bay", "7" = "Upper Laguna Madre", "8" = "Lower Laguna Madre")
if (!all(BAYS %in% names(BAY_NAMES))) stop("A bay in the settings has no name in BAY_NAMES", call. = FALSE)
GEARS <- c(BS = "Bag seine", BT = "Bay trawl")
out_docx <- file.path(tab_dir, "Table_S1_sampling_design.docx")
out_design <- file.path(dirs$data, "sampling_design_by_bay_gear.csv")
out_div <- file.path(dirs$data, "diversity_median_cohens_d_by_bay_gear.csv")
old <- list.files(dirs$data, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)

# -----------------------------------------------------------------------------
# 2) Sampling design
# -----------------------------------------------------------------------------
station <- readRDS(in_station)
design_one <- function(df, g) {
  df %>% mutate(MAJOR_AREA = as.character(MAJOR_AREA), YM = YEAR_i * 100L + MONTH_i) %>% group_by(MAJOR_AREA) %>%
    summarise(GEAR = g, year_start = min(YEAR_i), year_end = max(YEAR_i), months = n_distinct(YM),
              visits = n(), per_month = n() / n_distinct(YM), .groups = "drop")
}
design <- bind_rows(design_one(station$station_bs_core, "BS"), design_one(station$station_bt_core, "BT"))
if (!setequal(design$MAJOR_AREA, BAYS) || nrow(design) != 2L * length(BAYS)) {
  stop("Station data do not cover exactly the settings' estuaries for both sampling methods", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 3) Diversity: median and IQR per bay and gear; Cohen's d against the other bays
# -----------------------------------------------------------------------------
alpha <- readRDS(in_alpha) %>%
  filter(!is.na(q0), !is.na(q1), !is.na(q2), N_total > 0, !coverage_fail) %>%
  mutate(BAY_ID = as.character(BAY_ID), GEAR = as.character(GEAR))
if (!setequal(unique(alpha$BAY_ID), BAYS)) stop("Biodiversity cells do not cover exactly the settings' estuaries", call. = FALSE)
div <- alpha %>% group_by(BAY_ID, GEAR) %>%
  summarise(cells = n(), across(c(q0, q1, q2), list(med = stats::median, iqr = stats::IQR)), .groups = "drop")

cohens_d_boot <- function(x, y, B) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  psd <- function(a, b) sqrt(((length(a) - 1) * stats::var(a) + (length(b) - 1) * stats::var(b)) / (length(a) + length(b) - 2))
  d <- (mean(x) - mean(y)) / psd(x, y)
  bd <- vapply(seq_len(B), function(i) { a <- sample(x, replace = TRUE); b <- sample(y, replace = TRUE); (mean(a) - mean(b)) / psd(a, b) }, numeric(1))
  c(d = d, lo = unname(stats::quantile(bd, 0.025)), hi = unname(stats::quantile(bd, 0.975)))
}
set.seed(settings$seed)
eff <- do.call(rbind, lapply(names(GEARS), function(g) do.call(rbind, lapply(BAYS, function(b) {
  a <- alpha %>% filter(GEAR == g)
  r <- lapply(c("q0", "q1", "q2"), function(q) cohens_d_boot(a[[q]][a$BAY_ID == b], a[[q]][a$BAY_ID != b], B_BOOT))
  data.frame(BAY_ID = b, GEAR = g, d_q0 = r[[1]][["d"]], lo_q0 = r[[1]][["lo"]], hi_q0 = r[[1]][["hi"]],
             d_q1 = r[[2]][["d"]], lo_q1 = r[[2]][["lo"]], hi_q1 = r[[2]][["hi"]],
             d_q2 = r[[3]][["d"]], lo_q2 = r[[3]][["lo"]], hi_q2 = r[[3]][["hi"]])
}))))
div_eff <- div %>% left_join(eff, by = c("BAY_ID", "GEAR")) %>% mutate(Estuary = BAY_NAMES[BAY_ID]) %>% arrange(BAY_ID, GEAR)
if (nrow(div_eff) != 2L * length(BAYS) || anyNA(div_eff)) stop("Biodiversity summary incomplete", call. = FALSE)

# -----------------------------------------------------------------------------
# 4) Table
# -----------------------------------------------------------------------------
f2 <- function(x) sub("^-0[.]00$", "0.00", sprintf("%.2f", x))
tab <- design %>% left_join(div_eff, by = c("MAJOR_AREA" = "BAY_ID", "GEAR")) %>%
  mutate(order = match(MAJOR_AREA, BAYS)) %>% arrange(order, GEAR) %>%
  transmute(Estuary = BAY_NAMES[MAJOR_AREA], Gear = GEARS[GEAR], Years = paste(year_start, "to", year_end),
            `Sampling months` = as.character(months), `Station visits` = formatC(visits, format = "d", big.mark = ","),
            `Stations per month` = sprintf("%.1f", per_month),
            `q = 0` = paste0(f2(q0_med), " (", f2(q0_iqr), ")"), `q = 1` = paste0(f2(q1_med), " (", f2(q1_iqr), ")"),
            `q = 2` = paste0(f2(q2_med), " (", f2(q2_iqr), ")"),
            `d, q = 0` = paste0(f2(d_q0), " [", f2(lo_q0), ", ", f2(hi_q0), "]"),
            `d, q = 1` = paste0(f2(d_q1), " [", f2(lo_q1), ", ", f2(hi_q1), "]"),
            `d, q = 2` = paste0(f2(d_q2), " [", f2(lo_q2), ", ", f2(hi_q2), "]"))
print(as.data.frame(tab), row.names = FALSE)
n_col <- ncol(tab)
ft <- flextable(as.data.frame(tab) |> setNames(sub("^Gear$", "Sampling method", names(tab)))) %>%
  add_header_row(values = c("", "Sampling design", "Median biodiversity (IQR)", "Cohen's d against the other estuaries [95% range]"),
                 colwidths = c(2, 4, 3, 3)) %>%
  width(j = seq_len(n_col), width = c(1.10, 0.55, 0.75, 0.55, 0.55, 0.55, 0.72, 0.72, 0.72, 1.05, 1.05, 1.05)) %>%   # 9.4 inches; page margins 0.75 inch
  font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 7.5, part = "body") %>% fontsize(size = 8, part = "header") %>%
  bold(part = "header") %>% align(j = 1:2, align = "left", part = "all") %>% align(j = 3:n_col, align = "center", part = "all") %>%
  merge_v(j = 1) %>% valign(j = 1, valign = "top", part = "body") %>% border_remove() %>%
  hline_top(border = fp_border(width = 1.5), part = "header") %>% hline(i = 1, border = fp_border(width = 0.5, style = "dotted"), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.0), part = "header") %>% hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  height_all(height = 0.20) %>% hrule(rule = "atleast") %>%
  padding(padding.top = 1, padding.bottom = 1, padding.left = 3, padding.right = 3, part = "all") %>% line_spacing(space = 1, part = "all")
### Cohen's d in italics (statistical symbol).
ft <- compose(ft, part = "header", i = 1, j = n_col - 2, value = as_paragraph("Cohen's ", as_i("d"), " against the other estuaries [95% range]"))
for (k in 0:2) ft <- compose(ft, part = "header", i = 2, j = n_col - 2 + k, value = as_paragraph(as_i("d"), paste0(", q = ", k)))
ft <- ft %>% font(fontname = "Times New Roman", part = "header") %>% fontsize(size = 8, part = "header") %>% bold(part = "header")
caption_txt <- paste0(
  "Sampling design and community biodiversity of the seven estuaries analyzed, Texas Parks and Wildlife ",
  "Department Coastal Fisheries Marine Resource Monitoring Program. Sampling months: calendar months with at least one ",
  "haul or tow. Station visits: all hauls or tows over the record. Stations per month: mean number of hauls or tows per ",
  "sampled month. Median biodiversity (interquartile range): Hill numbers of each estuary, sampling method, and month at the common ",
  "sample coverage of 0.90 (orders q = 0, q = 1 and q = 2). Cohen's ")
caption_d <- "d"
caption_end <- paste0(
  " compares each estuary's monthly values with those of the other six estuaries for the same sampling method (positive: higher ",
  "biodiversity); the 95% range is from ", formatC(B_BOOT, format = "d", big.mark = ","), " bootstrap resamples of months ",
  "treated as independent, so it understates the uncertainty and is indicative only.")
footnote_txt <- paste0(
  "Sabine Lake was excluded from all analyses, at the estuary and the coastwide scale, and is not shown. No samples were ",
  "taken in April and May 2020. Bay trawl sampling in Matagorda Bay, Corpus Christi Bay and the two Laguna Madre ",
  "estuaries began in May 1982.")
landscape <- block_section(prop_section(page_size = page_size(orient = "landscape"),
  page_margins = page_mar(top = 0.75, bottom = 0.75, left = 0.75, right = 0.75), type = "continuous"))
if (file.exists(out_docx)) archive_files(out_docx)
doc <- read_docx() %>%
  body_add_fpar(fpar(ftext(paste0(TABLE_LABEL, " "), fp_text(font.size = 10, bold = TRUE, font.family = "Times New Roman")),
                     ftext(caption_txt, fp_text(font.size = 10, font.family = "Times New Roman")),
                     ftext(caption_d, fp_text(font.size = 10, italic = TRUE, font.family = "Times New Roman")),
                     ftext(caption_end, fp_text(font.size = 10, font.family = "Times New Roman")))) %>%
  body_add_flextable(ft) %>%
  body_add_fpar(fpar(ftext(footnote_txt, fp_text(font.size = 9, font.family = "Times New Roman")))) %>%
  body_end_block_section(value = landscape)
print(doc, target = out_docx)
readr::write_csv(design %>% mutate(per_month = round(per_month, 3)), out_design)
readr::write_csv(div_eff %>% mutate(across(where(is.double), ~ round(.x, 4))), out_div)
for (f in c(out_docx, out_design, out_div)) cat("Saved:", f, "\n")
end_log("finished")
cat("h12_sampling_design.R completed\n")
