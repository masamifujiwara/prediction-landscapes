# =============================================================================
# Script:  h13_missing_months.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds Table S2 (Online Resource 2): months without sampling in the monthly
#   series used by the forecasts. For each data set (biodiversity, abundance,
#   water temperature and salinity), scale and sampling method: the number of
#   series, the first month, the number of months and the number of months
#   inside the record without sampling.
#
# Run after: b08, b09, c02, e03
#
# Inputs:
#   01_Analysis/Data/c2_deseasonalize_alpha_gamma_scaled/alpha_q012_deseas_scaled.rds,
#     gamma_q012_deseas_scaled.rds (c02)
#   01_Analysis/Data/e3_deseasonalize_abundance/alpha_abundance_deseas.rds,
#     gamma_abundance_deseas.rds (e03)
#   01_Analysis/Data/b8_build_environment_timeseries/env_alpha_index.rds,
#     env_gamma_index.rds (b08)
#   01_Analysis/Data/b9_forecast_settings/settings.rds (b09: gap months)
#
# Outputs:
#   01_Analysis/tables/Table_S2_missing_months.docx
#   01_Analysis/Data/h13_missing_months/missing_months_by_series.csv
#   01_Analysis/Data/h13_missing_months/missing_months_summary.csv
#   01_Analysis/Logs/h13_missing_months/run_<timestamp>.txt
#   Earlier CSV files and the earlier Word file are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   The facts come from the series the models use: the c02 and e03 departure
#   files (targets, with the `filled` flag) and the b08 covariate files.
#   The script stops if any series is not a complete monthly sequence ending
#   in December 2024, if the filled months of any series differ from the gap
#   months in the settings (April and May 2020), or if the number of series
#   differs from the design (48 biodiversity, 128 abundance, 16 environmental).
#   Series that begin later (for example bay trawl series in some estuaries,
#   and the coastwide bay trawl series, because a coastwide month needs every
#   estuary) are not extended backwards.
#   The note states how the gap was handled: targets filled by interpolating
#   the departure from each calendar month's mean (log scale for biodiversity,
#   log(CPUE + 1) for abundance; monthly means from months up to June 2020
#   only), covariates by linear interpolation; filled months are never
#   forecast targets or forecast dates.
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
need <- c("dplyr", "readr", "officer", "flextable")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(dplyr); library(officer); library(flextable)})

TABLE_LABEL <- "Table S2"   # Online Resource 2 number
script_stem <- "h13_missing_months"
D <- PATHS$data_dir
inputs <- c(div_a = file.path(D, "c2_deseasonalize_alpha_gamma_scaled", "alpha_q012_deseas_scaled.rds"),
            div_g = file.path(D, "c2_deseasonalize_alpha_gamma_scaled", "gamma_q012_deseas_scaled.rds"),
            abd_a = file.path(D, "e3_deseasonalize_abundance", "alpha_abundance_deseas.rds"),
            abd_g = file.path(D, "e3_deseasonalize_abundance", "gamma_abundance_deseas.rds"),
            env_a = file.path(D, "b8_build_environment_timeseries", "env_alpha_index.rds"),
            env_g = file.path(D, "b8_build_environment_timeseries", "env_gamma_index.rds"),
            settings = file.path(D, "b9_forecast_settings", "settings.rds"))
dirs <- script_dirs(script_stem, c("data", "logs"))
tab_dir <- file.path(PATHS$analysis_dir, "tables"); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
start_log(dirs$logs, script_stem)
stop_if_missing(inputs, c(rep("c02_deseasonalize_alpha_gamma_scaled.R", 2), rep("e03_deseasonalize_abundance.R", 2),
                          rep("b08_build_environment_timeseries.R", 2), "b09_forecast_settings.R"))
settings <- load_settings(inputs[["settings"]])
out_docx <- file.path(tab_dir, "Table_S2_missing_months.docx")
out_series <- file.path(dirs$data, "missing_months_by_series.csv"); out_sum <- file.path(dirs$data, "missing_months_summary.csv")
old <- list.files(dirs$data, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)

# -----------------------------------------------------------------------------
# 1) One row per series: first and last month, months, filled months
# -----------------------------------------------------------------------------
std <- function(x, data, scale, id) {
  x <- as.data.frame(x)
  if (!"date" %in% names(x)) x$date <- x$DATE
  if (!"BAY_ID" %in% names(x)) x$BAY_ID <- NA_integer_
  x %>% transmute(data = data, scale = scale, BAY_ID = as.integer(BAY_ID), GEAR, target = as.character(.data[[id]]), date = as.Date(date), filled)
}
s <- bind_rows(
  std(readRDS(inputs[["div_a"]]), "Biodiversity", "bay", "q"), std(readRDS(inputs[["div_g"]]), "Biodiversity", "coastwide", "q"),
  std(readRDS(inputs[["abd_a"]]), "Abundance", "bay", "SPECIES_CODE"), std(readRDS(inputs[["abd_g"]]), "Abundance", "coastwide", "SPECIES_CODE"),
  std(transform(readRDS(inputs[["env_a"]]), target = "TEMP and SAL"), "Water temperature and salinity", "bay", "target"),
  std(transform(readRDS(inputs[["env_g"]]), target = "TEMP and SAL"), "Water temperature and salinity", "coastwide", "target"))
by_series <- s %>% group_by(data, scale, BAY_ID, GEAR, target) %>% arrange(date, .by_group = TRUE) %>%
  summarise(first = min(date), last = max(date), months = n(),
            complete = all(diff(as.numeric(format(date, "%Y")) * 12 + as.numeric(format(date, "%m"))) == 1),
            filled_months = paste(format(date[filled], "%Y-%m"), collapse = ";"), filled = sum(filled), .groups = "drop")
gap <- paste(format(settings$gap_months, "%Y-%m"), collapse = ";")
if (any(by_series$last != as.Date("2024-12-01"))) stop("A series does not end in December 2024", call. = FALSE)
if (!all(by_series$complete)) stop("A series is not a complete monthly sequence", call. = FALSE)
bad <- by_series %>% filter(filled_months != gap)
if (nrow(bad)) { print(as.data.frame(bad)); stop("Filled months other than ", gap, " (listed above)", call. = FALSE) }
n_exp <- c(Biodiversity = 48L, Abundance = 128L, `Water temperature and salinity` = 16L)
n_got <- table(by_series$data)[names(n_exp)]
if (!identical(as.integer(n_got), unname(n_exp))) stop("Unexpected number of series: ", paste(names(n_got), n_got, collapse = ", "), call. = FALSE)
cat("Series:", nrow(by_series), "| every series complete; the only months without sampling inside any series:", gap, "\n")

# -----------------------------------------------------------------------------
# 2) Summary by data set, scale and gear
# -----------------------------------------------------------------------------
GEARS <- c(BS = "Bag seine", BT = "Bay trawl")
fmt_first <- function(first, scale, bay) {
  t <- table(format(first, "%b %Y"))
  if (length(t) == 1L) return(names(t))
  # more than one start: say how many estuaries start when
  bt <- tapply(bay, format(first, "%b %Y"), function(b) length(unique(b)))
  paste0(names(bt), " (", bt, " estuaries)", collapse = "; ")
}
summ <- by_series %>% group_by(data, scale, GEAR) %>%
  summarise(series = n(), first = fmt_first(first, scale[1], BAY_ID), last = paste(unique(format(last, "%b %Y")), collapse = "; "),
            months = paste(sort(unique(months), decreasing = TRUE), collapse = " or "),
            not_sampled = paste0(unique(filled), " (April and May 2020)"), not_sampled_n = unique(filled), .groups = "drop") %>%
  mutate(data = factor(data, names(n_exp)), scale = factor(scale, c("bay", "coastwide"))) %>% arrange(data, scale, GEAR)
print(as.data.frame(summ), row.names = FALSE)

tab <- summ %>% transmute(`Data set` = as.character(data), Scale = ifelse(scale == "bay", "Estuaries", "Coastwide"), `Sampling method` = GEARS[GEAR],
                          `Number of series` = as.character(series), `First month` = gsub("; ", ";\n", first, fixed = TRUE), `Number of months` = months,   # one start per line in the table
                          `Number of months not sampled (April and May 2020)` = as.character(summ$not_sampled_n))
ft <- flextable(as.data.frame(tab)) %>% merge_v(j = 1:2) %>% valign(j = 1:2, valign = "top", part = "body") %>%
  font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 8.5, part = "all") %>% bold(part = "header") %>%
  align(j = 4:ncol(tab), align = "center", part = "all") %>% border_remove() %>%
  hline_top(border = fp_border(width = 1.5), part = "header") %>% hline_bottom(border = fp_border(width = 0.75), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  padding(padding.top = 1, padding.bottom = 1, padding.left = 3, padding.right = 3, part = "all") %>%
  width(width = c(1.3, 0.8, 0.75, 0.55, 1.5, 0.75, 0.85))
caption <- paste0(
  "Months without sampling in the monthly series used by the forecasts. Biodiversity: Hill numbers q = 0, 1, 2; ",
  "abundance: log(CPUE + 1) of eight taxa, CPUE being catch per unit effort; water temperature and salinity: monthly means used as covariates. Estuaries: seven ",
  "estuaries; coastwide: all estuaries pooled, a month counting only when every estuary was sampled. All series end in December 2024.")
note <- paste0(
  "No samples were taken in April and May 2020, the only months without sampling inside any series. For model fitting, ",
  "these two months were filled: biodiversity and abundance by interpolating the departure from each calendar month's mean ",
  "(log scale for biodiversity, log(CPUE + 1) for abundance; monthly means from the months up to June 2020 only), water ",
  "temperature and salinity by linear interpolation. The filled months were never used as forecast targets or forecast ",
  "dates, so they enter no accuracy or skill measure. Bay trawl sampling in Matagorda Bay, Corpus Christi Bay and the two ",
  "Laguna Madre estuaries began in May 1982; these series, and the coastwide bay trawl series, start then and are not ",
  "extended backwards. Forecast dates run from April 1992 to November 2024.")
if (file.exists(out_docx)) archive_files(out_docx)
doc <- read_docx() %>% body_add_fpar(fpar(ftext(paste0(TABLE_LABEL, " "), fp_text(font.size = 10, bold = TRUE, font.family = "Times New Roman")),
                                              ftext(caption, fp_text(font.size = 10, font.family = "Times New Roman")))) %>%
  body_add_flextable(ft) %>% body_add_fpar(fpar(ftext(note, fp_text(font.size = 9, font.family = "Times New Roman"))))
print(doc, target = out_docx)
readr::write_csv(by_series %>% select(-complete), out_series); readr::write_csv(summ, out_sum)
for (f in c(out_docx, out_series, out_sum)) cat("Saved:", f, "\n")
end_log("finished")
cat("h13_missing_months.R completed\n")
