# =============================================================================
# Script:  h17_environment_trends.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Builds Table S9 (Online Resource 2): Mann-Kendall test and Sen slope for
#   water temperature and salinity in every estuary and sampling method and
#   coastwide, the two variables side by side.
#
# Run after: f04
#
# Inputs:
#   01_Analysis/Data/f4_environment_trends/f4_environment_trends.rds (f04)
#
# Outputs:
#   01_Analysis/tables/Table_S9_environment_trends.docx
#   01_Analysis/Data/h17_environment_trends/environment_trends_table.csv
#   01_Analysis/Logs/h17_environment_trends/run_<timestamp>.txt
#   Earlier CSV files and the earlier Word file are moved to ARCHIVE_<timestamp>.
#
# Method notes:
#   The values are not recalculated: the table shows the f04 results
#   (Mann-Kendall test and Sen slope on the yearly means of the departures
#   from the seasonal curve, years with at least six sampled months;
#   Benjamini-Hochberg adjustment over the seven estuaries within each
#   variable and sampling method). The script stops if f04 does not give
#   exactly 16 series per variable, or if the series differ in years.
#   Columns: estuary, sampling method, then for each variable Kendall's tau
#   (header "Tau", defined in the note), the Sen slope per decade, p and
#   adjusted p. The coastwide series are single tests and have no adjusted p;
#   the cell is left empty and the note says so. p below 0.001 is shown as
#   "< 0.001". The prewhitened check of f04 (Yue et al. 2002) is not shown.
#   A test copy with the stem ending in _TEST writes a Word file ending in
#   _TEST. h20 joins this table with the others.
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

TABLE_LABEL <- "Table S9"   # Online Resource 2 number
script_stem <- "h17_environment_trends"    # a test copy uses its own stem, e.g. "h17_environment_trends_TEST"
in_f04 <- file.path(PATHS$data_dir, "f4_environment_trends", "f4_environment_trends.rds")
dirs <- script_dirs(script_stem, c("data", "logs"))
tab_dir <- file.path(PATHS$analysis_dir, "tables"); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
start_log(dirs$logs, script_stem)
stop_if_missing(in_f04, "f04_environment_trends.R")
out_docx <- file.path(tab_dir, paste0(sub(" ", "_", TABLE_LABEL), "_environment_trends",
                                      if (grepl("_TEST$", script_stem)) "_TEST" else "", ".docx"))
out_csv <- file.path(dirs$data, "environment_trends_table.csv")
old <- list.files(dirs$data, pattern = "[.]csv$", full.names = TRUE); if (length(old)) archive_files(old)

# -----------------------------------------------------------------------------
# 1) f04 results, checked
# -----------------------------------------------------------------------------
tr <- as.data.frame(readRDS(in_f04))
expect_cols(tr, c("variable", "scale", "BAY_ID", "area", "GEAR", "n_years", "mk_tau", "mk_p", "sen_per_decade", "p_bh"), "f4_environment_trends")
n_per <- table(tr$variable)
if (!setequal(names(n_per), c("temperature", "salinity")) || any(n_per != 16))
  stop("Expected 16 series for each of temperature and salinity, found: ", paste(names(n_per), n_per, collapse = ", "), call. = FALSE)
if (length(unique(tr$n_years)) != 1) stop("The series differ in the number of years", call. = FALSE)
n_years <- unique(tr$n_years)
cat("Series:", nrow(tr), "| years per series:", n_years, "\n")

fmt2 <- function(x) { s <- formatC(x, format = "f", digits = 2); ifelse(s == "-0.00", "0.00", s) }   # no "-0.00"
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3)))
one <- function(v) tr %>% filter(variable == v) %>%
  transmute(BAY_ID, GEAR, tau = fmt2(mk_tau), sen = fmt2(sen_per_decade),
            p = fmt_p(mk_p), p_adj = fmt_p(ifelse(scale == "coastwide", NA_real_, p_bh)))
keys <- tr %>% filter(variable == "temperature") %>% distinct(scale, BAY_ID, area, GEAR) %>%
  mutate(o = ifelse(scale == "coastwide", 99L, BAY_ID)) %>% arrange(o, GEAR)
tab <- keys %>% left_join(one("temperature"), by = c("BAY_ID", "GEAR")) %>%
  left_join(one("salinity"), by = c("BAY_ID", "GEAR"), suffix = c("_t", "_s")) %>%
  transmute(Estuary = area, `Sampling method` = ifelse(GEAR == "BS", "Bag seine", "Bay trawl"),
            tau_t, sen_t, p_t, p_adj_t, tau_s, sen_s, p_s, p_adj_s)
if (anyNA(tab) || nrow(tab) != 16) stop("Temperature and salinity rows do not match", call. = FALSE)
readr::write_csv(tab, out_csv)

# -----------------------------------------------------------------------------
# 2) Word table
# -----------------------------------------------------------------------------
ft <- flextable(tab) %>%
  set_header_labels(tau_t = "Tau", sen_t = "Sen slope (\u00b0C per decade)", p_t = "p", p_adj_t = "Adjusted p",
                    tau_s = "Tau", sen_s = "Sen slope (per decade)", p_s = "p", p_adj_s = "Adjusted p") %>%
  add_header_row(values = c("", "Water temperature", "Salinity"), colwidths = c(2, 4, 4)) %>%
  merge_v(j = 1) %>% valign(j = 1, valign = "top", part = "body") %>%
  font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 8.5, part = "all") %>% bold(part = "header") %>%
  align(j = 3:10, align = "center", part = "all") %>% align(i = 1, align = "center", part = "header") %>% border_remove() %>%
  hline_top(border = fp_border(width = 1.5), part = "header") %>% hline(i = 1, j = 3:10, border = fp_border(width = 0.5), part = "header") %>%
  hline_bottom(border = fp_border(width = 0.75), part = "header") %>% hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  width(width = c(1.2, 0.7, rep(c(0.6, 0.75, 0.5, 0.6), 2)))
caption <- paste0(
  "Trends in water temperature and salinity from 1982 to 2024 in each estuary and sampling method and coastwide. ",
  "Mann-Kendall test and Sen slope on the yearly means of the departures from the seasonal curve (", n_years,
  " years in every series, each with at least six sampled months).")
note <- paste0(
  "Tau: Kendall's tau of the Mann-Kendall test. Adjusted p: p adjusted for multiple tests over the seven estuaries within each variable and sampling method ",
  "(Benjamini and Hochberg 1995). The coastwide series are single tests and have no adjusted p. ",
  "Salinity has no unit. Coastwide: effort-weighted monthly means over all estuaries (Online Resource 1, Sect. 1.3).")
if (file.exists(out_docx)) archive_files(out_docx)
doc <- read_docx() %>% body_add_fpar(fpar(ftext(paste0(TABLE_LABEL, " "), fp_text(font.size = 10, bold = TRUE, font.family = "Times New Roman")),
                                              ftext(caption, fp_text(font.size = 10, font.family = "Times New Roman")))) %>%
  body_add_flextable(ft) %>% body_add_fpar(fpar(ftext(note, fp_text(font.size = 9, font.family = "Times New Roman"))))
print(doc, target = out_docx)
print(tab, row.names = FALSE)
for (f in c(out_docx, out_csv)) cat("Saved:", f, "\n")
end_log("finished")
cat("h17_environment_trends.R completed\n")
