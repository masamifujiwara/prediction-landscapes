# =============================================================================
# Script: c1_inputs_check.R
#
# Purpose:
#   Load ALL c1 inputs and validate they are usable:
#     - gamma_q012.RDS  (wide q0/q1/q2 OR long q/value)
#     - alpha_q012.RDS  (bay-level long; may be list-based)
#   For each:
#     - confirm required columns exist
#     - check month range
#     - detect duplicate cells (GEAR,q,YEAR,MONTH) after coercion to long
#     - report NA rate
#     - write CSV diagnostics (+ basic plots when ggplot2 available)
#
# Inputs (per project conventions):
#   - 01_Analysis/Data/01_TimeSeries/gamma_q012.RDS
#   - 01_Analysis/Data/b6_compute_hill_alpha_q012/alpha_q012.rds
#
# Outputs (per project conventions):
#   Derived data:
#     - 01_Analysis/Data/c1_inputs_check/coerced_long/gamma_q012_long.rds
#     - 01_Analysis/Data/c1_inputs_check/coerced_long/alpha_q012_long.rds (if alpha exists)
#   Diagnostics:
#     - 02_Results/03_Diagnostics/c1_inputs_check/*.csv
#     - 02_Results/03_Diagnostics/c1_inputs_check/*.pdf + *.tiff (plots; if ggplot2 available)
# =============================================================================

rm(list = ls())

# -----------------------------
# Required bootstrap snippet (must include exactly)
# -----------------------------
get_script_dir <- function() {
  # 1) RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  
  # 2) Rscript / R CMD BATCH: --file=...
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  
  # 3) source("script.R")
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  
  # 4) fallback: assume you launched from 01_Analysis
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
setwd(script_dir)

# -----------------------------
# PATHS.RData usage (required)
# -----------------------------
project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file) # loads PATHS
stopifnot(exists("PATHS"))

# Validate PATHS project_dir matches computed project_dir
if (!is.null(PATHS$project_dir)) {
  paths_proj <- normalizePath(PATHS$project_dir, winslash = "/", mustWork = TRUE)
  if (!identical(paths_proj, project_dir)) {
    stop(
      "PATHS$project_dir does not match computed project_dir.\n",
      "Computed: ", project_dir, "\n",
      "PATHS:    ", paths_proj
    )
  }
}

# -----------------------------
# Minimal dependencies
# -----------------------------
need_pkgs <- c("dplyr", "tidyr", "readr", "stringr")
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

# -----------------------------
# Script stem + directory convention
# -----------------------------
script_stem <- tools::file_path_sans_ext(basename(sys.frames()[[1]]$ofile %||% "c1_inputs_check.R"))
# If sys.frames()[[1]]$ofile is NULL under some launch patterns, keep stable stem:
if (is.null(script_stem) || !nzchar(script_stem) || script_stem == "NULL") script_stem <- "c1_inputs_check"

# Derived data directory
dir_data_script <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_data_coerced <- file.path(dir_data_script, "coerced_long")
dir.create(dir_data_coerced, recursive = TRUE, showWarnings = FALSE)

# Logs directory (optional; not used heavily here but kept consistent)
dir_logs_script <- file.path(PATHS$analysis_dir, "Logs", script_stem)
dir.create(dir_logs_script, recursive = TRUE, showWarnings = FALSE)

# Diagnostics directory
dir_diag_script <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)
dir.create(dir_diag_script, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Inputs (authoritative)
# -----------------------------
IN_GAMMA <- file.path(PATHS$analysis_dir, "Data", "01_TimeSeries", "gamma_q012.RDS")
IN_ALPHA <- file.path(PATHS$analysis_dir, "Data", "b6_compute_hill_alpha_q012", "alpha_q012.rds")

# -----------------------------
# Utility functions
# -----------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

expect_cols <- function(df, cols, name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(name, " missing required columns: ", paste(miss, collapse = ", "))
}

safe_write_csv <- function(x, path) {
  readr::write_csv(x, path)
  invisible(path)
}

save_plot_pdf_tiff <- function(p, file_stem, width = 10, height = 6) {
  if (!has_ggplot2) return(invisible(NULL))
  pdf_path  <- file.path(dir_diag_script, paste0(file_stem, ".pdf"))
  tiff_path <- file.path(dir_diag_script, paste0(file_stem, ".tiff"))
  ggplot2::ggsave(filename = pdf_path, plot = p, width = width, height = height, units = "in", dpi = 300)
  ggplot2::ggsave(filename = tiff_path, plot = p, width = width, height = height, units = "in", dpi = 300, compression = "lzw")
  invisible(list(pdf = pdf_path, tiff = tiff_path))
}

# Coerce ANY (gamma/alpha) object into long: GEAR,q,YEAR,MONTH,value (+ optional BAY_ID)
coerce_to_long_ts <- function(obj, source_name = "input") {
  
  df <- obj
  if (is.list(obj) && !is.data.frame(obj)) {
    idx <- which(vapply(obj, is.data.frame, logical(1)))
    if (length(idx) < 1) stop("Input is a list but contains no data.frame: ", source_name)
    df <- obj[[idx[1]]]
  }
  stopifnot(is.data.frame(df))
  
  # Standardize common name variants
  if (!("GEAR" %in% names(df))) {
    gear_like <- intersect(names(df), c("gear", "Gear"))
    if (length(gear_like) >= 1) df <- df %>% rename(GEAR = all_of(gear_like[1]))
  }
  if (!("YEAR" %in% names(df))) {
    yr_like <- intersect(names(df), c("year", "Year"))
    if (length(yr_like) >= 1) df <- df %>% rename(YEAR = all_of(yr_like[1]))
  }
  if (!("MONTH" %in% names(df))) {
    mo_like <- intersect(names(df), c("month", "Month"))
    if (length(mo_like) >= 1) df <- df %>% rename(MONTH = all_of(mo_like[1]))
  }
  
  # If already long: q + value
  if (all(c("GEAR", "q", "YEAR", "MONTH", "value") %in% names(df))) {
    out <- df %>%
      mutate(
        GEAR  = as.character(GEAR),
        q     = suppressWarnings(as.numeric(q)),
        YEAR  = suppressWarnings(as.integer(YEAR)),
        MONTH = suppressWarnings(as.integer(MONTH)),
        value = suppressWarnings(as.numeric(value))
      )
    if ("BAY_ID" %in% names(out)) out <- out %>% mutate(BAY_ID = suppressWarnings(as.integer(BAY_ID)))
    return(out %>% select(any_of(c("BAY_ID", "GEAR", "q", "YEAR", "MONTH", "value"))))
  }
  
  # If long but q named differently (iNEXT-ish)
  q_alt <- intersect(names(df), c("Order.q", "Order", "order", "q_d", "q.d"))
  if (length(q_alt) == 1 && "value" %in% names(df)) {
    out <- df %>%
      rename(q = all_of(q_alt)) %>%
      mutate(
        GEAR  = as.character(GEAR),
        q     = suppressWarnings(as.numeric(q)),
        YEAR  = suppressWarnings(as.integer(YEAR)),
        MONTH = suppressWarnings(as.integer(MONTH)),
        value = suppressWarnings(as.numeric(value))
      )
    if ("BAY_ID" %in% names(out)) out <- out %>% mutate(BAY_ID = suppressWarnings(as.integer(BAY_ID)))
    return(out %>% select(any_of(c("BAY_ID", "GEAR", "q", "YEAR", "MONTH", "value"))))
  }
  
  # Otherwise: WIDE q0/q1/q2
  q_cols <- intersect(names(df), c("q0", "q1", "q2"))
  if (!all(c("GEAR", "YEAR", "MONTH") %in% names(df)) || length(q_cols) == 0) {
    stop(
      "Could not coerce input to long format.\nSource: ", source_name, "\nColumns:\n  ",
      paste(names(df), collapse = ", ")
    )
  }
  
  long <- df %>%
    mutate(
      GEAR  = as.character(GEAR),
      YEAR  = suppressWarnings(as.integer(YEAR)),
      MONTH = suppressWarnings(as.integer(MONTH))
    ) %>%
    pivot_longer(all_of(q_cols), names_to = "q_raw", values_to = "value") %>%
    mutate(
      q = suppressWarnings(as.numeric(str_extract(q_raw, "[012]"))),
      value = suppressWarnings(as.numeric(value))
    )
  
  # If BAY_ID exists in wide table, keep it
  if ("BAY_ID" %in% names(df)) long <- long %>% mutate(BAY_ID = suppressWarnings(as.integer(.data$BAY_ID)))
  
  long %>% select(any_of(c("BAY_ID", "GEAR", "q", "YEAR", "MONTH", "value")))
}

# Summaries + write diagnostics
check_one_input <- function(path, label) {
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Checking: ", label, "\n", sep = "")
  cat("Path: ", normalizePath(path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  
  if (!file.exists(path)) stop("Missing file: ", path)
  
  raw_obj <- readRDS(path)
  ts_long <- coerce_to_long_ts(raw_obj, source_name = basename(path))
  
  expect_cols(ts_long, c("GEAR", "q", "YEAR", "MONTH", "value"), label)
  
  # Basic summary
  basic <- ts_long %>%
    summarise(
      n_rows = n(),
      n_gears = n_distinct(GEAR),
      n_q = n_distinct(q),
      min_year = min(YEAR, na.rm = TRUE),
      max_year = max(YEAR, na.rm = TRUE),
      min_month = min(MONTH, na.rm = TRUE),
      max_month = max(MONTH, na.rm = TRUE),
      na_rate_value = mean(is.na(value)),
      n_infinite_value = sum(is.infinite(value), na.rm = TRUE)
    )
  
  # Month range check
  month_ok <- !(basic$min_month < 1 || basic$max_month > 12)
  
  # Duplicate cells AFTER coercion
  dup_cells <- ts_long %>%
    count(GEAR, q, YEAR, MONTH, name = "n") %>%
    filter(n > 1) %>%
    arrange(desc(n), GEAR, q, YEAR, MONTH)
  
  dup_summary <- dup_cells %>%
    summarise(
      n_dup_cells = n(),
      max_n = ifelse(n() == 0, 1, max(n)),
      median_n = ifelse(n() == 0, 1, median(n)),
      prop_cells_duplicated = n_dup_cells / nrow(distinct(ts_long, GEAR, q, YEAR, MONTH))
    )
  
  # Missing month diagnostics (per GEAR×q)
  miss_months <- ts_long %>%
    mutate(date = as.Date(sprintf("%04d-%02d-01", YEAR, MONTH))) %>%
    group_by(GEAR, q) %>%
    summarise(
      min_date = min(date, na.rm = TRUE),
      max_date = max(date, na.rm = TRUE),
      n_unique_months = n_distinct(date),
      n_expected = as.integer(
        (12 * (as.integer(format(max_date, "%Y")) - as.integer(format(min_date, "%Y")))) +
          (as.integer(format(max_date, "%m")) - as.integer(format(min_date, "%m"))) + 1
      ),
      n_missing = n_expected - n_unique_months,
      .groups = "drop"
    )
  
  # NA-rate summary
  na_by_group <- ts_long %>%
    mutate(date = as.Date(sprintf("%04d-%02d-01", YEAR, MONTH))) %>%
    group_by(GEAR, q) %>%
    summarise(
      n_rows = n(),
      n_unique_months = n_distinct(date),
      na_rate_value = mean(is.na(value)),
      .groups = "drop"
    )
  
  # -----------------------------
  # Write diagnostics (CSV)
  # -----------------------------
  out_basic <- file.path(dir_diag_script, paste0("c1_inputcheck__", label, "__basic_summary.csv"))
  safe_write_csv(basic, out_basic)
  
  out_dup_summary <- file.path(dir_diag_script, paste0("c1_inputcheck__", label, "__dup_summary.csv"))
  safe_write_csv(dup_summary, out_dup_summary)
  
  out_miss <- file.path(dir_diag_script, paste0("c1_inputcheck__", label, "__missing_months_summary.csv"))
  safe_write_csv(miss_months, out_miss)
  
  out_na <- file.path(dir_diag_script, paste0("c1_inputcheck__", label, "__na_rates.csv"))
  safe_write_csv(na_by_group, out_na)
  
  if (nrow(dup_cells) > 0) {
    out_dup <- file.path(dir_diag_script, paste0("c1_inputcheck__", label, "__dup_cells.csv"))
    safe_write_csv(dup_cells, out_dup)
    
    # Show which columns vary within one duplicated cell (in the coerced long table)
    ex <- dup_cells[1, ]
    ex_rows <- ts_long %>%
      filter(GEAR == ex$GEAR, q == ex$q, YEAR == ex$YEAR, MONTH == ex$MONTH)
    
    varying <- ex_rows %>%
      summarise(across(everything(), ~ dplyr::n_distinct(.x))) %>%
      pivot_longer(everything(), names_to = "col", values_to = "n_distinct") %>%
      filter(n_distinct > 1) %>%
      arrange(desc(n_distinct))
    
    out_var <- file.path(dir_diag_script, paste0("c1_inputcheck__", label, "__varying_cols_example.csv"))
    safe_write_csv(varying, out_var)
  }
  
  # -----------------------------
  # Save derived (coerced) long version
  # -----------------------------
  out_long_rds <- file.path(dir_data_coerced, paste0(label, "_long.rds"))
  saveRDS(ts_long, out_long_rds)
  
  # -----------------------------
  # Plots (optional; ggplot2 only)
  # -----------------------------
  if (has_ggplot2) {
    # Missing months by group
    p_miss <- ggplot2::ggplot(miss_months, ggplot2::aes(x = interaction(GEAR, q, drop = TRUE), y = n_missing)) +
      ggplot2::geom_col() +
      ggplot2::labs(x = "GEAR × q", y = "Missing months (count)", title = paste0(label, ": missing months by GEAR×q")) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_plot_pdf_tiff(p_miss, paste0("c1_inputcheck__", label, "__missing_months_by_group"))
    
    # NA rate by group
    p_na <- ggplot2::ggplot(na_by_group, ggplot2::aes(x = interaction(GEAR, q, drop = TRUE), y = na_rate_value)) +
      ggplot2::geom_col() +
      ggplot2::labs(x = "GEAR × q", y = "NA rate in value", title = paste0(label, ": NA rate by GEAR×q")) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_plot_pdf_tiff(p_na, paste0("c1_inputcheck__", label, "__na_rate_by_group"))
  }
  
  # -----------------------------
  # Console summary
  # -----------------------------
  cat("\nSummary: ", label, "\n", sep = "")
  print(basic)
  cat("\nMonth range OK (1..12): ", month_ok, "\n", sep = "")
  cat("Duplicate cells (GEAR,q,YEAR,MONTH): ", nrow(dup_cells), "\n", sep = "")
  print(dup_summary)
  
  invisible(list(
    ts_long = ts_long,
    basic = basic,
    month_ok = month_ok,
    dup_cells = dup_cells,
    dup_summary = dup_summary,
    miss_months = miss_months,
    na_by_group = na_by_group,
    out_long_rds = out_long_rds
  ))
}

# =============================================================================
# Run checks for all c1 inputs
# =============================================================================
res_gamma <- check_one_input(IN_GAMMA, "gamma_q012")

# Alpha check is optional: if missing, we skip (gamma is preferred input for c1)
if (file.exists(IN_ALPHA)) {
  res_alpha <- check_one_input(IN_ALPHA, "alpha_q012")
} else {
  cat("\n(Info) alpha_q012.rds not found; skipped.\n")
}

cat(
  "\nDone: c1_inputs_check.R\n",
  "Derived (coerced) long RDS in:\n  - ", dir_data_coerced, "\n",
  "Diagnostics in:\n  - ", dir_diag_script, "\n",
  sep = ""
)