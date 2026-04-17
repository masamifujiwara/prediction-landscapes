# =============================================================================
# Script: c2_deseasonalize_alpha_gamma.R
#
# Purpose:
#   Produce deseasonalized monthly time series for α and γ diversity (q = 0,1,2),
#   by gear, using the harmonic seasonality fits from c1.
#
# Inputs (authoritative paths per project conventions):
#   - 01_Analysis/Data/b6_compute_hill_alpha_q012/alpha_q012.rds
#   - 01_Analysis/Data/01_TimeSeries/gamma_q012.RDS
#   - 01_Analysis/Data/c1_seasonality_characterization/seasonality_fits.RDS  (from c1)
#
# Outputs (per project conventions):
#   Derived data (RDS):
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma/alpha_q012_deseas.rds
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma/gamma_q012_deseas.rds
#   - 01_Analysis/Data/c2_deseasonalize_alpha_gamma/seasonal_map.rds
#
# Diagnostics:
#   - 02_Results/03_Diagnostics/c2_deseasonalize_alpha_gamma/*.csv
#   - 02_Results/01_Figures/c2_deseasonalize_alpha_gamma/*.(pdf|tiff)
#
# Notes:
#   - q = 2 is dominance-weighted diversity in labels.
#   - USE_UNCENTERED_SEASONAL=TRUE removes the full fitted month mean (incl. intercept).
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
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# Script stem + directories
# -----------------------------
script_stem <- tools::file_path_sans_ext(basename(sys.frames()[[1]]$ofile %||% "c2_deseasonalize_alpha_gamma.R"))
if (is.null(script_stem) || !nzchar(script_stem) || script_stem == "NULL") script_stem <- "c2_deseasonalize_alpha_gamma"

dir_data_script <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_logs_script <- file.path(PATHS$analysis_dir, "Logs", script_stem)
dir.create(dir_data_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_logs_script, recursive = TRUE, showWarnings = FALSE)

dir_diag_script <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)
dir_fig_script  <- file.path(PATHS$results_dir, "01_Figures", script_stem)
dir.create(dir_diag_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig_script, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Inputs (authoritative)
# -----------------------------
IN_ALPHA <- file.path(PATHS$analysis_dir, "Data", "b6_compute_hill_alpha_q012", "alpha_q012.rds")
IN_GAMMA <- file.path(PATHS$analysis_dir, "Data", "01_TimeSeries", "gamma_q012.RDS")
IN_SEAS_FITS <- file.path(PATHS$analysis_dir, "Data", "c1_seasonality_characterization", "seasonality_fits.RDS")

# Outputs (derived)
OUT_ALPHA <- file.path(dir_data_script, "alpha_q012_deseas.rds")
OUT_GAMMA <- file.path(dir_data_script, "gamma_q012_deseas.rds")
OUT_SEASMAP <- file.path(dir_data_script, "seasonal_map.rds")

# Plot settings
ACF_MAX_LAG <- 36

# Deseasonalization option:
# TRUE: remove uncentered fitted month mean (includes intercept)
# FALSE: remove centered oscillation only (mean 0)
USE_UNCENTERED_SEASONAL <- TRUE

# -----------------------------
# Utility
# -----------------------------
q_label <- function(q) {
  dplyr::case_when(
    is.na(q)      ~ NA_character_,
    q == 0        ~ "q = 0 (richness)",
    q == 1        ~ "q = 1 (Shannon-weighted diversity)",
    q == 2        ~ "q = 2 (dominance-weighted diversity)",
    TRUE          ~ paste0("q = ", q)
  )
}

expect_cols <- function(df, required, name = "data") {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(name, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

# Coerce alpha/gamma objects to long format:
#   GEAR, q, YEAR, MONTH, value (+ optional BAY_ID)
coerce_to_long_ts <- function(obj, source_name = "input") {
  
  df <- obj
  if (is.list(obj) && !is.data.frame(obj)) {
    idx <- which(vapply(obj, is.data.frame, logical(1)))
    if (length(idx) >= 1) df <- obj[[idx[1]]]
  }
  stopifnot(is.data.frame(df))
  
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
  
  # Already long
  if (all(c("GEAR", "q", "YEAR", "MONTH", "value") %in% names(df))) {
    out <- df %>%
      mutate(
        GEAR  = as.character(GEAR),
        q     = suppressWarnings(as.numeric(q)),
        YEAR  = suppressWarnings(as.integer(YEAR)),
        MONTH = suppressWarnings(as.integer(MONTH)),
        value = suppressWarnings(as.numeric(value))
      ) %>%
      filter(MONTH >= 1, MONTH <= 12)
    return(out)
  }
  
  # Wide -> long
  q_cols <- intersect(c("q0", "q1", "q2"), names(df))
  if (!all(c("GEAR", "YEAR", "MONTH") %in% names(df)) || length(q_cols) == 0) {
    stop(
      "Could not coerce input to long format.\nSource: ", source_name,
      "\nColumns:\n  ", paste(names(df), collapse = ", ")
    )
  }
  
  df %>%
    mutate(
      GEAR  = as.character(GEAR),
      YEAR  = suppressWarnings(as.integer(YEAR)),
      MONTH = suppressWarnings(as.integer(MONTH))
    ) %>%
    pivot_longer(cols = all_of(q_cols), names_to = "q_raw", values_to = "value") %>%
    mutate(
      q = suppressWarnings(as.numeric(stringr::str_extract(q_raw, "[012]"))),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    filter(!is.na(q), MONTH >= 1, MONTH <= 12)
}

# Build month->seasonal mapping from c1 fits
# Returns: GEAR, q, MONTH, seasonal
build_seasonal_map <- function(seas_fits, use_uncentered = TRUE) {
  
  if (is.null(seas_fits) || !is.list(seas_fits)) stop("seas_fits is not a list. Did IN_SEAS_FITS load correctly?")
  
  nm <- names(seas_fits)
  if (is.null(nm) || any(!nzchar(nm))) stop("seas_fits must be a named list like 'BS__q0'.")
  
  rows <- lapply(seq_along(seas_fits), function(i) {
    fit <- seas_fits[[i]]
    nmi <- nm[[i]]
    if (!isTRUE(fit$ok)) return(NULL)
    
    parts <- stringr::str_split(nmi, "__q", simplify = TRUE)
    gear  <- parts[1]
    qq    <- suppressWarnings(as.numeric(parts[2]))
    
    sc <- fit$seasonal_curve
    if (!is.data.frame(sc) || !all(c("MONTH", "value_hat", "seasonal_centered") %in% names(sc))) {
      stop("seasonal_curve missing expected columns for ", nmi)
    }
    
    data.frame(
      GEAR = as.character(gear),
      q = as.numeric(qq),
      MONTH = as.integer(sc$MONTH),
      seasonal = if (use_uncentered) as.numeric(sc$value_hat) else as.numeric(sc$seasonal_centered),
      stringsAsFactors = FALSE
    )
  })
  
  out <- dplyr::bind_rows(rows)
  
  # Diagnostics: check duplicates in seasonal_map
  dup <- out %>% count(GEAR, q, MONTH, name = "n") %>% filter(n > 1)
  if (nrow(dup) > 0) {
    readr::write_csv(dup, file.path(dir_diag_script, "c2_seasonal_map__dup_GEAR_q_MONTH.csv"))
    stop("seasonal_map has duplicates (GEAR,q,MONTH). See c2_seasonal_map__dup_GEAR_q_MONTH.csv")
  }
  
  out
}

# Deseasonalize: join seasonal by GEAR×q×MONTH
deseasonalize_ts <- function(ts_df, seasonal_map) {
  required <- c("GEAR", "q", "YEAR", "MONTH", "value")
  expect_cols(ts_df, required, "time series")
  
  ts_df2 <- ts_df %>%
    mutate(
      GEAR  = as.character(GEAR),
      q     = suppressWarnings(as.numeric(q)),
      YEAR  = suppressWarnings(as.integer(YEAR)),
      MONTH = suppressWarnings(as.integer(MONTH)),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    filter(MONTH >= 1, MONTH <= 12) %>%
    mutate(
      date = as.Date(sprintf("%04d-%02d-01", YEAR, MONTH)),
      q_label = q_label(q)
    ) %>%
    arrange(GEAR, q, YEAR, MONTH)
  
  out <- ts_df2 %>%
    left_join(seasonal_map, by = c("GEAR", "q", "MONTH")) %>%
    mutate(
      seasonal_component = seasonal,
      residual_component = ifelse(is.na(seasonal_component), NA_real_, value - seasonal_component),
      value_deseas = residual_component
    ) %>%
    select(any_of("BAY_ID"), GEAR, q, q_label, YEAR, MONTH, date,
           value, seasonal_component, residual_component, value_deseas)
  
  out
}

# ACF/PACF plots (PDF + TIFF) (only if ggplot2 available)
plot_acf_pacf_group <- function(df_group, out_prefix) {
  if (!has_ggplot2) return(invisible(NULL))
  
  x <- df_group$residual_component
  x <- x[!is.na(x)]
  if (length(x) < 20) return(invisible(NULL))
  
  acf_obj  <- stats::acf(x, plot = FALSE, lag.max = ACF_MAX_LAG)
  pacf_obj <- stats::pacf(x, plot = FALSE, lag.max = ACF_MAX_LAG)
  
  acf_df  <- tibble(lag = as.numeric(acf_obj$lag),  acf  = as.numeric(acf_obj$acf))
  pacf_df <- tibble(lag = as.numeric(pacf_obj$lag), pacf = as.numeric(pacf_obj$acf))
  
  p1 <- ggplot2::ggplot(acf_df, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0)) +
    ggplot2::labs(x = "Lag", y = "ACF", title = "Residual ACF") +
    ggplot2::theme_bw()
  
  p2 <- ggplot2::ggplot(pacf_df, ggplot2::aes(x = lag, y = pacf)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0)) +
    ggplot2::labs(x = "Lag", y = "PACF", title = "Residual PACF") +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(paste0(out_prefix, "__acf.pdf"),  p1, width = 6, height = 4, units = "in", dpi = 300)
  ggplot2::ggsave(paste0(out_prefix, "__acf.tiff"), p1, width = 6, height = 4, units = "in", dpi = 300, compression = "lzw")
  ggplot2::ggsave(paste0(out_prefix, "__pacf.pdf"),  p2, width = 6, height = 4, units = "in", dpi = 300)
  ggplot2::ggsave(paste0(out_prefix, "__pacf.tiff"), p2, width = 6, height = 4, units = "in", dpi = 300, compression = "lzw")
  
  invisible(NULL)
}

# Decomposition plot (PDF + TIFF) (only if ggplot2 available)
plot_decomposition <- function(df, out_file_base, title_prefix, use_uncentered) {
  if (!has_ggplot2) return(invisible(NULL))
  
  df_long <- df %>%
    select(any_of(c("BAY_ID")), GEAR, q_label, date,
           value, seasonal_component, residual_component) %>%
    pivot_longer(
      cols = c(value, seasonal_component, residual_component),
      names_to = "component",
      values_to = "y"
    ) %>%
    mutate(component = recode(component,
                              value = "original",
                              seasonal_component = "seasonal",
                              residual_component = "residual"))
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = date, y = y)) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(q_label ~ GEAR, scales = "free_y") +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = paste0(title_prefix, ": original vs seasonal vs residual"),
      subtitle = if (use_uncentered) "Seasonal component = uncentered fitted month mean"
      else "Seasonal component = centered (mean 0) oscillation"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(paste0(out_file_base, ".pdf"),  p, width = 13, height = 9, units = "in", dpi = 300)
  ggplot2::ggsave(paste0(out_file_base, ".tiff"), p, width = 13, height = 9, units = "in", dpi = 300, compression = "lzw")
}

# =============================================================================
# Load inputs
# =============================================================================
stopifnot(file.exists(IN_ALPHA))
stopifnot(file.exists(IN_GAMMA))
stopifnot(file.exists(IN_SEAS_FITS))

alpha_ts  <- coerce_to_long_ts(readRDS(IN_ALPHA),  source_name = basename(IN_ALPHA))
gamma_ts  <- coerce_to_long_ts(readRDS(IN_GAMMA),  source_name = basename(IN_GAMMA))
seas_fits <- readRDS(IN_SEAS_FITS)

seasonal_map <- build_seasonal_map(seas_fits, use_uncentered = USE_UNCENTERED_SEASONAL)
saveRDS(seasonal_map, OUT_SEASMAP)

# Diagnostics: seasonal_map coverage (should be 12 months per gear×q)
diag_seasmap <- seasonal_map %>%
  count(GEAR, q, name = "n_months") %>%
  arrange(GEAR, q)
readr::write_csv(diag_seasmap, file.path(dir_diag_script, "c2_seasonal_map__n_months_by_group.csv"))

# Diagnostics: check duplicates in gamma_ts key
dup_gamma <- gamma_ts %>% count(GEAR, q, YEAR, MONTH, name = "n") %>% filter(n > 1)
readr::write_csv(dup_gamma, file.path(dir_diag_script, "c2_gamma__dup_cells.csv"))
if (nrow(dup_gamma) > 0) stop("gamma_ts has duplicated (GEAR,q,YEAR,MONTH). See c2_gamma__dup_cells.csv")

# =============================================================================
# Deseasonalize
# =============================================================================
alpha_deseas <- deseasonalize_ts(alpha_ts, seasonal_map)
gamma_deseas <- deseasonalize_ts(gamma_ts, seasonal_map)

saveRDS(alpha_deseas, OUT_ALPHA)
saveRDS(gamma_deseas, OUT_GAMMA)

cat("Saved:\n  - ", OUT_ALPHA, "\n  - ", OUT_GAMMA, "\n", sep = "")

# =============================================================================
# Diagnostics: decomposition plots + ACF/PACF
# =============================================================================
# Alpha: if BAY_ID present, sample a few bays for plots only
alpha_plot_df <- alpha_deseas
if ("BAY_ID" %in% names(alpha_deseas)) {
  set.seed(1)
  bays_tbl <- alpha_deseas %>% distinct(BAY_ID)
  n_bays <- nrow(bays_tbl)
  bays_keep <- bays_tbl %>% slice_sample(n = min(3, n_bays)) %>% pull(BAY_ID)
  alpha_plot_df <- alpha_deseas %>% filter(BAY_ID %in% bays_keep)
}

plot_decomposition(
  alpha_plot_df,
  out_file_base = file.path(dir_fig_script, "alpha_decomposition"),
  title_prefix = "Alpha diversity",
  use_uncentered = USE_UNCENTERED_SEASONAL
)
plot_decomposition(
  gamma_deseas,
  out_file_base = file.path(dir_fig_script, "gamma_decomposition"),
  title_prefix = "Gamma diversity",
  use_uncentered = USE_UNCENTERED_SEASONAL
)

run_acf_pacf_all <- function(df, prefix) {
  if (!has_ggplot2) return(invisible(NULL))
  
  keys <- df %>% distinct(GEAR, q, q_label) %>% arrange(GEAR, q)
  for (i in seq_len(nrow(keys))) {
    gear <- keys$GEAR[i]
    qq   <- keys$q[i]
    g <- df %>% filter(GEAR == gear, q == qq)
    
    out_prefix <- file.path(dir_fig_script, sprintf("%s__%s__q%s", prefix, gear, qq))
    plot_acf_pacf_group(g, out_prefix)
  }
}

run_acf_pacf_all(alpha_deseas, "alpha_resid")
run_acf_pacf_all(gamma_deseas, "gamma_resid")

message(
  "Done: c2_deseasonalize_alpha_gamma.R\n",
  "Derived outputs:\n",
  "  - ", OUT_ALPHA, "\n",
  "  - ", OUT_GAMMA, "\n",
  "  - ", OUT_SEASMAP, "\n",
  "Diagnostics:\n",
  "  - ", dir_diag_script, "\n",
  "Figures:\n",
  "  - ", dir_fig_script, "\n"
)