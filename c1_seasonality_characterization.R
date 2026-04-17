# =============================================================================
# Script: c1_seasonality_characterization.R
#
# Purpose:
#   Quantify seasonality strength for each gear × q (q = 0, 1, 2) using
#   harmonic regression on a monthly time series, and generate diagnostics
#   and figures for seasonal cycles.
#
# Inputs (authoritative paths per project conventions):
#   - 01_Analysis/Data/01_TimeSeries/gamma_q012.RDS (preferred): coastwide monthly γ per gear × q
#   - 01_Analysis/Data/b6_compute_hill_alpha_q012/alpha_q012.rds (fallback): bay-level monthly α; aggregated to coastwide
#
# Outputs (per project conventions):
#   Derived data:
#     - 01_Analysis/Data/c1_seasonality_characterization/seasonality_summary.csv
#     - 01_Analysis/Data/c1_seasonality_characterization/seasonality_fits.RDS
#     - 01_Analysis/Data/c1_seasonality_characterization/seasonal_curves_long.rds
#     - 01_Analysis/Data/c1_seasonality_characterization/seasonality_summary__text_table.txt
#   Diagnostics + figures:
#     - 02_Results/03_Diagnostics/c1_seasonality_characterization/*.csv
#     - 02_Results/01_Figures/c1_seasonality_characterization/*.pdf and *.tiff
#
# Notes:
#   - q = 2 is dominance-weighted diversity in labels.
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
script_stem <- tools::file_path_sans_ext(basename(sys.frames()[[1]]$ofile %||% "c1_seasonality_characterization.R"))
if (is.null(script_stem) || !nzchar(script_stem) || script_stem == "NULL") script_stem <- "c1_seasonality_characterization"

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
INFILE_GAMMA <- file.path(PATHS$analysis_dir, "Data", "01_TimeSeries", "gamma_q012.RDS")
INFILE_ALPHA <- file.path(PATHS$analysis_dir, "Data", "b6_compute_hill_alpha_q012", "alpha_q012.rds")

# -----------------------------
# Options (unchanged logic)
# -----------------------------
USE_SECOND_HARMONIC <- TRUE
DO_SPLIT_HALF       <- TRUE
MIN_N_MONTHS        <- 60
LJUNGBOX_LAG        <- 12
ACF_MAX_LAG         <- 24

# -----------------------------
# Labels + harmonics
# -----------------------------
q_label <- function(q) {
  if (is.na(q)) return(NA_character_)
  if (q == 0) return("q = 0 (richness)")
  if (q == 1) return("q = 1 (Shannon-weighted diversity)")
  if (q == 2) return("q = 2 (dominance-weighted diversity)")
  paste0("q = ", q)
}

make_harmonics <- function(df, use_second = TRUE) {
  df %>%
    mutate(
      m = as.integer(MONTH),
      sin1 = sin(2 * pi * m / 12),
      cos1 = cos(2 * pi * m / 12),
      sin2 = if (use_second) sin(4 * pi * m / 12) else 0,
      cos2 = if (use_second) cos(4 * pi * m / 12) else 0
    )
}

# -----------------------------
# Coerce input to: GEAR, q, YEAR, MONTH, value
# -----------------------------
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
  
  # If already long
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
    
    # If bay exists (alpha), aggregate to coastwide within gear×q×time
    if ("BAY_ID" %in% names(out)) {
      out <- out %>%
        group_by(GEAR, q, YEAR, MONTH) %>%
        summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
    }
    return(out)
  }
  
  # Otherwise pivot wide q0/q1/q2 into long
  q_cols <- intersect(names(df), c("q0", "q1", "q2"))
  if (!all(c("GEAR", "YEAR", "MONTH") %in% names(df)) || length(q_cols) == 0) {
    stop(
      "Could not coerce input to long format.\n",
      "Source: ", source_name, "\nColumns:\n  ", paste(names(df), collapse = ", ")
    )
  }
  
  long <- df %>%
    mutate(
      GEAR  = as.character(GEAR),
      YEAR  = suppressWarnings(as.integer(YEAR)),
      MONTH = suppressWarnings(as.integer(MONTH))
    ) %>%
    pivot_longer(cols = all_of(q_cols), names_to = "q_raw", values_to = "value") %>%
    mutate(
      q = suppressWarnings(as.numeric(str_extract(q_raw, "[012]"))),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    filter(!is.na(q), MONTH >= 1, MONTH <= 12)
  
  if ("BAY_ID" %in% names(df)) {
    long <- long %>%
      group_by(GEAR, q, YEAR, MONTH) %>%
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  }
  
  long %>% select(GEAR, q, YEAR, MONTH, value)
}

# -----------------------------
# Harmonic fit
# -----------------------------
fit_harmonic <- function(df, use_second = TRUE) {
  
  df2 <- df %>%
    filter(!is.na(value), !is.na(MONTH)) %>%
    make_harmonics(use_second = use_second) %>%
    arrange(YEAR, MONTH)  
  
  n_obs <- nrow(df2)
  n_months <- n_distinct(paste(df2$YEAR, df2$MONTH, sep = "-"))
  
  if (n_months < MIN_N_MONTHS) {
    return(list(ok = FALSE,
                reason = paste0("Insufficient months: ", n_months, " < ", MIN_N_MONTHS),
                n_obs = n_obs, n_months = n_months))
  }
  
  fml <- if (use_second) value ~ sin1 + cos1 + sin2 + cos2 else value ~ sin1 + cos1
  mod <- stats::lm(fml, data = df2)
  
  co <- stats::coef(mod)
  A <- unname(co["sin1"]); B <- unname(co["cos1"])
  amp_annual <- sqrt(A^2 + B^2)
  peak_to_trough <- 2 * amp_annual
  
  resid <- stats::residuals(mod)
  var_y <- stats::var(df2$value)
  var_r <- stats::var(resid)
  var_expl <- if (is.na(var_y) || var_y == 0) NA_real_ else (1 - (var_r / var_y))
  
  ac <- stats::acf(resid, plot = FALSE, lag.max = ACF_MAX_LAG)
  ac_vals <- as.numeric(ac$acf)[-1]
  lag1 <- if (length(ac_vals) >= 1) ac_vals[1] else NA_real_
  max_abs <- if (length(ac_vals) >= 1) max(abs(ac_vals), na.rm = TRUE) else NA_real_
  
  lb_p <- tryCatch({
    stats::Box.test(resid, lag = LJUNGBOX_LAG, type = "Ljung-Box")$p.value
  }, error = function(e) NA_real_)
  
  months_df <- tibble(MONTH = 1:12) %>%
    make_harmonics(use_second = use_second) %>%
    mutate(value_hat = stats::predict(mod, newdata = .)) %>%
    mutate(seasonal_centered = value_hat - mean(value_hat, na.rm = TRUE))
  
  list(
    ok = TRUE,
    n_obs = n_obs,
    n_months = n_months,
    model = mod,
    metrics = tibble(
      seasonal_amplitude = amp_annual,
      seasonal_peak_to_trough = peak_to_trough,
      variance_explained_season = var_expl,
      resid_acf_lag1 = lag1,
      resid_acf_max_abs_lag1_24 = max_abs,
      ljungbox_p_lag12 = lb_p
    ),
    seasonal_curve = months_df
  )
}

fit_split_half <- function(df, use_second = TRUE) {
  yrs <- sort(unique(df$YEAR))
  if (length(yrs) < 4) return(tibble(
    split_half_ok = FALSE, split_year = NA_integer_,
    amp_early = NA_real_, amp_late = NA_real_, amp_diff_late_minus_early = NA_real_
  ))
  
  split_year <- yrs[ceiling(length(yrs) / 2)]
  df_early <- df %>% filter(YEAR <= split_year)
  df_late  <- df %>% filter(YEAR >  split_year)
  
  fe <- fit_harmonic(df_early, use_second = use_second)
  fl <- fit_harmonic(df_late,  use_second = use_second)
  
  tibble(
    split_half_ok = isTRUE(fe$ok) & isTRUE(fl$ok),
    split_year = split_year,
    amp_early = if (isTRUE(fe$ok)) fe$metrics$seasonal_amplitude else NA_real_,
    amp_late  = if (isTRUE(fl$ok)) fl$metrics$seasonal_amplitude else NA_real_,
    amp_diff_late_minus_early =
      (if (isTRUE(fl$ok)) fl$metrics$seasonal_amplitude else NA_real_) -
      (if (isTRUE(fe$ok)) fe$metrics$seasonal_amplitude else NA_real_)
  )
}

# -----------------------------
# Load input (gamma preferred)
# -----------------------------
infile <- if (file.exists(INFILE_GAMMA)) INFILE_GAMMA else INFILE_ALPHA
if (!file.exists(infile)) {
  stop("No input found.\nTried:\n  - ", INFILE_GAMMA, "\n  - ", INFILE_ALPHA)
}
cat("Using input: ", infile, "\n", sep = "")

raw_obj <- readRDS(infile)
ts_dat <- coerce_to_long_ts(raw_obj, source_name = basename(infile))

stopifnot(all(c("GEAR", "q", "YEAR", "MONTH", "value") %in% names(ts_dat)))
ts_dat <- ts_dat %>% filter(MONTH >= 1, MONTH <= 12)

# Diagnostics: duplicates in ts_dat (should be none after coastwide aggregation)
dup_cells <- ts_dat %>%
  count(GEAR, q, YEAR, MONTH, name = "n") %>%
  filter(n > 1) %>%
  arrange(desc(n), GEAR, q, YEAR, MONTH)

out_dup <- file.path(dir_diag_script, "c1_seasonality__dup_cells_after_coercion.csv")
readr::write_csv(dup_cells, out_dup)

# -----------------------------
# Fit seasonality per gear × q
# -----------------------------
gq_keys <- ts_dat %>% distinct(GEAR, q) %>% arrange(GEAR, q)

fit_one_gq <- function(gear, qq) {
  df_gq <- ts_dat %>% filter(GEAR == gear, q == qq)
  
  fit <- fit_harmonic(df_gq, use_second = USE_SECOND_HARMONIC)
  
  split_metrics <- if (DO_SPLIT_HALF) {
    fit_split_half(df_gq, use_second = USE_SECOND_HARMONIC)
  } else {
    tibble(split_half_ok = NA, split_year = NA_integer_,
           amp_early = NA_real_, amp_late = NA_real_, amp_diff_late_minus_early = NA_real_)
  }
  
  if (!isTRUE(fit$ok)) {
    tibble(GEAR = gear, q = qq, q_label = q_label(qq),
           ok = FALSE, reason = fit$reason, n_obs = fit$n_obs, n_months = fit$n_months) %>%
      bind_cols(
        tibble(seasonal_amplitude = NA_real_, seasonal_peak_to_trough = NA_real_,
               variance_explained_season = NA_real_, resid_acf_lag1 = NA_real_,
               resid_acf_max_abs_lag1_24 = NA_real_, ljungbox_p_lag12 = NA_real_),
        split_metrics
      )
  } else {
    tibble(GEAR = gear, q = qq, q_label = q_label(qq),
           ok = TRUE, reason = NA_character_, n_obs = fit$n_obs, n_months = fit$n_months) %>%
      bind_cols(fit$metrics, split_metrics)
  }
}

fits_tbl <- purrr::pmap_dfr(
  list(gq_keys$GEAR, gq_keys$q),
  fit_one_gq
) %>%
  arrange(GEAR, q)

# Save summary CSV (derived data)
outfile_csv <- file.path(dir_data_script, "seasonality_summary.csv")
readr::write_csv(fits_tbl, outfile_csv)

# Save also to diagnostics bucket (handy for quick browsing)
outfile_csv_diag <- file.path(dir_diag_script, "seasonality_summary.csv")
readr::write_csv(fits_tbl, outfile_csv_diag)

# Save fits RDS (derived data)
gq_list <- ts_dat %>%
  group_by(GEAR, q) %>%
  group_split()

gq_names <- vapply(
  gq_list,
  function(x) paste0(unique(x$GEAR), "__q", unique(x$q)),
  character(1)
)

gq_list <- stats::setNames(gq_list, gq_names)

fits_full <- lapply(
  gq_list,
  function(x) fit_harmonic(x, use_second = USE_SECOND_HARMONIC)
)

outfile_rds <- file.path(dir_data_script, "seasonality_fits.RDS")
saveRDS(fits_full, outfile_rds)

# -----------------------------
# Figures: seasonal cycle curves (PDF + TIFF) (optional if ggplot2 installed)
# -----------------------------
seasonal_curves <- purrr::imap_dfr(fits_full, function(fit, nm) {
  if (!isTRUE(fit$ok)) return(NULL)
  parts <- stringr::str_split(nm, "__q", simplify = TRUE)
  gear <- parts[1]
  qq   <- as.numeric(parts[2])
  
  fit$seasonal_curve %>%
    transmute(GEAR = gear, q = qq, q_label = q_label(qq),
              MONTH, seasonal_centered)
})

# Save seasonal curves as derived data
saveRDS(seasonal_curves, file.path(dir_data_script, "seasonal_curves_long.rds"))

if (has_ggplot2 && nrow(seasonal_curves) > 0) {
  p <- ggplot2::ggplot(seasonal_curves, ggplot2::aes(x = MONTH, y = seasonal_centered, group = 1)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(breaks = 1:12) +
    ggplot2::labs(
      x = "Month",
      y = "Seasonal component (centered)",
      title = "Seasonal cycle (harmonic regression; centered)",
      subtitle = if (USE_SECOND_HARMONIC) "Annual + semiannual harmonics" else "Annual harmonic only"
    ) +
    ggplot2::facet_grid(q_label ~ GEAR, scales = "free_y") +
    ggplot2::theme_bw()
  
  fig_pdf  <- file.path(dir_fig_script, "seasonal_cycle_curves__by_q_by_gear.pdf")
  fig_tiff <- file.path(dir_fig_script, "seasonal_cycle_curves__by_q_by_gear.tiff")
  
  ggplot2::ggsave(fig_pdf,  plot = p, width = 12, height = 8, units = "in", dpi = 300)
  ggplot2::ggsave(fig_tiff, plot = p, width = 12, height = 8, units = "in", dpi = 300, compression = "lzw")
}

# -----------------------------
# Optional: clean text table (copy/paste into Word) (derived data)
# -----------------------------
txt_table <- fits_tbl %>%
  mutate(
    ljungbox_sig = case_when(
      is.na(ljungbox_p_lag12) ~ NA_character_,
      ljungbox_p_lag12 < 0.01 ~ "<0.01",
      ljungbox_p_lag12 < 0.05 ~ "<0.05",
      TRUE ~ "n.s."
    )
  ) %>%
  select(
    GEAR, q_label,
    seasonal_amplitude,
    seasonal_peak_to_trough,
    variance_explained_season,
    resid_acf_lag1,
    resid_acf_max_abs_lag1_24,
    ljungbox_p_lag12,
    ljungbox_sig,
    split_year, amp_early, amp_late, amp_diff_late_minus_early,
    ok, reason
  )

outfile_txt <- file.path(dir_data_script, "seasonality_summary__text_table.txt")
sink(outfile_txt)
cat("Seasonality summary (harmonic regression)\n")
cat("q = 2 is dominance-weighted diversity\n")
cat("Input used: ", infile, "\n\n", sep = "")
print(txt_table, n = nrow(txt_table), width = Inf)
sink()

message(
  "Done: c1_seasonality_characterization.R\n",
  "Input used:\n  - ", infile, "\n",
  "Derived data:\n",
  "  - ", outfile_csv, "\n",
  "  - ", outfile_rds, "\n",
  "  - ", outfile_txt, "\n",
  "Diagnostics:\n",
  "  - ", dir_diag_script, "\n",
  "Figures:\n",
  "  - ", dir_fig_script, "\n"
)