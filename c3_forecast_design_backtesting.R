# =============================================================================
# Script: c3_forecast_design_backtesting.R
#
# Purpose:
#   Define rolling-origin backtesting folds for forecasting (per series × gear × q),
#   using deseasonalized series from c2.
#
# Inputs (authoritative; from c2 outputs under 01_Analysis/Data/c2_deseasonalize_alpha_gamma/):
#   - alpha_q012_deseas.rds
#   - gamma_q012_deseas.rds
#
# Outputs (per project conventions):
#   Derived:
#   - 01_Analysis/Data/c3_forecast_design_backtesting/backtest_folds.rds
#
# Diagnostics:
#   - 02_Results/03_Diagnostics/c3_forecast_design_backtesting/backtest_groups_check.csv
#   - 02_Results/01_Figures/c3_forecast_design_backtesting/folds_timeaxis__h*.pdf + *.tiff (if ggplot2 installed)
#
# Notes:
#   - If alpha series includes BAY_ID, we aggregate to coastwide within gear×q×month
#     for this forecast design (one series per gear×q).
#   - Folds are built on a COMPLETE monthly grid so “12 months ahead” is true months.
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
need_pkgs <- c("dplyr", "tidyr", "readr")
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# Script stem + directories
# -----------------------------
script_stem <- tools::file_path_sans_ext(basename(sys.frames()[[1]]$ofile %||% "c3_forecast_design_backtesting.R"))
if (is.null(script_stem) || !nzchar(script_stem) || script_stem == "NULL") script_stem <- "c3_forecast_design_backtesting"

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
IN_ALPHA <- file.path(PATHS$analysis_dir, "Data", "c2_deseasonalize_alpha_gamma", "alpha_q012_deseas.rds")
IN_GAMMA <- file.path(PATHS$analysis_dir, "Data", "c2_deseasonalize_alpha_gamma", "gamma_q012_deseas.rds")

# Output (derived)
OUT_FOLDS <- file.path(dir_data_script, "backtest_folds.rds")

# -----------------------------
# Forecast design parameters (unchanged logic)
# -----------------------------
HORIZONS <- c(1, 12)
TRAINING_RULE <- "expanding"    # "expanding" or "fixed"
MIN_TRAIN_MONTHS <- 120         # counted as NON-NA observations
FIXED_TRAIN_MONTHS <- 120       # used only if TRAINING_RULE == "fixed"
ORIGIN_STEP <- 1
HOLDOUT_LAST_N <- 0
TARGET_COL <- "value_deseas"

# -----------------------------
# Helpers
# -----------------------------
expect_cols <- function(df, required, name = "data") {
  miss <- setdiff(required, names(df))
  if (length(miss) > 0) stop(name, " missing required columns: ", paste(miss, collapse = ", "))
}

make_group_id <- function(series_name, gear, q) {
  paste(series_name, gear, paste0("q", q), sep = "__")
}

complete_monthly <- function(df) {
  df <- df %>% arrange(date)
  d0 <- min(df$date, na.rm = TRUE)
  d1 <- max(df$date, na.rm = TRUE)
  grid <- tibble(date = seq.Date(d0, d1, by = "month"))
  grid %>% left_join(df, by = "date") %>% arrange(date)
}

build_folds_one_series <- function(df, series_name, horizons,
                                   training_rule = "expanding",
                                   min_train_months = 120,
                                   fixed_train_months = 120,
                                   origin_step = 1,
                                   holdout_last_n = 0) {
  
  df <- df %>% arrange(date)
  n <- nrow(df)
  
  n_effective <- if (holdout_last_n > 0) max(0, n - holdout_last_n) else n
  
  if (n_effective < (min_train_months + max(horizons) + 1)) {
    return(list(ok = FALSE,
                reason = "Too short for requested min_train + max_horizon",
                folds = list(), n = n, n_effective = n_effective))
  }
  
  folds_out <- list()
  
  for (h in horizons) {
    origins <- seq(from = 1, to = (n_effective - h), by = origin_step)
    
    folds_h <- list()
    k <- 0L
    
    for (t0 in origins) {
      
      train_idx <- if (training_rule == "expanding") {
        1:t0
      } else {
        start <- max(1, t0 - fixed_train_months + 1)
        start:t0
      }
      
      test_idx <- t0 + h
      
      # Require test target present
      if (is.na(df$y[test_idx])) next
      
      # Require enough NON-NA training observations
      n_train_ok <- sum(!is.na(df$y[train_idx]))
      if (n_train_ok < min_train_months) next
      
      k <- k + 1L
      folds_h[[k]] <- list(
        series = series_name,
        horizon = h,
        origin_i = t0,
        origin_date = df$date[t0],
        test_i = test_idx,
        test_date = df$date[test_idx],
        train_idx = train_idx,
        test_idx = test_idx,
        n_train = length(train_idx),
        n_train_non_na = n_train_ok,
        n_test = 1
      )
    }
    
    folds_out[[paste0("h", h)]] <- folds_h
  }
  
  list(ok = TRUE, reason = NA_character_, folds = folds_out, n = n, n_effective = n_effective)
}

folds_to_tibble <- function(groups_folds) {
  rows <- lapply(seq_along(groups_folds), function(i) {
    g <- groups_folds[[i]]
    gid <- g$group_id
    if (!isTRUE(g$ok)) return(NULL)
    
    hnames <- names(g$folds)
    out_h <- lapply(hnames, function(hname) {
      folds_h <- g$folds[[hname]]
      if (length(folds_h) == 0) return(NULL)
      
      do.call(rbind, lapply(folds_h, function(f) {
        data.frame(
          group_id = gid,
          series = g$series,
          GEAR = g$GEAR,
          q = g$q,
          horizon = f$horizon,
          origin_date = as.Date(f$origin_date),
          test_date = as.Date(f$test_date),
          train_start = g$dates[min(f$train_idx)],
          train_end = g$dates[max(f$train_idx)],
          n_train = f$n_train,
          n_train_non_na = f$n_train_non_na,
          stringsAsFactors = FALSE
        )
      }))
    })
    dplyr::bind_rows(out_h)
  })
  dplyr::bind_rows(rows)
}

# -----------------------------
# Load time series
# -----------------------------
stopifnot(file.exists(IN_ALPHA))
stopifnot(file.exists(IN_GAMMA))

alpha <- readRDS(IN_ALPHA)
gamma <- readRDS(IN_GAMMA)

required <- c("GEAR", "q", "date", TARGET_COL)
expect_cols(alpha, required, "alpha")
expect_cols(gamma, required, "gamma")

# Alpha: if BAY_ID exists, aggregate to coastwide per gear×q×date
alpha2 <- alpha %>%
  transmute(
    series = "alpha",
    BAY_ID = if ("BAY_ID" %in% names(alpha)) as.integer(.data$BAY_ID) else NA_integer_,
    GEAR = as.character(GEAR),
    q = as.numeric(q),
    date = as.Date(date),
    y = as.numeric(.data[[TARGET_COL]])
  )

if ("BAY_ID" %in% names(alpha)) {
  alpha2 <- alpha2 %>%
    group_by(series, GEAR, q, date) %>%
    summarise(y = mean(y, na.rm = TRUE), .groups = "drop")
} else {
  alpha2 <- alpha2 %>% select(series, GEAR, q, date, y)
}

gamma2 <- gamma %>%
  transmute(
    series = "gamma",
    GEAR = as.character(GEAR),
    q = as.numeric(q),
    date = as.Date(date),
    y = as.numeric(.data[[TARGET_COL]])
  )

ts_all <- bind_rows(alpha2, gamma2) %>%
  filter(!is.na(date)) %>%
  arrange(series, GEAR, q, date)

# -----------------------------
# Build folds per (series × gear × q) on complete monthly grids
# -----------------------------
groups_split <- ts_all %>%
  group_by(series, GEAR, q) %>%
  group_split()

group_ids <- vapply(
  groups_split,
  function(x) make_group_id(unique(x$series), unique(x$GEAR), unique(x$q)),
  character(1)
)

groups_named <- stats::setNames(groups_split, group_ids)

groups_folds <- lapply(names(groups_named), function(gid) {
  
  df_g <- groups_named[[gid]]
  
  series_name <- unique(df_g$series)
  gear <- unique(df_g$GEAR)
  qq <- unique(df_g$q)
  
  df_g2 <- df_g %>%
    select(date, y) %>%
    complete_monthly()
  
  res <- build_folds_one_series(
    df = df_g2,
    series_name = series_name,
    horizons = HORIZONS,
    training_rule = TRAINING_RULE,
    min_train_months = MIN_TRAIN_MONTHS,
    fixed_train_months = FIXED_TRAIN_MONTHS,
    origin_step = ORIGIN_STEP,
    holdout_last_n = HOLDOUT_LAST_N
  )
  
  list(
    ok = res$ok,
    reason = res$reason,
    series = series_name,
    GEAR = gear,
    q = qq,
    group_id = gid,
    dates = df_g2$date,
    y = df_g2$y,
    n = res$n,
    n_effective = res$n_effective,
    folds = res$folds
  )
})

# attach names so downstream can index by group_id
names(groups_folds) <- names(groups_named)

backtest_folds <- list(
  created = Sys.time(),
  target_col = TARGET_COL,
  horizons = HORIZONS,
  training_rule = TRAINING_RULE,
  min_train_months = MIN_TRAIN_MONTHS,
  fixed_train_months = if (TRAINING_RULE == "fixed") FIXED_TRAIN_MONTHS else NA_integer_,
  origin_step = ORIGIN_STEP,
  holdout_last_n = HOLDOUT_LAST_N,
  groups = groups_folds
)

saveRDS(backtest_folds, OUT_FOLDS)

# -----------------------------
# Diagnostics: group checks + fold plot
# -----------------------------
check_tbl <- dplyr::bind_rows(lapply(names(backtest_folds$groups), function(gid) {
  g <- backtest_folds$groups[[gid]]
  
  if (!isTRUE(g$ok)) {
    return(tibble(
      group_id = gid, series = g$series, GEAR = g$GEAR, q = g$q,
      ok = FALSE, reason = g$reason, n = g$n, n_effective = g$n_effective,
      n_folds = 0L
    ))
  }
  
  n_folds <- sum(vapply(g$folds, length, integer(1)))
  
  min_train_non_na <- suppressWarnings({
    mins <- unlist(lapply(g$folds, function(folds_h) {
      if (length(folds_h) == 0) return(numeric(0))
      vapply(folds_h, function(f) f$n_train_non_na, numeric(1))
    }))
    if (length(mins) == 0) NA_real_ else min(mins, na.rm = TRUE)
  })
  
  tibble(
    group_id = gid, series = g$series, GEAR = g$GEAR, q = g$q,
    ok = TRUE, reason = NA_character_, n = g$n, n_effective = g$n_effective,
    n_folds = n_folds, min_n_train_non_na = min_train_non_na
  )
}))

readr::write_csv(check_tbl, file.path(dir_diag_script, "backtest_groups_check.csv"))

folds_diag <- folds_to_tibble(backtest_folds$groups)

if (has_ggplot2 && nrow(folds_diag) > 0) {
  
  folds_diag <- folds_diag %>%
    mutate(
      q_label = case_when(
        q == 0 ~ "q = 0",
        q == 1 ~ "q = 1",
        q == 2 ~ "q = 2",
        TRUE ~ paste0("q = ", q)
      )
    )
  
  for (h in HORIZONS) {
    d <- folds_diag %>% filter(horizon == h)
    if (nrow(d) == 0) next
    
    p <- ggplot2::ggplot(d, ggplot2::aes(y = interaction(series, GEAR, q_label, sep = " | "))) +
      ggplot2::geom_segment(ggplot2::aes(x = train_start, xend = train_end,
                                         yend = interaction(series, GEAR, q_label, sep = " | "))) +
      ggplot2::geom_point(ggplot2::aes(x = test_date), size = 0.8) +
      ggplot2::labs(
        x = "Date",
        y = "Series | Gear | q",
        title = paste0("Rolling-origin folds (horizon = ", h, " month)"),
        subtitle = paste0(
          "Training: ", TRAINING_RULE,
          " | Min train (non-NA months): ", MIN_TRAIN_MONTHS,
          " | Step: ", ORIGIN_STEP
        )
      ) +
      ggplot2::theme_bw()
    
    ht <- max(6, 0.22 * dplyr::n_distinct(d$group_id))
    
    ggplot2::ggsave(file.path(dir_fig_script, paste0("folds_timeaxis__h", h, ".pdf")),
                    p, width = 14, height = ht, units = "in", dpi = 300)
    
    ggplot2::ggsave(file.path(dir_fig_script, paste0("folds_timeaxis__h", h, ".tiff")),
                    p, width = 14, height = ht, units = "in", dpi = 300, compression = "lzw")
  }
}

message(
  "Done: c3_forecast_design_backtesting.R\n",
  "Saved folds:\n  - ", OUT_FOLDS, "\n",
  "Diagnostics:\n  - ", file.path(dir_diag_script, "backtest_groups_check.csv"), "\n",
  "Figures:\n  - ", dir_fig_script, "\n"
)