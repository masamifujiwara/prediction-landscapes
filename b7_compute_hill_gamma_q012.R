# =============================================================================
# Script: b7_compute_hill_gamma_q012.R
#
# Purpose:
#   Compute γ diversity (q = 0, 1, 2) per (YEAR × MONTH × GEAR) coastwide
#   (pooling across all bays within each gear), using COVERAGE-BASED
#   standardization consistent with b6 (α).
#
# Key requirements:
#   - Sentinel rows: SPECIES_CODE = NA, catch = 0 (stations with no organisms)
#     are EXCLUDED from abundance vectors.
#   - Time cells with sampling effort but zero organisms yield γ=0, Chat_obs=0.
#
# Directory convention (this project):
#   - Load PATHS from: 01_Analysis/PATHS.RData (but refresh to current machine)
#   - Inputs (derived data):
#       01_Analysis/Data/b3_harmonize_station_metadata/station_core.RDS
#       01_Analysis/Data/b4_join_organisms_to_station/bs_events.RDS
#       01_Analysis/Data/b4_join_organisms_to_station/bt_events.RDS
#       01_Analysis/Data/b5_define_spatiotemporal_units/index_keys.RDS
#   - Outputs (derived data):
#       01_Analysis/Data/01_TimeSeries/gamma_q012.RDS
#       01_Analysis/Data/01_TimeSeries/gamma_q012.csv
#   - Diagnostics (figures + csv):
#       02_Results/03_Diagnostics/b7_compute_hill_gamma_q012/...
#
# Notes (Linux stability):
#   - iNEXT::estimateD(datatype="abundance") assumes INTEGER abundances.
#     Non-integer pooled catch can crash during bootstrap (rmultinom).
#   - This script rounds pooled abundances after filtering finite/positive values
#     and can optionally stop early if non-integers are detected before rounding.
# =============================================================================

rm(list = ls())

set.seed(102)

# ─────────────────────────────────────────────────────────────────────────────
# 0) Script directory (RStudio / Rscript / source / nohup) + setwd()
# ─────────────────────────────────────────────────────────────────────────────
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
cat("Working directory:", normalizePath(getwd(), winslash = "/"), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1) Load PATHS (authoritative), then REFRESH to current machine/project_dir
# ─────────────────────────────────────────────────────────────────────────────
get_project_root <- function() {
  here <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  normalizePath(file.path(here, ".."), winslash = "/", mustWork = TRUE)
}
project_dir <- get_project_root()

paths_file <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file) # loads PATHS

`%||%` <- function(a, b) if (!is.null(a)) a else b

refresh_paths <- function(PATHS, project_dir) {
  # infer script name/stem robustly
  script_file <- PATHS$script_file %||% NA_character_
  if (is.na(script_file) || !nzchar(script_file) || script_file == "interactive_session.R") {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) == 1) script_file <- basename(sub("^--file=", "", file_arg))
    if (is.na(script_file) || !nzchar(script_file)) script_file <- "interactive_session.R"
  }
  script_stem <- sub("\\.[Rr]$", "", script_file)
  
  PATHS$project_dir  <- project_dir
  PATHS$analysis_dir <- file.path(project_dir, "01_Analysis")
  PATHS$logs_dir     <- file.path(PATHS$analysis_dir, "Logs")
  PATHS$data_dir     <- file.path(PATHS$analysis_dir, "Data")
  
  PATHS$results_dir  <- file.path(project_dir, "02_Results")
  PATHS$fig_root     <- file.path(PATHS$results_dir, "01_Figures")
  PATHS$tab_root     <- file.path(PATHS$results_dir, "02_Tables")
  PATHS$diag_root    <- file.path(PATHS$results_dir, "03_Diagnostics")
  PATHS$model_root   <- file.path(PATHS$results_dir, "04_Models")
  
  PATHS$script_file  <- script_file
  PATHS$script_stem  <- script_stem
  
  PATHS$fig_dir   <- file.path(PATHS$fig_root, script_stem)
  PATHS$tab_dir   <- file.path(PATHS$tab_root, script_stem)
  PATHS$diag_dir  <- file.path(PATHS$diag_root, script_stem)
  PATHS$model_dir <- file.path(PATHS$model_root, script_stem)
  
  # create fixed structure + allowed script subfolders
  dir.create(PATHS$logs_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(PATHS$data_dir,  recursive = TRUE, showWarnings = FALSE)
  
  dir.create(PATHS$fig_root,  recursive = TRUE, showWarnings = FALSE)
  dir.create(PATHS$tab_root,  recursive = TRUE, showWarnings = FALSE)
  dir.create(PATHS$diag_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(PATHS$model_root,recursive = TRUE, showWarnings = FALSE)
  
  dir.create(PATHS$fig_dir,   recursive = TRUE, showWarnings = FALSE)
  dir.create(PATHS$tab_dir,   recursive = TRUE, showWarnings = FALSE)
  dir.create(PATHS$diag_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(PATHS$model_dir, recursive = TRUE, showWarnings = FALSE)
  
  PATHS
}

PATHS <- refresh_paths(PATHS, project_dir)

cat("Project directory (refreshed):", PATHS$project_dir, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2) Options + logging helpers (nohup-safe; no emoji by default)
# ─────────────────────────────────────────────────────────────────────────────
options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

use_emoji <- FALSE
msg <- function(..., .nl = TRUE) {
  txt <- paste0(..., collapse = "")
  if (.nl) cat(txt, "\n") else cat(txt)
}
icon <- function(x) if (use_emoji) x else ""

# ─────────────────────────────────────────────────────────────────────────────
# 3) Packages
# ─────────────────────────────────────────────────────────────────────────────
required_pkgs <- c("dplyr", "tidyr", "purrr", "readr", "ggplot2", "iNEXT")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(iNEXT)
})

# ─────────────────────────────────────────────────────────────────────────────
# 4) Parameters
# ─────────────────────────────────────────────────────────────────────────────
nboot <- 100
conf  <- 0.95

C_cap      <- 0.90
C_quantile <- 0.10
C_min      <- 0.50  # diagnostics flag only

run_preflight_probe <- TRUE
nboot_probe <- 1

strict_noninteger_check <- FALSE

# ─────────────────────────────────────────────────────────────────────────────
# 5) Paths + load inputs (from 01_Analysis/Data/<script_stem>/...)
# ─────────────────────────────────────────────────────────────────────────────
in_station_core <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.RDS")
in_bs_events    <- file.path(PATHS$data_dir, "b4_join_organisms_to_station", "bs_events.RDS")
in_bt_events    <- file.path(PATHS$data_dir, "b4_join_organisms_to_station", "bt_events.RDS")
in_index_keys   <- file.path(PATHS$data_dir, "b5_define_spatiotemporal_units", "index_keys.RDS")

stopifnot(file.exists(in_station_core))
stopifnot(file.exists(in_bs_events))
stopifnot(file.exists(in_bt_events))
stopifnot(file.exists(in_index_keys))

station_core <- readRDS(in_station_core)
bs_events    <- readRDS(in_bs_events)
bt_events    <- readRDS(in_bt_events)
index_keys   <- readRDS(in_index_keys)

station_bs_core <- station_core$station_bs_core
station_bt_core <- station_core$station_bt_core

valid_bays <- index_keys$keys$valid_bays %||% 1:8

# Derived output folder for time series
ts_dir <- file.path(PATHS$data_dir, "01_TimeSeries")
dir.create(ts_dir, recursive = TRUE, showWarnings = FALSE)

# Diagnostics: MUST go under 02_Results/03_Diagnostics/<script_stem>/
diag_dir <- PATHS$diag_dir
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 6) Build γ “effort grid” (gear×year×month) from station_core (keeps effort-only cells)
# ─────────────────────────────────────────────────────────────────────────────
station_all <- bind_rows(station_bs_core, station_bt_core) %>%
  mutate(
    BAY_ID = as.integer(MAJOR_AREA),
    YEAR   = as.integer(YEAR_i),
    MONTH  = as.integer(MONTH_i)
  ) %>%
  filter(BAY_ID %in% valid_bays) %>%
  select(SAMPLE_ID, GEAR, YEAR, MONTH)

gamma_grid <- station_all %>%
  count(GEAR, YEAR, MONTH, name = "n_samples") %>%
  arrange(GEAR, YEAR, MONTH) %>%
  mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH))

# ─────────────────────────────────────────────────────────────────────────────
# 7) Utilities
# ─────────────────────────────────────────────────────────────────────────────
safe_chat <- function(x) {
  if (length(x) == 0) return(0)
  x <- as.numeric(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(0)
  
  N <- sum(x)
  if (!is.finite(N) || N <= 0) return(0)
  
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  
  if (f1 == 0) return(1)
  
  A <- if (f2 > 0) ((N - 1) * f1) / (((N - 1) * f1) + (2 * f2)) else 1
  C_hat <- 1 - (f1 / N) * A
  max(min(C_hat, 1), 0)
}

sanitize_abundance <- function(x, strict_noninteger_check = FALSE) {
  if (length(x) == 0) return(setNames(numeric(0), character(0)))
  
  nm <- names(x)
  x <- as.numeric(x)
  
  ok <- is.finite(x) & (x > 0)
  x  <- x[ok]
  if (!is.null(nm)) nm <- nm[ok]
  
  if (length(x) == 0) return(setNames(numeric(0), character(0)))
  
  if (strict_noninteger_check) {
    if (any(abs(x - round(x)) > 1e-8)) {
      stop("Non-integer abundances detected before rounding. ",
           "If 'catch' is not a count, datatype='abundance' may be inappropriate.",
           call. = FALSE)
    }
  }
  
  x <- round(x)
  keep <- x > 0
  x <- x[keep]
  if (!is.null(nm)) nm <- nm[keep]
  
  if (length(x) == 0) return(setNames(numeric(0), character(0)))
  
  if (sum(x) > .Machine$integer.max) {
    stop("sum(x) exceeds integer max; cannot bootstrap with rmultinom safely.", call. = FALSE)
  }
  
  if (!is.null(nm)) x <- setNames(x, nm)
  x
}

estimate_q012_at_coverage <- function(x, C_target, nboot, conf, strict_noninteger_check = FALSE) {
  x <- sanitize_abundance(x, strict_noninteger_check = strict_noninteger_check)
  
  if (length(x) == 0 || sum(x) <= 0) {
    return(tibble(
      q0 = 0, q1 = 0, q2 = 0,
      se_q0 = NA_real_, se_q1 = NA_real_, se_q2 = NA_real_,
      lcl_q0 = NA_real_, lcl_q1 = NA_real_, lcl_q2 = NA_real_,
      ucl_q0 = NA_real_, ucl_q1 = NA_real_, ucl_q2 = NA_real_
    ))
  }
  
  est <- iNEXT::estimateD(
    x = x,
    datatype = "abundance",
    base = "coverage",
    level = C_target,
    conf = conf,
    nboot = nboot
  )
  
  cn <- names(est)
  
  order_candidates <- c("Order.q", "q", "order", "Order", "q.d")
  order_col <- order_candidates[order_candidates %in% cn]
  if (length(order_col) == 0) {
    stop("estimateD() output has no order column. Columns:\n  - ",
         paste(cn, collapse = "\n  - "), call. = FALSE)
  }
  order_col <- order_col[1]
  
  est_candidates <- c("qD", "Estimate", "Diversity", "qD.Estimate", "Estimate.qD")
  est_col <- est_candidates[est_candidates %in% cn]
  if (length(est_col) == 0) {
    stop("estimateD() output has no estimate column. Columns:\n  - ",
         paste(cn, collapse = "\n  - "), call. = FALSE)
  }
  est_col <- est_col[1]
  
  se_col  <- intersect(c("s.e.", "se", "Std.Error", "StdError"), cn)
  lcl_col <- intersect(c("qD.LCL", "LCL", "LCI", "lcl", "Lower"), cn)
  ucl_col <- intersect(c("qD.UCL", "UCL", "UCI", "ucl", "Upper"), cn)
  
  se_col  <- if (length(se_col)  > 0) se_col[1]  else NA_character_
  lcl_col <- if (length(lcl_col) > 0) lcl_col[1] else NA_character_
  ucl_col <- if (length(ucl_col) > 0) ucl_col[1] else NA_character_
  
  ord <- suppressWarnings(as.numeric(est[[order_col]]))
  
  get_one <- function(q) {
    idx <- which(ord == q)
    if (length(idx) == 0) return(list(estimate = NA_real_, se = NA_real_, lcl = NA_real_, ucl = NA_real_))
    idx <- idx[1]
    list(
      estimate = as.numeric(est[[est_col]][idx]),
      se  = if (!is.na(se_col))  as.numeric(est[[se_col]][idx])  else NA_real_,
      lcl = if (!is.na(lcl_col)) as.numeric(est[[lcl_col]][idx]) else NA_real_,
      ucl = if (!is.na(ucl_col)) as.numeric(est[[ucl_col]][idx]) else NA_real_
    )
  }
  
  r0 <- get_one(0); r1 <- get_one(1); r2 <- get_one(2)
  
  tibble(
    q0 = r0$estimate, q1 = r1$estimate, q2 = r2$estimate,
    se_q0 = r0$se, se_q1 = r1$se, se_q2 = r2$se,
    lcl_q0 = r0$lcl, lcl_q1 = r1$lcl, lcl_q2 = r2$lcl,
    ucl_q0 = r0$ucl, ucl_q1 = r1$ucl, ucl_q2 = r2$ucl
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 8) Prepare events and pool abundances coastwide by (gear×year×month)
#     IMPORTANT: remove sentinel NA species rows before pooling
# ─────────────────────────────────────────────────────────────────────────────
prep_events <- function(events_df, gear_label, valid_bays) {
  events_df %>%
    mutate(
      GEAR = gear_label,
      BAY_ID = as.integer(MAJOR_AREA),
      YEAR  = as.integer(YEAR_i),
      MONTH = as.integer(MONTH_i),
      SPECIES_CODE = as.character(SPECIES_CODE),
      catch = as.numeric(catch)
    ) %>%
    filter(BAY_ID %in% valid_bays)
}

bs_events2 <- prep_events(bs_events, "BS", valid_bays)
bt_events2 <- prep_events(bt_events, "BT", valid_bays)
events_all <- bind_rows(bs_events2, bt_events2)

gamma_species <- events_all %>%
  filter(!is.na(SPECIES_CODE)) %>%  # sentinel removed
  group_by(GEAR, YEAR, MONTH, SPECIES_CODE) %>%
  summarise(abund = sum(catch, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(abund) & abund > 0)

gamma_vectors <- gamma_species %>%
  group_by(GEAR, YEAR, MONTH) %>%
  summarise(
    abund_vec = list(setNames(abund, SPECIES_CODE)),
    N_total   = sum(abund),
    S_obs     = n(),
    .groups   = "drop"
  )

gamma_cells <- gamma_grid %>%
  left_join(gamma_vectors, by = c("GEAR", "YEAR", "MONTH")) %>%
  mutate(
    N_total = coalesce(N_total, 0),
    S_obs   = coalesce(S_obs, 0)
  ) %>%
  mutate(
    abund_vec = map2(abund_vec, N_total, ~{
      if (is.null(.x) || is.na(.y) || .y == 0) setNames(numeric(0), character(0)) else .x
    })
  )

gamma_cells <- gamma_cells %>%
  mutate(Chat_obs = map_dbl(abund_vec, safe_chat))

chat_nonempty <- gamma_cells$Chat_obs[gamma_cells$N_total > 0]
if (length(chat_nonempty) == 0) stop("No non-empty gamma cells found. Check input events.", call. = FALSE)

C_target <- min(C_cap, as.numeric(stats::quantile(chat_nonempty, probs = C_quantile, na.rm = TRUE)))
C_target <- max(min(C_target, 0.999), 0.001)

msg(icon("🎯 "), "Coverage target (C_target): ", sprintf("%.3f", C_target),
    " (cap=", C_cap, ", quantile=", C_quantile, ")")

gamma_cells <- gamma_cells %>%
  mutate(
    extrapolated  = (N_total > 0) & (C_target > Chat_obs),
    coverage_fail = (N_total == 0) | (Chat_obs < C_min),
    C_target      = C_target
  )

# ─────────────────────────────────────────────────────────────────────────────
# 9) Diagnostics: non-integer check + optional preflight probe
# ─────────────────────────────────────────────────────────────────────────────
nonint_summary <- gamma_vectors %>%
  mutate(any_nonint = map_lgl(abund_vec, ~{
    v <- as.numeric(.x)
    v <- v[is.finite(v) & v > 0]
    if (length(v) == 0) return(FALSE)
    any(abs(v - round(v)) > 1e-8)
  })) %>%
  summarise(
    n_cells = n(),
    n_any_nonint = sum(any_nonint),
    prop_any_nonint = mean(any_nonint)
  )

msg("\n=== Abundance integrality check (before rounding) ===")
print(nonint_summary)

if (isTRUE(strict_noninteger_check) && nonint_summary$n_any_nonint > 0) {
  stop("Non-integer pooled abundances detected. Set strict_noninteger_check=FALSE to allow rounding, ",
       "or revisit what 'catch' represents in events.", call. = FALSE)
}

if (isTRUE(run_preflight_probe)) {
  msg("\n=== Preflight probe: estimateD per cell (nboot_probe=", nboot_probe, ") ===", sep = "")
  
  ok_vec <- logical(nrow(gamma_cells))
  
  for (i in seq_len(nrow(gamma_cells))) {
    
    if (gamma_cells$N_total[i] <= 0) {
      ok_vec[i] <- TRUE
      next
    }
    
    x <- sanitize_abundance(gamma_cells$abund_vec[[i]],
                            strict_noninteger_check = strict_noninteger_check)
    
    res <- try(iNEXT::estimateD(
      x = x,
      datatype = "abundance",
      base = "coverage",
      level = C_target,
      conf = conf,
      nboot = nboot_probe
    ), silent = TRUE)
    
    ok_vec[i] <- !inherits(res, "try-error")
    
    if (!ok_vec[i]) {
      msg("Preflight FAILED at row ", i,
          " (", gamma_cells$GEAR[i], " ",
          gamma_cells$YEAR[i], "-", sprintf("%02d", gamma_cells$MONTH[i]), "): ",
          as.character(res)[1])
    }
    
    if (i %% 200 == 0) msg("... preflight progress: ", i, "/", nrow(gamma_cells))
  }
  
  if (any(!ok_vec)) {
    bad <- gamma_cells[!ok_vec, ] %>%
      select(GEAR, YEAR, MONTH, YEAR_MONTH, n_samples, N_total, S_obs, Chat_obs, C_target)
    
    bad_path <- file.path(diag_dir, "b7_gamma_cells_that_fail_estimateD.csv")
    write_csv(bad, bad_path)
    msg("Saved failing-cell list to: ", bad_path)
    
    stop("Preflight probe found cells that fail estimateD() (see CSV above). ",
         "Fix those cells or adjust parameters before running full bootstrap.",
         call. = FALSE)
  } else {
    msg("Preflight probe passed for all non-empty cells.")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 10) Estimate γ q0/q1/q2 per cell
# ─────────────────────────────────────────────────────────────────────────────
msg("\n=== Estimating gamma diversity per cell (nboot=", nboot, ") ===", sep = "")

est_list <- vector("list", nrow(gamma_cells))
for (i in seq_len(nrow(gamma_cells))) {
  est_list[[i]] <- estimate_q012_at_coverage(
    x = gamma_cells$abund_vec[[i]],
    C_target = C_target,
    nboot = nboot,
    conf = conf,
    strict_noninteger_check = strict_noninteger_check
  )
  if (i %% 200 == 0) msg("... estimate progress: ", i, "/", nrow(gamma_cells))
}

gamma_est <- bind_cols(gamma_cells, bind_rows(est_list))

# ─────────────────────────────────────────────────────────────────────────────
# 11) Diagnostics outputs (ALL under 02_Results/03_Diagnostics/<script_stem>/)
# ─────────────────────────────────────────────────────────────────────────────
diag_cells_path <- file.path(diag_dir, "b7_gamma_coverage_diagnostics__cells.csv")
write_csv(
  gamma_est %>%
    select(GEAR, YEAR, MONTH, YEAR_MONTH, n_samples, N_total, S_obs, Chat_obs,
           C_target, extrapolated, coverage_fail,
           q0, q1, q2, se_q0, se_q1, se_q2, lcl_q0, lcl_q1, lcl_q2, ucl_q0, ucl_q1, ucl_q2),
  diag_cells_path
)
msg(icon("🧾 "), "Saved: ", diag_cells_path)

diag_summary <- gamma_est %>%
  summarise(
    n_cells       = n(),
    n_empty       = sum(N_total == 0),
    prop_empty    = mean(N_total == 0),
    n_fail_cov    = sum(coverage_fail),
    prop_fail_cov = mean(coverage_fail),
    median_chat   = median(Chat_obs),
    p10_chat      = as.numeric(quantile(Chat_obs, 0.10)),
    p90_chat      = as.numeric(quantile(Chat_obs, 0.90))
  )

msg("\n=== Coverage diagnostics summary (gamma cells) ===")
print(diag_summary)

plot_df <- gamma_est %>%
  select(GEAR, YEAR_MONTH, q0, q1, q2) %>%
  pivot_longer(cols = c(q0, q1, q2), names_to = "q", values_to = "gamma") %>%
  mutate(date = as.Date(paste0(YEAR_MONTH, "-01")))

p <- ggplot(plot_df, aes(x = date, y = gamma)) +
  geom_line() +
  facet_grid(q ~ GEAR, scales = "free_y") +
  labs(
    x = "Year-month",
    y = "Coverage-standardized gamma diversity",
    title = paste0("Gamma diversity time series (C_target=", sprintf("%.3f", C_target), ")")
  )

plot_path <- file.path(diag_dir, "b7_gamma_q012_timeseries__by_gear.png")
ggsave(plot_path, p, width = 12, height = 7, dpi = 300)
msg(icon("📈 "), "Saved: ", plot_path)

# ─────────────────────────────────────────────────────────────────────────────
# 12) Save outputs (derived data)
# ─────────────────────────────────────────────────────────────────────────────
out_rds <- file.path(ts_dir, "gamma_q012.RDS")
out_csv <- file.path(ts_dir, "gamma_q012.csv")

saveRDS(gamma_est, out_rds)
write_csv(gamma_est, out_csv)

msg(icon("💾 "), "Saved: ", out_rds)
msg(icon("💾 "), "Saved: ", out_csv)

msg("\nNOTE:")
msg("  - Sentinel rows (SPECIES_CODE = NA, catch = 0) are excluded from abundance vectors.")
msg("  - Time cells with effort but N_total = 0 are kept and assigned q0=q1=q2=0, Chat_obs=0.")
msg("  - Abundances are rounded to integers before iNEXT::estimateD() to avoid Linux bootstrap crashes.")
msg("\nOK: b7_compute_hill_gamma_q012.R completed successfully")