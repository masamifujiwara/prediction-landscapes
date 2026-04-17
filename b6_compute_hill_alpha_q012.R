# =============================================================================
# Script: b6_compute_hill_alpha_q012.R
#
# Purpose (project directory convention):
#   Compute α diversity (q = 0, 1, 2) per (BAY_ID × YEAR × MONTH × GEAR)
#   using COVERAGE-BASED standardization (Hill numbers).
#
# Key requirements:
#   - Sentinel rows: SPECIES_CODE = NA, catch = 0 are EXCLUDED from abundance vectors.
#   - Cells with sampling effort but zero organisms yield α=0, Chat_obs=0.
#
# Linux stability notes:
#   - iNEXT::estimateD(datatype="abundance") assumes INTEGER abundances.
#     Non-integer pooled catch can crash during bootstrap (rmultinom).
#   - This script sanitizes pooled abundances (finite, >0, round->integer) and
#     avoids ifelse() on list-columns.
#   - Optional preflight probe identifies cells that fail estimateD() before long runs.
#
# Inputs (derived data):
#   - 01_Analysis/Data/b4_join_organisms_to_station/bs_events.rds
#   - 01_Analysis/Data/b4_join_organisms_to_station/bt_events.rds
#   - 01_Analysis/Data/b3_harmonize_station_metadata/station_core.rds
#   - 01_Analysis/Data/b5_define_spatiotemporal_units/index_keys.rds
#
# Outputs (derived data):
#   - 01_Analysis/Data/b6_compute_hill_alpha_q012/alpha_q012.rds
#   - 01_Analysis/Data/b6_compute_hill_alpha_q012/alpha_q012.csv
#
# Diagnostics:
#   - 02_Results/03_Diagnostics/b6_compute_hill_alpha_q012/
#       * b6_alpha_coverage_diagnostics__cells.csv
#       * b6_alpha_cells_that_fail_estimateD.csv (if any)
#       * b6_alpha_q012_timeseries_examples.png
# =============================================================================

rm(list = ls())

set.seed(101)

# ─────────────────────────────────────────────────────────────────────────────
# 0) Minimal messaging (avoid emoji in nohup logs)
# ─────────────────────────────────────────────────────────────────────────────
use_emoji <- FALSE
msg <- function(..., .nl = TRUE) {
  txt <- paste0(..., collapse = "")
  if (.nl) cat(txt, "\n") else cat(txt)
}
icon <- function(x) if (use_emoji) x else ""

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ─────────────────────────────────────────────────────────────────────────────
# 1) Bootstrap: setwd() to script folder + load PATHS
# ─────────────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  
  # 1) RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) {
      return(dirname(normalizePath(p)))
    }
  }
  
  # 2) Rscript (works with nohup)
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) {
      return(dirname(normalizePath(p)))
    }
  }
  
  # 3) Fallback: current working directory
  # (Safe for nohup if you cd into 01_Analysis first)
  return(normalizePath(getwd()))
}

script_dir <- get_script_dir()
setwd(script_dir)

cat("Working directory:", normalizePath(getwd()), "\n")
msg(icon("📁 "), "Working directory: ", normalizePath(getwd(), winslash = "/"))

get_project_root <- function() {
  here <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  normalizePath(file.path(here, ".."), winslash = "/", mustWork = TRUE)
}
project_dir <- get_project_root()

paths_file <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)  # loads PATHS
msg("Project directory: ", PATHS$project_dir)

options(scipen = 999, warn = 1)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Packages
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
# 3) Parameters
# ─────────────────────────────────────────────────────────────────────────────
nboot <- 100
conf  <- 0.95

C_cap      <- 0.90
C_quantile <- 0.10
C_min      <- 0.50

example_bays <- c(1, 4, 8)

# Preflight probe (recommended on Linux)
run_preflight_probe <- TRUE
nboot_probe <- 1

# If TRUE, stop if any non-integer abundances are detected before rounding
strict_noninteger_check <- FALSE

# ─────────────────────────────────────────────────────────────────────────────
# 4) Paths + load derived inputs
# ─────────────────────────────────────────────────────────────────────────────
in_bs_events <- file.path(PATHS$data_dir, "b4_join_organisms_to_station", "bs_events.rds")
in_bt_events <- file.path(PATHS$data_dir, "b4_join_organisms_to_station", "bt_events.rds")
in_station   <- file.path(PATHS$data_dir, "b3_harmonize_station_metadata", "station_core.rds")
in_keys      <- file.path(PATHS$data_dir, "b5_define_spatiotemporal_units", "index_keys.rds")

stopifnot(file.exists(in_bs_events), file.exists(in_bt_events), file.exists(in_station), file.exists(in_keys))

bs_events <- readRDS(in_bs_events)
bt_events <- readRDS(in_bt_events)
station_core <- readRDS(in_station)
index_keys <- readRDS(in_keys)

station_bs_core <- station_core$station_bs_core
station_bt_core <- station_core$station_bt_core

valid_bays <- index_keys$keys$valid_bays %||% 1:8

# Station effort table for alpha grid (includes zero-catch stations, excludes no-effort months)
station_all <- bind_rows(station_bs_core, station_bt_core) %>%
  mutate(
    BAY_ID = as.integer(MAJOR_AREA),
    YEAR   = as.integer(YEAR_i),
    MONTH  = as.integer(MONTH_i)
  ) %>%
  filter(BAY_ID %in% valid_bays) %>%
  select(SAMPLE_ID, GEAR, BAY_ID, YEAR, MONTH)

stopifnot(!anyDuplicated(station_bs_core$SAMPLE_ID))
stopifnot(!anyDuplicated(station_bt_core$SAMPLE_ID))

# Output folders (allowed via PATHS)
out_dir  <- file.path(PATHS$data_dir, "b6_compute_hill_alpha_q012")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

diag_dir <- PATHS$diag_dir  # 02_Results/03_Diagnostics/<script_stem>/ already created by PATHS
# (If you prefer b6 to write diagnostics into its own folder name rather than script_stem,
#  tell me and we’ll adjust PATHS logic. For now we follow the rule.)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Build pooled abundance vectors per α cell (BAY_ID×GEAR×YEAR×MONTH)
#    IMPORTANT: remove sentinel rows (SPECIES_CODE = NA) before pooling
# ─────────────────────────────────────────────────────────────────────────────
prep_events <- function(events_df, gear_label, valid_bays) {
  events_df %>%
    mutate(
      GEAR  = gear_label,
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

cell_species <- events_all %>%
  filter(!is.na(SPECIES_CODE)) %>%   # sentinel removed here
  group_by(BAY_ID, GEAR, YEAR, MONTH, SPECIES_CODE) %>%
  summarise(abund = sum(catch, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(abund) & abund > 0)

cell_n_samples <- station_all %>%
  count(BAY_ID, GEAR, YEAR, MONTH, name = "n_samples")

alpha_grid <- cell_n_samples %>%
  mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH)) %>%
  arrange(BAY_ID, GEAR, YEAR, MONTH)

# ─────────────────────────────────────────────────────────────────────────────
# 6) Utilities
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
           "If catch is not a count, datatype='abundance' may be inappropriate.",
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

estimate_q012_at_coverage <- function(x, C_target, nboot, conf,
                                      strict_noninteger_check = FALSE) {
  
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
    stop("estimateD() output does not contain an order column. Columns are:\n  - ",
         paste(cn, collapse = "\n  - "), call. = FALSE)
  }
  order_col <- order_col[1]
  
  est_candidates <- c("qD", "Estimate", "Diversity", "qD.Estimate", "Estimate.qD")
  est_col <- est_candidates[est_candidates %in% cn]
  if (length(est_col) == 0) {
    stop("estimateD() output does not contain an estimate column. Columns are:\n  - ",
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
    if (length(idx) == 0) {
      return(list(estimate = NA_real_, se = NA_real_, lcl = NA_real_, ucl = NA_real_))
    }
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
# 7) Compute observed coverage per cell, define global C_target
# ─────────────────────────────────────────────────────────────────────────────
cell_vectors <- cell_species %>%
  group_by(BAY_ID, GEAR, YEAR, MONTH) %>%
  summarise(
    abund_vec = list(setNames(abund, SPECIES_CODE)),
    N_total   = sum(abund),
    S_obs     = n(),
    .groups   = "drop"
  )

alpha_cells <- alpha_grid %>%
  left_join(cell_vectors, by = c("BAY_ID", "GEAR", "YEAR", "MONTH")) %>%
  mutate(
    N_total = coalesce(N_total, 0),
    S_obs   = coalesce(S_obs, 0)
  ) %>%
  mutate(
    abund_vec = map2(abund_vec, N_total, ~{
      if (is.null(.x) || is.na(.y) || .y == 0) setNames(numeric(0), character(0)) else .x
    })
  ) %>%
  mutate(Chat_obs = map_dbl(abund_vec, safe_chat))

chat_nonempty <- alpha_cells$Chat_obs[alpha_cells$N_total > 0]
if (length(chat_nonempty) == 0) stop("No non-empty cells found. Check input events.", call. = FALSE)

C_target <- min(C_cap, as.numeric(stats::quantile(chat_nonempty, probs = C_quantile, na.rm = TRUE)))
C_target <- max(min(C_target, 0.999), 0.001)

msg(icon("🎯 "), "Coverage target (C_target): ", sprintf("%.3f", C_target),
    " (cap=", C_cap, ", quantile=", C_quantile, ")")

alpha_cells <- alpha_cells %>%
  mutate(
    extrapolated  = (N_total > 0) & (C_target > Chat_obs),
    coverage_fail = (N_total == 0) | (Chat_obs < C_min),
    C_target      = C_target
  )

# ─────────────────────────────────────────────────────────────────────────────
# 7b) Diagnostics: non-integer check (informational)
# ─────────────────────────────────────────────────────────────────────────────
nonint_summary <- cell_vectors %>%
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
       "or revisit what 'catch' represents.", call. = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
# 7c) Optional preflight probe
# ─────────────────────────────────────────────────────────────────────────────
if (isTRUE(run_preflight_probe)) {
  
  msg("\n=== Preflight probe: estimateD per cell (nboot_probe=", nboot_probe, ") ===", sep = "")
  
  ok_vec <- logical(nrow(alpha_cells))
  
  for (i in seq_len(nrow(alpha_cells))) {
    
    if (alpha_cells$N_total[i] <= 0) {
      ok_vec[i] <- TRUE
      next
    }
    
    x <- sanitize_abundance(alpha_cells$abund_vec[[i]],
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
          " (BAY_ID=", alpha_cells$BAY_ID[i],
          ", GEAR=", alpha_cells$GEAR[i],
          ", ", alpha_cells$YEAR[i], "-", sprintf("%02d", alpha_cells$MONTH[i]), "): ",
          as.character(res)[1])
    }
    
    if (i %% 500 == 0) msg("... preflight progress: ", i, "/", nrow(alpha_cells))
  }
  
  if (any(!ok_vec)) {
    bad <- alpha_cells[!ok_vec, ] %>%
      select(BAY_ID, GEAR, YEAR, MONTH, YEAR_MONTH, n_samples, N_total, S_obs, Chat_obs, C_target)
    
    bad_path <- file.path(diag_dir, "b6_alpha_cells_that_fail_estimateD.csv")
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
# 8) Estimate α q0/q1/q2 at C_target per cell
# ─────────────────────────────────────────────────────────────────────────────
msg("\n=== Estimating alpha diversity per cell (nboot=", nboot, ") ===", sep = "")

est_list <- vector("list", nrow(alpha_cells))
for (i in seq_len(nrow(alpha_cells))) {
  est_list[[i]] <- estimate_q012_at_coverage(
    x = alpha_cells$abund_vec[[i]],
    C_target = C_target,
    nboot = nboot,
    conf = conf,
    strict_noninteger_check = strict_noninteger_check
  )
  if (i %% 500 == 0) msg("... estimate progress: ", i, "/", nrow(alpha_cells))
}

alpha_est <- bind_cols(alpha_cells, bind_rows(est_list)) %>%
  mutate(YEAR_MONTH = sprintf("%04d-%02d", YEAR, MONTH))

# ─────────────────────────────────────────────────────────────────────────────
# 9) Diagnostics summary + cell CSV
# ─────────────────────────────────────────────────────────────────────────────
diag_summary <- alpha_est %>%
  summarise(
    n_cells       = n(),
    n_empty       = sum(N_total == 0),
    prop_empty    = mean(N_total == 0),
    n_fail_cov    = sum(coverage_fail),
    prop_fail_cov = mean(coverage_fail),
    median_chat   = stats::median(Chat_obs),
    p10_chat      = as.numeric(stats::quantile(Chat_obs, 0.10)),
    p90_chat      = as.numeric(stats::quantile(Chat_obs, 0.90))
  )

msg("\n=== Coverage diagnostics summary (alpha cells) ===")
print(diag_summary)

diag_cells_path <- file.path(diag_dir, "b6_alpha_coverage_diagnostics__cells.csv")
write_csv(
  alpha_est %>%
    select(BAY_ID, GEAR, YEAR, MONTH, YEAR_MONTH, n_samples, N_total, S_obs, Chat_obs,
           C_target, extrapolated, coverage_fail,
           q0, q1, q2, se_q0, se_q1, se_q2, lcl_q0, lcl_q1, lcl_q2, ucl_q0, ucl_q1, ucl_q2),
  diag_cells_path
)
msg(icon("🧾 "), "Saved: ", diag_cells_path)

# ─────────────────────────────────────────────────────────────────────────────
# 10) Quick plot examples for a few bays (q0/q1/q2)
# ─────────────────────────────────────────────────────────────────────────────
plot_df <- alpha_est %>%
  filter(BAY_ID %in% example_bays) %>%
  select(BAY_ID, GEAR, YEAR_MONTH, q0, q1, q2) %>%
  pivot_longer(cols = c(q0, q1, q2), names_to = "q", values_to = "alpha") %>%
  mutate(date = as.Date(paste0(YEAR_MONTH, "-01")))

p <- ggplot(plot_df, aes(x = date, y = alpha)) +
  geom_line() +
  facet_grid(q ~ BAY_ID + GEAR, scales = "free_y") +
  labs(
    x = "Year-month",
    y = "Coverage-standardized alpha diversity",
    title = paste0("Alpha diversity time series examples (C_target=", sprintf("%.3f", C_target), ")")
  )

plot_path <- file.path(diag_dir, "b6_alpha_q012_timeseries_examples.png")
ggsave(plot_path, p, width = 12, height = 8, dpi = 300)
msg(icon("📈 "), "Saved: ", plot_path)

# ─────────────────────────────────────────────────────────────────────────────
# 11) Save outputs -> 01_Analysis/Data/<script_stem>/
# ─────────────────────────────────────────────────────────────────────────────
out_rds <- file.path(out_dir, "alpha_q012.rds")
out_csv <- file.path(out_dir, "alpha_q012.csv")

saveRDS(alpha_est, out_rds)
write_csv(alpha_est, out_csv)

msg(icon("💾 "), "Saved: ", out_rds)
msg(icon("💾 "), "Saved: ", out_csv)

msg("\nNOTE:")
msg("  - Sentinel rows (SPECIES_CODE = NA, catch = 0) are excluded from abundance vectors.")
msg("  - Cells with sampling effort but N_total = 0 are kept and assigned q0=q1=q2=0, Chat_obs=0.")
msg("  - Abundances are rounded to integers before iNEXT::estimateD() to avoid Linux bootstrap crashes.")
msg("\nOK: b6_compute_hill_alpha_q012.R completed successfully")