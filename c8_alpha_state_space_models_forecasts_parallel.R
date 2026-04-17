# =============================================================================
# Script: c8_alpha_state_space_models_forecasts_parallel.R  (REVISED to match c7)
#
# Purpose:
#   KFAS SSM forecasts for alpha (q=0,1,2) by BAY×GEAR×q with rolling-origin folds
#   on a GLOBAL complete monthly grid (to align fold indices with c7/c9).
#
# Output conventions (match revised c7):
#   - Save RDS as a list: predictions + eval_target_col + target_scale + script_stem + created
#   - predictions columns include: BAY,GEAR,q,horizon,origin_i,test_i,origin_date,test_date,
#     model,y_obs,y_hat,lower80,upper80,lower95,upper95,ok,reason
#   - Save to BOTH 01_Analysis/Data/<stem>/ and 02_Results/04_Models/<stem>/
# =============================================================================

rm(list = ls())

# -----------------------------
# Required bootstrap snippet
# -----------------------------
get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 1) {
    p <- sub("^--file=", "", file_arg)
    if (nzchar(p)) return(dirname(normalizePath(p, winslash = "/")))
  }
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) return(dirname(normalizePath(ofile, winslash = "/")))
  normalizePath(getwd(), winslash = "/")
}

script_dir <- get_script_dir()
setwd(script_dir)

# -----------------------------
# PATHS
# -----------------------------
project_dir <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
paths_file  <- file.path(project_dir, "01_Analysis", "PATHS.RData")
stopifnot(file.exists(paths_file))
load(paths_file)
stopifnot(exists("PATHS"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------
# Dependencies
# -----------------------------
need_pkgs <- c("dplyr", "tidyr", "readr", "KFAS", "future", "future.apply", "parallel")
missing_pkgs <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(KFAS)
  library(future)
  library(future.apply)
})

# -----------------------------
# Config (match c7 philosophy)
# -----------------------------
script_stem <- "c8_alpha_state_space_models_forecasts_parallel"

MODEL_NAME <- "ssm_ll_seas12"
FREQ <- 12
HORIZONS <- c(1, 12)

TRAINING_RULE <- "expanding"
MIN_TRAIN_MONTHS <- 120
FIXED_TRAIN_MONTHS <- 120
ORIGIN_STEP <- 1
HOLDOUT_LAST_N <- 0
MIN_TRAIN_NON_NA <- 24

TARGET_COL_PRIMARY <- "value_deseas"
TARGET_COL_FALLBACKS <- c("residual_component", "value")

TARGET_SCALE <- "deseasonalized"
PI_LEVELS <- c(80, 95)

MAXIT <- 200
QUIET <- TRUE

# Exclude Sabine Lake (Bay 1) going forward
DROP_BAYS <- c("1")
# -----------------------------
# Parallel plan (like your other scripts)
# -----------------------------
N_WORKERS <- min(6, max(1, parallel::detectCores(logical = FALSE) - 2))
future::plan(future::multisession, workers = N_WORKERS)
# -----------------------------
# Directories
# -----------------------------
dir_data_script  <- file.path(PATHS$analysis_dir, "Data", script_stem)
dir_logs_script  <- file.path(PATHS$analysis_dir, "Logs", script_stem)
dir_diag_script  <- file.path(PATHS$results_dir, "03_Diagnostics", script_stem)
dir_model_script <- file.path(PATHS$results_dir, "04_Models", script_stem)

dir.create(dir_data_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_logs_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_diag_script, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_model_script, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(dir_logs_script, paste0("runlog_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, split = TRUE)
cat("Script:", script_stem, "\nStarted:", format(Sys.time()), "\nWorkers:", N_WORKERS, "\n\n")

# Inputs
IN_ALPHA <- file.path(PATHS$analysis_dir, "Data", "c2_deseasonalize_alpha_gamma", "alpha_q012_deseas.rds")
stopifnot(file.exists(IN_ALPHA))

# Outputs
OUT_FC_DERIVED  <- file.path(dir_data_script,  "ssm_forecasts_alpha_by_bay.rds")
OUT_MET_DERIVED <- file.path(dir_data_script,  "ssm_metrics_alpha_by_bay.csv")
OUT_FC_MODEL    <- file.path(dir_model_script, "ssm_forecasts_alpha_by_bay.rds")
OUT_MET_MODEL   <- file.path(dir_model_script, "ssm_metrics_alpha_by_bay.csv")
OUT_FITLOG      <- file.path(dir_diag_script,  "ssm_fit_log.csv")





# -----------------------------
# Helpers
# -----------------------------
expect_cols <- function(df, required, name = "data") {
  miss <- setdiff(required, names(df))
  if (length(miss) > 0) stop(name, " missing required columns: ", paste(miss, collapse = ", "))
}

pick_target_col <- function(df) {
  cands <- c(TARGET_COL_PRIMARY, TARGET_COL_FALLBACKS)
  hit <- cands[cands %in% names(df)]
  if (length(hit) == 0) stop("No target column found. Tried: ", paste(cands, collapse = ", "))
  hit[1]
}

standardize_alpha_names <- function(alpha) {
  if (!("GEAR" %in% names(alpha))) {
    gcol <- intersect(names(alpha), c("gear", "Gear"))
    if (length(gcol) >= 1) alpha <- dplyr::rename(alpha, GEAR = all_of(gcol[1]))
  }
  if (!("q" %in% names(alpha))) {
    qcol <- intersect(names(alpha), c("Q", "hill_q", "order_q"))
    if (length(qcol) >= 1) alpha <- dplyr::rename(alpha, q = all_of(qcol[1]))
  }
  if (!("date" %in% names(alpha))) {
    if (all(c("YEAR", "MONTH") %in% names(alpha))) {
      alpha <- alpha %>% mutate(date = as.Date(sprintf("%04d-%02d-01", as.integer(YEAR), as.integer(MONTH))))
    } else if (all(c("year", "month") %in% names(alpha))) {
      alpha <- alpha %>% mutate(date = as.Date(sprintf("%04d-%02d-01", as.integer(.data$year), as.integer(.data$month))))
    } else {
      stop("alpha input has no date or (YEAR, MONTH) columns.")
    }
  }
  if (!("BAY" %in% names(alpha))) {
    bcol <- intersect(names(alpha), c("BAY_ID", "bay", "Bay", "BAYID"))
    if (length(bcol) >= 1) alpha <- dplyr::rename(alpha, BAY = all_of(bcol[1]))
  }
  alpha
}

collapse_dups_alpha <- function(df, value_col) {
  df %>%
    group_by(BAY, GEAR, q, date) %>%
    summarise(
      value = if (all(is.na(.data[[value_col]]))) NA_real_ else stats::median(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

complete_monthly_grid_global <- function(df_obs, global_start, global_end) {
  grid <- tibble(date = seq.Date(global_start, global_end, by = "month"))
  grid %>% left_join(df_obs, by = "date") %>% arrange(date)
}

build_folds_one_series <- function(n, horizons,
                                   training_rule = "expanding",
                                   min_train_months = 120,
                                   fixed_train_months = 120,
                                   origin_step = 1,
                                   holdout_last_n = 0) {
  
  n_eff <- if (holdout_last_n > 0) max(0, n - holdout_last_n) else n
  if (n_eff < (min_train_months + max(horizons) + 1)) {
    return(list(ok = FALSE, reason = "Too short for min_train + max_horizon", folds = list()))
  }
  
  folds <- list()
  for (h in horizons) {
    origins <- seq(from = min_train_months, to = (n_eff - h), by = origin_step)
    folds[[paste0("h", h)]] <- lapply(origins, function(t0) {
      train_idx <- if (training_rule == "expanding") 1:t0 else max(1, t0 - fixed_train_months + 1):t0
      list(horizon = h, origin_i = t0, test_i = t0 + h, train_idx = train_idx)
    })
  }
  list(ok = TRUE, reason = NA_character_, folds = folds)
}

fit_forecast_ssm <- function(y_train, h, freq = 12, pi_levels = c(80, 95)) {
  tryCatch({
    y_ts <- stats::ts(as.numeric(y_train), frequency = freq)
    
    model <- SSModel(
      y_ts ~ SSMtrend(degree = 1, Q = NA) +
        SSMseasonal(period = freq, sea.type = "dummy", Q = NA),
      H = NA
    )
    
    v0 <- stats::var(y_train, na.rm = TRUE)
    if (!is.finite(v0) || v0 <= 0) v0 <- 1
    init <- rep(log(v0), 3)
    
    fit <- fitSSM(model, inits = init, method = "BFGS",
                  control = list(maxit = MAXIT), quiet = QUIET)
    mod_f <- fit$model
    
    pred_mean <- predict(mod_f, n.ahead = h, interval = "none")
    mu_h <- as.numeric(pred_mean[h])
    
    get_int <- function(level_pct) {
      lev <- level_pct / 100
      pr <- predict(mod_f, n.ahead = h, interval = "prediction", level = lev)
      c(lwr = as.numeric(pr[h, "lwr"]), upr = as.numeric(pr[h, "upr"]))
    }
    
    int80 <- if (80 %in% pi_levels) get_int(80) else c(lwr = NA_real_, upr = NA_real_)
    int95 <- if (95 %in% pi_levels) get_int(95) else c(lwr = NA_real_, upr = NA_real_)
    
    list(ok = TRUE, mean = mu_h,
         lower80 = unname(int80["lwr"]), upper80 = unname(int80["upr"]),
         lower95 = unname(int95["lwr"]), upper95 = unname(int95["upr"]),
         reason = NA_character_)
  }, error = function(e) {
    list(ok = FALSE, mean = NA_real_,
         lower80 = NA_real_, upper80 = NA_real_,
         lower95 = NA_real_, upper95 = NA_real_,
         reason = conditionMessage(e))
  })
}

rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
mae  <- function(e) mean(abs(e), na.rm = TRUE)

coverage_rate <- function(y, lo, hi) {
  ok <- is.finite(y) & is.finite(lo) & is.finite(hi)
  if (!any(ok)) return(NA_real_)
  mean(y[ok] >= lo[ok] & y[ok] <= hi[ok])
}

# -----------------------------
# Load + prep alpha
# -----------------------------
alpha_obj <- readRDS(IN_ALPHA)
alpha <- alpha_obj
if (is.list(alpha_obj) && !is.data.frame(alpha_obj)) {
  idx <- which(vapply(alpha_obj, is.data.frame, logical(1)))
  if (length(idx) >= 1) alpha <- alpha_obj[[idx[1]]]
}
stopifnot(is.data.frame(alpha))

alpha <- standardize_alpha_names(alpha)
expect_cols(alpha, c("BAY", "GEAR", "q", "date"), "alpha")

target_col <- pick_target_col(alpha)
EVAL_TARGET_COL <- target_col
cat("Using target column:", target_col, "\n")

alpha2 <- alpha %>%
  mutate(
    BAY  = as.character(BAY),
    GEAR = as.character(GEAR),
    q    = as.character(q),
    date = as.Date(date)
  ) %>%
  filter(!(BAY %in% DROP_BAYS))

alpha_collapsed <- collapse_dups_alpha(alpha2, value_col = target_col)

GLOBAL_START <- min(alpha_collapsed$date, na.rm = TRUE)
GLOBAL_END   <- max(alpha_collapsed$date, na.rm = TRUE)

series_keys <- alpha_collapsed %>%
  distinct(BAY, GEAR, q) %>%
  arrange(BAY, GEAR, q)

# -----------------------------
# Build series tables + folds
# -----------------------------
series_list <- vector("list", nrow(series_keys))
folds_list  <- vector("list", nrow(series_keys))

for (i in seq_len(nrow(series_keys))) {
  key <- series_keys[i, ]
  df_i <- alpha_collapsed %>%
    filter(BAY == key$BAY, GEAR == key$GEAR, q == key$q) %>%
    arrange(date)
  
  df_g <- complete_monthly_grid_global(df_i, GLOBAL_START, GLOBAL_END)
  
  folds_i <- build_folds_one_series(
    n = nrow(df_g),
    horizons = HORIZONS,
    training_rule = TRAINING_RULE,
    min_train_months = MIN_TRAIN_MONTHS,
    fixed_train_months = FIXED_TRAIN_MONTHS,
    origin_step = ORIGIN_STEP,
    holdout_last_n = HOLDOUT_LAST_N
  )
  
  series_list[[i]] <- df_g
  folds_list[[i]]  <- folds_i
}

# -----------------------------
# Build job table (series×folds×h)
# -----------------------------
jobs <- list()
for (i in seq_len(nrow(series_keys))) {
  key <- series_keys[i, ]
  folds_i <- folds_list[[i]]
  if (!isTRUE(folds_i$ok)) next
  
  for (h_name in names(folds_i$folds)) {
    for (f in folds_i$folds[[h_name]]) {
      jobs[[length(jobs) + 1]] <- list(
        series_i = i,
        BAY = key$BAY,
        GEAR = key$GEAR,
        q = key$q,
        horizon = f$horizon,
        origin_i = f$origin_i,
        test_i = f$test_i,
        train_idx = f$train_idx
      )
    }
  }
}

cat("Jobs:", length(jobs), "\n")

# -----------------------------
# Run one job (standardized row)
# -----------------------------
run_one_job <- function(job) {
  df_g <- series_list[[job$series_i]]
  y <- df_g$value
  
  y_train <- y[job$train_idx]
  n_non_na <- sum(is.finite(y_train))
  y_obs <- y[job$test_i]
  
  # defaults
  ok_row <- TRUE
  reason <- NA_character_
  fc <- list(ok = FALSE, mean = NA_real_,
             lower80 = NA_real_, upper80 = NA_real_,
             lower95 = NA_real_, upper95 = NA_real_,
             reason = NA_character_)
  
  if (!is.finite(y_obs)) {
    ok_row <- FALSE
    reason <- "missing_target_month"
  }
  
  if (n_non_na < MIN_TRAIN_NON_NA) {
    ok_row <- FALSE
    reason <- paste0("insufficient_training_nonNA<", MIN_TRAIN_NON_NA)
  }
  
  # Fit only if we can score + have enough training
  if (ok_row && n_non_na >= MIN_TRAIN_NON_NA) {
    fc <- fit_forecast_ssm(y_train = y_train, h = job$horizon, freq = FREQ, pi_levels = PI_LEVELS)
    if (!isTRUE(fc$ok)) {
      ok_row <- FALSE
      reason <- fc$reason %||% "model_fit_failed"
    }
  }
  
  tibble(
    BAY = job$BAY, GEAR = job$GEAR, q = job$q,
    model = MODEL_NAME,
    horizon = job$horizon,
    origin_i = job$origin_i,
    test_i = job$test_i,
    origin_date = df_g$date[job$origin_i],
    test_date   = df_g$date[job$test_i],
    y_obs = y_obs,
    y_hat = fc$mean,
    lower80 = fc$lower80, upper80 = fc$upper80,
    lower95 = fc$lower95, upper95 = fc$upper95,
    ok = ok_row,
    reason = reason
  )
}

res <- future.apply::future_lapply(jobs, run_one_job, future.seed = TRUE)
pred_tbl <- bind_rows(res)

# Fit log
fit_log <- pred_tbl %>%
  transmute(BAY, GEAR, q, horizon, origin_i, test_i, test_date, ok, reason)
readr::write_csv(fit_log, OUT_FITLOG)

# Metrics (optional; ok to keep)
metrics_tbl <- pred_tbl %>%
  group_by(BAY, GEAR, q, horizon, model) %>%
  summarise(
    n = sum(is.finite(y_obs) & is.finite(y_hat) & ok),
    rmse = rmse(y_hat - y_obs),
    mae  = mae(y_hat - y_obs),
    coverage80 = coverage_rate(y_obs, lower80, upper80),
    coverage95 = coverage_rate(y_obs, lower95, upper95),
    frac_failed = mean(!ok),
    .groups = "drop"
  )

# Wrap like c7
save_obj <- list(
  predictions     = pred_tbl,
  eval_target_col = EVAL_TARGET_COL,
  target_scale    = TARGET_SCALE,
  script_stem     = script_stem,
  created         = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

saveRDS(save_obj, OUT_FC_DERIVED)
saveRDS(save_obj, OUT_FC_MODEL)
readr::write_csv(metrics_tbl, OUT_MET_DERIVED)
readr::write_csv(metrics_tbl, OUT_MET_MODEL)

cat("\nWrote:\n")
cat("  -", OUT_FC_DERIVED, "\n")
cat("  -", OUT_FC_MODEL, "\n")
cat("  -", OUT_MET_DERIVED, "\n")
cat("  -", OUT_MET_MODEL, "\n")
cat("  -", OUT_FITLOG, "\n")

future::plan(future::sequential)
cat("\nDone.\n")

sink()
cat("Log written to:", log_file, "\n")