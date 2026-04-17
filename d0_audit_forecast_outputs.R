# =============================================================================
# d0_audit_forecast_outputs.R
#
# Phase 1 audit: schema + fold alignment across forecast outputs
# Base R only.
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- PATHS loader (your standard approach) ----
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

find_paths_file <- function(start_dir, max_up = 8) {
  d <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  for (i in 0:max_up) {
    candidate <- file.path(d, "01_Analysis", "PATHS.RData")
    if (file.exists(candidate)) return(candidate)
    d2 <- dirname(d); if (identical(d2, d)) break
    d <- d2
  }
  NA_character_
}

load_PATHS <- function() {
  pf <- find_paths_file(get_script_dir(), max_up = 8)
  if (is.na(pf)) stop("Could not find 01_Analysis/PATHS.RData walking up from script dir.")
  e <- new.env(parent = emptyenv())
  load(pf, envir = e)
  if (!exists("PATHS", envir = e)) stop("PATHS not found inside PATHS.RData: ", pf)
  get("PATHS", envir = e)
}

PATHS <- load_PATHS()

# ---- Helpers ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

extract_pred_df <- function(obj) {
  if (is.data.frame(obj)) return(obj)
  if (is.list(obj)) {
    if (!is.null(obj$predictions) && is.data.frame(obj$predictions)) return(obj$predictions)
    # fallback: first data.frame in list
    idx <- which(vapply(obj, is.data.frame, logical(1)))
    if (length(idx) >= 1) return(obj[[idx[1]]])
  }
  NULL
}

# Standardize just enough for auditing keys
ensure_cols <- function(df) {
  nms <- names(df)
  
  if (!("y_obs" %in% nms) && ("y_true" %in% nms)) df$y_obs <- df$y_true
  if (!("y_hat" %in% nms) && ("y_pred" %in% nms)) df$y_hat <- df$y_pred
  
  if (!("horizon" %in% nms) && ("h" %in% nms)) df$horizon <- df$h
  if (!("test_date" %in% nms) && ("date" %in% nms)) df$test_date <- df$date
  
  if (!("model" %in% nms) && ("model_id" %in% nms)) df$model <- df$model_id
  
  # PI names
  if (!("lower80" %in% nms) && ("lo80" %in% nms)) df$lower80 <- df$lo80
  if (!("upper80" %in% nms) && ("hi80" %in% nms)) df$upper80 <- df$hi80
  if (!("lower95" %in% nms) && ("lo95" %in% nms)) df$lower95 <- df$lo95
  if (!("upper95" %in% nms) && ("hi95" %in% nms)) df$upper95 <- df$hi95
  
  # ok/reason defaults
  if (!("ok" %in% names(df))) df$ok <- TRUE
  if (!("reason" %in% names(df))) df$reason <- NA_character_
  
  df
}

date_coerce <- function(x) {
  if (inherits(x, "Date")) return(x)
  suppressWarnings(as.Date(x))
}

summarize_pred <- function(df) {
  df <- ensure_cols(df)
  
  df$test_date <- date_coerce(df$test_date)
  if ("origin_date" %in% names(df)) df$origin_date <- date_coerce(df$origin_date)
  
  # core ranges
  horizons <- sort(unique(df$horizon))
  models   <- sort(unique(as.character(df$model)))
  gears    <- sort(unique(as.character(df$GEAR %||% df$gear)))
  qvals    <- sort(unique(as.character(df$q %||% df$Q)))
  
  # alpha/gamma
  has_bay <- "BAY" %in% names(df)
  bays <- if (has_bay) sort(unique(as.character(df$BAY))) else "ALL"
  
  ok_rate <- mean(as.logical(df$ok), na.rm = TRUE)
  n_ok <- sum(as.logical(df$ok), na.rm = TRUE)
  n_total <- nrow(df)
  
  # missingness on scoring fields
  y_obs_na <- sum(!is.finite(df$y_obs))
  y_hat_na <- sum(!is.finite(df$y_hat))
  
  # PI availability
  pi80_avail <- sum(is.finite(df$lower80) & is.finite(df$upper80))
  pi95_avail <- sum(is.finite(df$lower95) & is.finite(df$upper95))
  
  # duplicates on key fields (per model)
  key_fields <- c(if (has_bay) "BAY" else NULL, "GEAR", "q", "horizon",
                  if ("origin_date" %in% names(df)) "origin_date" else NULL,
                  if ("origin_i" %in% names(df)) "origin_i" else NULL,
                  "test_date",
                  if ("test_i" %in% names(df)) "test_i" else NULL,
                  "model")
  key_fields <- key_fields[key_fields %in% names(df)]
  key <- do.call(paste, c(df[key_fields], sep="||"))
  n_dup <- sum(duplicated(key))
  
  list(
    n_total = n_total,
    horizons = paste(horizons, collapse=","),
    models = paste(models, collapse=","),
    gears = paste(gears, collapse=","),
    q = paste(qvals, collapse=","),
    bays = if (has_bay) paste(head(bays, 10), collapse=",") else "ALL",
    ok_rate = ok_rate,
    n_ok = n_ok,
    y_obs_na = y_obs_na,
    y_hat_na = y_hat_na,
    pi80_avail = pi80_avail,
    pi95_avail = pi95_avail,
    n_dup = n_dup,
    key_fields = paste(key_fields, collapse=",")
  )
}

make_alignment_key <- function(df, family) {
  df <- ensure_cols(df)
  
  # Prefer index-based keys if present (best for fold-grid identity)
  has_idx <- all(c("origin_i", "test_i") %in% names(df))
  
  # Coerce to stable types for paste()
  to_chr <- function(x) as.character(x)
  
  if (has_idx) {
    if (family == "gamma") {
      fields <- c("GEAR", "q", "horizon", "origin_i", "test_i")
    } else {
      fields <- c("BAY", "GEAR", "q", "horizon", "origin_i", "test_i")
    }
    fields <- fields[fields %in% names(df)]
    tmp <- df[, fields, drop = FALSE]
    tmp[] <- lapply(tmp, to_chr)
    return(do.call(paste, c(tmp, sep = "||")))
  }
  
  # Fallback: date-based keys (older objects)
  df$test_date <- date_coerce(df$test_date)
  if ("origin_date" %in% names(df)) df$origin_date <- date_coerce(df$origin_date)
  
  if (family == "gamma") {
    fields <- c("GEAR", "q", "horizon",
                if ("origin_date" %in% names(df)) "origin_date" else NULL,
                "test_date")
  } else {
    fields <- c("BAY", "GEAR", "q", "horizon",
                if ("origin_date" %in% names(df)) "origin_date" else NULL,
                "test_date")
  }
  fields <- fields[fields %in% names(df)]
  tmp <- df[, fields, drop = FALSE]
  tmp[] <- lapply(tmp, to_chr)
  do.call(paste, c(tmp, sep = "||"))
}

# ---- File registry (canonical locations) ----
files <- list(
  gamma_baseline = list(
    family="gamma",
    path=file.path(PATHS$results_dir, "04_Models", "c4p_baseline_forecasts_gamma_parallel", "baseline_forecasts_gamma.rds")
  ),
  gamma_ssm = list(
    family="gamma",
    path=file.path(PATHS$results_dir, "04_Models", "c5_state_space_models_forecasts", "ssm_fits_gamma.rds")
  ),
  gamma_arimax = list(
    family="gamma",
    path=file.path(PATHS$results_dir, "04_Models", "c6_gamma_arimax_from_env_index", "gamma_arimax_forecasts.rds")
  ),
  gamma_ssm_xreg = list(
    family="gamma",
    path=file.path(PATHS$results_dir, "04_Models",
                   "c11_gamma_state_space_with_env_covariates_parallel",
                   "ssm_xreg_forecasts_gamma_by_gear_q.rds")
                   #"ssm_xreg_forecasts_gamma.rds")
  ),
  
  alpha_baseline = list(
    family="alpha",
    path=file.path(PATHS$results_dir, "04_Models",
                   "c7_alpha_forecasting_multilevel_or_bay_replication",
                   "alpha_forecasts_by_bay.rds")
  ),
  alpha_ssm = list(
    family="alpha",
    path=file.path(PATHS$results_dir, "04_Models",
                   "c8_alpha_state_space_models_forecasts_parallel",
                   "ssm_forecasts_alpha_by_bay.rds")
  ),
  alpha_ssm_xreg = list(
    family="alpha",
    path=file.path(PATHS$results_dir, "04_Models",
                   "c10_alpha_state_space_with_env_covariates_parallel",
                   "ssm_xreg_forecasts_alpha_by_bay.rds")
  ),
  alpha_arimax = list(
    family="alpha",
    path=file.path(PATHS$results_dir, "04_Models",
                   "c9_alpha_arimax_from_env_index_parallel",
                   "arimax_forecasts_alpha_by_bay.rds")
  )
)





# ---- Load + summarize ----
objs <- list()
preds <- list()
meta <- data.frame(
  name=character(), family=character(), path=character(),
  exists=logical(), notes=character(),
  stringsAsFactors=FALSE
)

for (nm in names(files)) {
  p <- files[[nm]]$path
  fam <- files[[nm]]$family
  ex <- file.exists(p)
  
  if (!ex) {
    meta <- rbind(meta, data.frame(name=nm, family=fam, path=p, exists=FALSE, notes="MISSING", stringsAsFactors=FALSE))
    next
  }
  
  obj <- readRDS(p)
  df <- extract_pred_df(obj)
  if (is.null(df)) {
    meta <- rbind(meta, data.frame(name=nm, family=fam, path=p, exists=TRUE, notes="No data.frame found", stringsAsFactors=FALSE))
    next
  }
  
  objs[[nm]] <- obj
  preds[[nm]] <- df
  
  s <- summarize_pred(df)
  meta <- rbind(meta, data.frame(
    name=nm, family=fam, path=p, exists=TRUE,
    notes=paste0(
      "n=", s$n_total,
      " | horizons=", s$horizons,
      " | models=", s$models,
      " | ok_rate=", round(s$ok_rate, 4),
      " | y_obs_na=", s$y_obs_na,
      " | y_hat_na=", s$y_hat_na,
      " | PI80=", s$pi80_avail,
      " | PI95=", s$pi95_avail,
      " | dup=", s$n_dup
    ),
    stringsAsFactors=FALSE
  ))
}

cat("\n============================\nPhase 1A: File presence + quick summaries\n============================\n")
print(meta, row.names = FALSE)

# ---- Check target provenance tags if present ----
cat("\n============================\nPhase 1B: Target provenance tags (eval_target_col / target_scale)\n============================\n")
for (nm in names(objs)) {
  obj <- objs[[nm]]
  tags <- c(
    eval_target_col = if (is.list(obj)) obj$eval_target_col %||% NA else NA,
    target_scale    = if (is.list(obj)) obj$target_scale %||% NA else NA
  )
  cat(nm, ":", paste(names(tags), unname(tags), sep="=", collapse=" | "), "\n")
}

# ---- Fold alignment across models (within family) ----
cat("\n============================\nPhase 1C: Fold-grid alignment checks\n============================\n")

align_family <- function(names_vec, family) {
  keys <- lapply(names_vec, function(nm) unique(make_alignment_key(preds[[nm]], family)))
  names(keys) <- names_vec
  # pairwise differences vs first
  ref <- keys[[1]]
  out <- data.frame(model=names_vec, n_keys=sapply(keys, length), diff_vs_ref=NA_integer_)
  for (i in seq_along(keys)) {
    out$diff_vs_ref[i] <- length(setdiff(keys[[i]], ref)) + length(setdiff(ref, keys[[i]]))
  }
  out
}

gamma_names <- intersect(names(preds), c("gamma_baseline","gamma_ssm","gamma_arimax"))
alpha_names <- intersect(names(preds), c("alpha_baseline","alpha_ssm","alpha_ssm_xreg","alpha_arimax"))

if (length(gamma_names) >= 2) {
  cat("\n-- Gamma alignment (key excludes model) --\n")
  print(align_family(gamma_names, "gamma"), row.names = FALSE)
}

if (length(alpha_names) >= 2) {
  cat("\n-- Alpha alignment (key excludes model) --\n")
  print(align_family(alpha_names, "alpha"), row.names = FALSE)
}

# ---- Report PI availability by model within each file ----
cat("\n============================\nPhase 1D: PI availability by model (within each file)\n============================\n")
for (nm in names(preds)) {
  df <- ensure_cols(preds[[nm]])
  mod <- as.character(df$model)
  pi80 <- is.finite(df$lower80) & is.finite(df$upper80)
  pi95 <- is.finite(df$lower95) & is.finite(df$upper95)
  tab <- aggregate(cbind(pi80=pi80, pi95=pi95) ~ mod, data=data.frame(mod=mod, pi80=pi80, pi95=pi95), FUN=function(x) mean(x, na.rm=TRUE))
  cat("\n", nm, "\n", sep="")
  print(tab, row.names = FALSE)
}

cat("\nDONE.\n")