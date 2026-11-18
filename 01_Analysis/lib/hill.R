# =============================================================================
# File:    lib/hill.R
# Author:  Masami Fujiwara
#
# Purpose:
#   Shared functions for coverage-standardized Hill numbers (q = 0, 1, 2) of
#   one sample of counts by species. Loaded after the standard opening block
#   and lib/helpers.R:
#     source(file.path(PATHS$analysis_dir, "lib", "hill.R"))
#
# Used by: b06 (estuary scale), b07 (coastwide), f07
#
# Functions:
#   safe_chat()                  estimated sample coverage of one abundance vector
#   sanitize_abundance()         keep finite positive counts, round to integers
#   estimate_q012_at_coverage()  iNEXT::estimateD at a given coverage; q0, q1, q2
#
# Method notes:
#   - safe_chat(): when no species is caught exactly twice (f2 = 0), the
#     Chao and Jost (2012) correction (n - 1)(f1 - 1) / ((n - 1)(f1 - 1) + 2)
#     is used instead of 1. This is the estimator iNEXT uses itself.
#   - estimate_q012_at_coverage(): with nboot = 0 it returns q0, q1 and q2
#     only. The confidence limits come from the bootstrap; the point
#     estimates do not depend on it. No standard-error columns are returned
#     (iNEXT 3.0.2 does not report a standard error).
#   - estimate_q012_at_coverage() also returns m_at_C_target, the sample size
#     (individuals) iNEXT needs to reach the target coverage. When it exceeds
#     twice the observed catch, iNEXT warns that q = 0 may be biased; b06
#     keeps those cells and flags them.
#   Needs iNEXT, tibble and dplyr (called with ::); nothing is attached.
#
# AI assistance: This file was written with the assistance of Claude Opus 5.5
# (Anthropic). The author reviewed and approved all code and comments.
# =============================================================================

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
  
  # Chao and Jost (2012); the f2 = 0 branch is their bias-corrected form.
  A <- if (f2 > 0) ((N - 1) * f1) / (((N - 1) * f1) + (2 * f2)) else
    ((N - 1) * (f1 - 1)) / (((N - 1) * (f1 - 1)) + 2)
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
    out <- tibble::tibble(q0 = 0, q1 = 0, q2 = 0, m_at_C_target = NA_real_)
    if (nboot > 0) {
      out <- dplyr::bind_cols(out, tibble::tibble(
        lcl_q0 = NA_real_, lcl_q1 = NA_real_, lcl_q2 = NA_real_,
        ucl_q0 = NA_real_, ucl_q1 = NA_real_, ucl_q2 = NA_real_))
    }
    return(out)
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
  
  lcl_col <- intersect(c("qD.LCL", "LCL", "LCI", "lcl", "Lower"), cn)
  ucl_col <- intersect(c("qD.UCL", "UCL", "UCI", "ucl", "Upper"), cn)
  
  lcl_col <- if (length(lcl_col) > 0) lcl_col[1] else NA_character_
  ucl_col <- if (length(ucl_col) > 0) ucl_col[1] else NA_character_
  
  ord <- suppressWarnings(as.numeric(est[[order_col]]))
  
  get_one <- function(q) {
    idx <- which(ord == q)
    if (length(idx) == 0) {
      return(list(estimate = NA_real_, lcl = NA_real_, ucl = NA_real_))
    }
    idx <- idx[1]
    list(
      estimate = as.numeric(est[[est_col]][idx]),
      lcl = if (!is.na(lcl_col)) as.numeric(est[[lcl_col]][idx]) else NA_real_,
      ucl = if (!is.na(ucl_col)) as.numeric(est[[ucl_col]][idx]) else NA_real_
    )
  }
  
  r0 <- get_one(0); r1 <- get_one(1); r2 <- get_one(2)
  
  m_col <- if ("m" %in% cn) as.numeric(est[["m"]][1]) else NA_real_
  out <- tibble::tibble(q0 = r0$estimate, q1 = r1$estimate, q2 = r2$estimate,
                        m_at_C_target = m_col)
  if (nboot > 0) {
    out <- dplyr::bind_cols(out, tibble::tibble(
      lcl_q0 = r0$lcl, lcl_q1 = r1$lcl, lcl_q2 = r2$lcl,
      ucl_q0 = r0$ucl, ucl_q1 = r1$ucl, ucl_q2 = r2$ucl))
  }
  out
}
