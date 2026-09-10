# Parallel execution and progress reporting for evinf (round 5).
#
# The package runs its bootstrap / refit loops through one helper, evinf_pmap(),
# built on furrr + future. Parallelism is controlled by the user's
# future::plan(); the multicore / ncores arguments are a convenience wrapper
# (evinf_with_plan()) that sets a plan temporarily.

#' Run a function over a list in parallel, with a progress bar
#'
#' The single point through which every parallel loop in the package goes.
#' Wraps \code{furrr::future_map()} with reproducible per-element RNG streams
#' (\code{seed}) and a \pkg{progressr} progress bar that updates once per
#' element.
#'
#' @param .x A list or vector to iterate over.
#' @param .f A function applied to each element of \code{.x} (plus \code{...}).
#' @param ... Additional arguments passed on to \code{.f} on every call.
#' @param seed An integer seed. \code{furrr} derives an independent L'Ecuyer
#'   stream for each element from it, identical regardless of the number of
#'   workers or the chunk size.
#' @param label Progress-bar message.
#' @param verbose If \code{TRUE}, always show a progress bar. Otherwise one is
#'   shown only when the user has enabled a global \pkg{progressr} handler.
#' @param chunk_size Number of elements handed to a worker at a time.
#'
#' @return A list, the element-wise results of \code{.f}.
#' @keywords internal
evinf_pmap <- function(.x, .f, ..., seed, label = "bootstrap", verbose = FALSE,
                       chunk_size = 1L) {
  if (!is.numeric(seed) || length(seed) != 1L) {
    stop("`seed` must be a single integer.", call. = FALSE)
  }
  opts <- furrr::furrr_options(
    seed = as.integer(seed), chunk_size = chunk_size, packages = "evinf"
  )
  run <- function(p) {
    furrr::future_map(
      .x,
      function(.e, ...) {
        out <- .f(.e, ...)
        p(message = label)
        out
      },
      ...,
      .options = opts
    )
  }
  if (!evinf_progress_active(verbose)) {
    return(run(function(...) invisible()))
  }
  res <- NULL
  progressr::with_progress({
    p <- progressr::progressor(along = .x)
    res <- run(p)
  })
  res
}

# Progress is shown when verbose = TRUE or the user has enabled a global
# progressr handler; otherwise the progressor is a no-op and there is no
# progress-bar overhead.
evinf_progress_active <- function(verbose) {
  if (isTRUE(verbose)) {
    return(TRUE)
  }
  isTRUE(tryCatch(progressr::handlers(global = NA), error = function(e) FALSE))
}

# One-per-session hint that setting a future::plan() is the idiomatic route.
evinf_plan_hint <- function() {
  if (isTRUE(getOption("evinf.plan_hint_shown"))) {
    return(invisible())
  }
  options(evinf.plan_hint_shown = TRUE)
  message(
    "evinf: `multicore = TRUE` sets a temporary future::plan(multisession). ",
    "For more control set the plan yourself once, e.g. ",
    "future::plan(future::multisession, workers = 8), and leave `multicore` ",
    "at its default."
  )
}

#' Evaluate an expression under a temporary \code{future::plan()}
#'
#' @param multicore \code{NULL} (default) leaves the current
#'   \code{future::plan()} untouched; \code{TRUE} sets
#'   \code{future::multisession} with \code{ncores} workers; \code{FALSE} sets
#'   \code{sequential}. The plan is always restored on exit.
#' @param ncores Number of workers for \code{multicore = TRUE}; default is one
#'   less than \code{parallelly::availableCores()}.
#' @param expr The expression to evaluate (lazily, after the plan is set).
#'
#' @return The value of \code{expr}.
#' @keywords internal
evinf_with_plan <- function(multicore = NULL, ncores = NULL, expr) {
  if (is.null(multicore)) {
    return(expr)
  }
  oplan <- future::plan()
  on.exit(future::plan(oplan), add = TRUE)
  if (isTRUE(multicore)) {
    workers <- if (is.null(ncores)) {
      max(1L, parallelly::availableCores() - 1L)
    } else {
      as.integer(ncores)
    }
    future::plan(future::multisession, workers = workers)
    evinf_plan_hint()
  } else {
    future::plan("sequential")
  }
  expr
}

# Flag a bootstrap replicate as degenerate -- a numerically valid fit that would
# nonetheless poison the bootstrap summaries. Checks, in order, and records the
# first matching reason (thresholds from evinf_control(): `alpha_floor`,
# `coef_limit`):
#   1. the EM did not converge;
#   2. a fitted linear-predictor coefficient or the NB dispersion is non-finite
#      or larger than `coef_limit` in absolute value;
#   3. the smallest fitted Pareto shape on the replicate's own resample is
#      non-finite or below `alpha_floor`.
# Sets $degenerate (a single logical) and $degenerate_reason (string or NA).
evinf_flag_degenerate <- function(boot, X.PL, control) {
  af <- control$alpha_floor %||% 0.001
  cl <- control$coef_limit %||% 50
  reason <- NA_character_

  if (isFALSE(boot$converge)) {
    reason <- "EM did not converge within max.no.em.steps"
  }

  if (is.na(reason)) {
    comps <- list(
      Beta.NB          = boot$coef$Beta.NB,
      Beta.multinom.ZC = boot$coef$Beta.multinom.ZC,   # NULL for evinb
      Beta.multinom.PL = boot$coef$Beta.multinom.PL,
      Beta.PL          = boot$coef$Beta.PL,
      Alpha.NB         = boot$coef$Alpha.NB
    )
    for (nm in names(comps)) {
      v <- comps[[nm]]
      if (is.null(v)) next
      if (any(!is.finite(v))) {
        reason <- sprintf("non-finite coefficient in %s", nm)
        break
      }
      hit <- which(abs(v) > cl)
      if (length(hit)) {
        lbl <- if (!is.null(names(v))) names(v)[hit[1]] else nm
        reason <- sprintf(
          "coefficient %s = %.3g in %s exceeds coef_limit (%g)",
          lbl, v[hit[1]], nm, cl)
        break
      }
    }
  }

  if (is.na(reason)) {
    min_a <- suppressWarnings(
      min(exp(as.numeric(cbind(1, X.PL) %*% boot$coef$Beta.PL))))
    if (!is.finite(min_a)) {
      reason <- "fitted Pareto shape is not finite"
    } else if (min_a < af) {
      reason <- sprintf(
        "smallest fitted Pareto shape %.3g < alpha_floor (%.3g)", min_a, af)
    }
  }

  boot$degenerate <- !is.na(reason)
  boot$degenerate_reason <- reason
  boot
}

# Trim a fitted evzinb / evinb object to the fields bootrun_*() reads, so that
# only a compact "spec" is shipped to parallel workers (not the raw data frame,
# the c-profile, the loglik trace or the fitted-value vectors).
evinf_boot_spec <- function(full_run) {
  spec <- list()
  spec$data <- list(
    y             = full_run$data$y,
    x.nb          = full_run$data$x.nb,
    x.pl          = full_run$data$x.pl,
    x.multinom.zc = full_run$data$x.multinom.zc,
    x.multinom.pl = full_run$data$x.multinom.pl
  )
  spec$offset_nb <- full_run$offset_nb
  spec$control   <- full_run$control
  spec$coef      <- full_run$coef
  spec$formulas  <- full_run$formulas
  spec$terms     <- full_run$terms
  spec$xlevels   <- full_run$xlevels
  spec$block_vec <- full_run$block_vec
  class(spec) <- class(full_run)
  spec
}

# Count the bootstrap replicates whose C_EV estimate sits on an endpoint of the
# C_EV candidate grid for the full sample (the unique observed response values
# inside control$c.lim). This is informational -- a warning, not a degeneracy
# criterion: excluding these replicates would bias the bootstrap distribution of
# C_EV inward.
evinf_c_boundary_count <- function(full_run, boots) {
  cl <- full_run$control$c.lim
  if (is.null(cl)) {
    return(0L)
  }
  uy <- sort(unique(full_run$data$y))
  grid <- uy[uy >= cl[1] & uy <= cl[2]]
  if (!length(grid)) {
    return(0L)
  }
  ends <- range(grid)
  c_boot <- vapply(boots, function(b) {
    if (inherits(b, "try-error")) NA_real_ else b$coef$C
  }, numeric(1))
  sum(c_boot %in% ends, na.rm = TRUE)
}
