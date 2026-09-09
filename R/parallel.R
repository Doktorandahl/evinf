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
