#' Add more bootstrap replicates to a fitted model
#'
#' Fits \code{n} additional bootstrap replicates with the same block structure
#' and appends them, so you don't have to refit everything when the original
#' number turns out to be too few.
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model (with bootstraps).
#' @param n Number of additional bootstrap replicates.
#' @param boot_seed RNG seed for this batch. Each batch needs its own seed:
#'   the bootstrap stream is fully determined by the seed, so reusing one that
#'   the model (or an earlier \code{add_bootstraps()} call) already used would
#'   silently duplicate those draws and is an error. When \code{NULL} a fresh
#'   seed is drawn and recorded in \code{object$boot_seeds}.
#' @inheritParams evzinb
#' @param verbose Show a progress bar.
#'
#' @inheritSection evzinb Parallel processing
#'
#' @return The model with \code{n} more replicates in \code{object$bootstraps}
#'   (names continue \code{bootstrap_<k>}).
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' model <- add_bootstraps(model, 5)
#' }
add_bootstraps <- function(object, n, boot_seed = NULL, multicore = NULL,
                           ncores = NULL, verbose = FALSE) {
  if (!inherits(object, c("evzinb", "evinb"))) {
    stop("`object` must be a fitted evzinb / evinb model.", call. = FALSE)
  }
  if (is.null(object$bootstraps)) {
    stop("`object` has no bootstraps; refit with bootstrap = TRUE.", call. = FALSE)
  }
  n <- as.integer(n)
  if (is.na(n) || n < 1L) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }

  # Reusing a seed replays the exact same L'Ecuyer stream, silently duplicating
  # the existing draws (audit N4). Draw one when NULL so it is always recorded.
  used <- unlist(object$boot_seeds)
  if (!is.null(boot_seed) && boot_seed %in% used) {
    stop("`boot_seed` = ", boot_seed, " was already used for this model's ",
         "bootstraps; reusing it would duplicate existing draws. Pass a ",
         "different seed.", call. = FALSE)
  }
  if (is.null(boot_seed)) {
    repeat {
      boot_seed <- sample.int(.Machine$integer.max, 1L)
      if (!boot_seed %in% used) break
    }
  }

  runner <- if (inherits(object, "evzinb")) bootrun_evzinb else bootrun_evinb
  block2 <- object$block_vec
  start <- length(object$bootstraps)
  boot_spec <- evinf_boot_spec(object)

  new_boots <- evinf_with_plan(multicore, ncores, {
    evinf_pmap(
      seq_len(n),
      function(i, spec, blk) try(runner(spec, blk)),
      spec = boot_spec, blk = block2,
      seed = boot_seed, label = "bootstrap", verbose = verbose
    )
  })
  names(new_boots) <- paste0("bootstrap_", start + seq_len(n))

  object$bootstraps <- c(object$bootstraps, new_boots)
  object$boot_seeds <- c(object$boot_seeds, list(boot_seed))
  n_c_bnd <- evinf_c_boundary_count(object, object$bootstraps)
  object$n_c_on_boundary <- n_c_bnd
  if (n_c_bnd > 0L) {
    warning("C_EV reached the boundary of the candidate range in ", n_c_bnd,
            " of ", length(object$bootstraps), " bootstrap replicates; ",
            "consider widening c.lim.", call. = FALSE)
  }
  object
}
