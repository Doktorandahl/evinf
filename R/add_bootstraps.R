#' Add more bootstrap replicates to a fitted model
#'
#' Fits \code{n} additional bootstrap replicates with the same block structure
#' and appends them, so you don't have to refit everything when the original
#' number turns out to be too few.
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model (with bootstraps).
#' @param n Number of additional bootstrap replicates.
#' @param boot_seed Optional RNG seed for this batch (recorded in
#'   \code{object$boot_seeds}).
#' @param multicore,ncores Passed to the parallel backend.
#' @param verbose Show a progress bar.
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
add_bootstraps <- function(object, n, boot_seed = NULL, multicore = FALSE,
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

  i <- "temp_iter"
  runner <- if (inherits(object, "evzinb")) bootrun_evzinb else bootrun_evinb
  block2 <- object$block_vec
  start <- length(object$bootstraps)

  be <- evinf_setup_backend(multicore, ncores)
  on.exit(be$stop(), add = TRUE)

  new_boots <- evinf_progress_run(n, verbose, function(p) {
    foreach::foreach(i = 1:n, .options.RNG = boot_seed, .packages = "evinf") %dorng% {
      res <- try(runner(object, block2))
      p()
      res
    }
  })
  names(new_boots) <- paste0("bootstrap_", start + seq_len(n))

  object$bootstraps <- c(object$bootstraps, new_boots)
  object$boot_seeds <- c(object$boot_seeds, list(boot_seed))
  object
}
