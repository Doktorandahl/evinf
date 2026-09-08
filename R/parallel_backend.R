# Internal parallel-backend management.
#
# Every exported function that may run a foreach loop calls evinf_setup_backend()
# once and registers its `stop` on exit. This guarantees that
#   * a parallel cluster is always torn down again (no leaked workers), and
#   * when multicore = FALSE a sequential backend is registered, so foreach's
#     "no parallel backend registered" warning never fires.
#
# Used by evzinb(), evinb(), compare_models(), lr_test(), predict.evzinb(),
# predict.evinb() and the quantile helpers.

evinf_setup_backend <- function(multicore = FALSE, ncores = NULL) {
  if (isTRUE(multicore)) {
    if (is.null(ncores)) {
      ncores <- max(1L, parallel::detectCores() - 1L)
    }
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      stopper <- function() {
        try(parallel::stopCluster(cl), silent = TRUE)
        foreach::registerDoSEQ()
      }
    } else {
      cl <- NULL
      doParallel::registerDoParallel(cores = ncores)
      stopper <- function() {
        try(doParallel::stopImplicitCluster(), silent = TRUE)
        foreach::registerDoSEQ()
      }
    }
    list(operator = foreach::`%dopar%`, stop = stopper, ncores = ncores,
         cl = cl, multicore = TRUE)
  } else {
    # Registering the sequential backend silences the %dopar%/%dorng%
    # "no parallel backend registered" warning without changing results.
    foreach::registerDoSEQ()
    list(operator = foreach::`%do%`, stop = function() invisible(NULL),
         ncores = 1L, cl = NULL, multicore = FALSE)
  }
}
