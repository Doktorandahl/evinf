# Internal parallel-backend management.
#
# Every exported function that may run a foreach loop calls evinf_setup_backend()
# once and registers its `stop` on exit. This guarantees that
#   * a parallel cluster we created is always torn down again, and
#   * when multicore = FALSE a backend is available, so foreach's "no parallel
#     backend registered" warning never fires,
# while never tearing down a backend the user registered themselves (audit R0.4):
# if the user has their own doParallel / doFuture backend and multicore = FALSE,
# we leave it alone and simply run the loop on it with %do% semantics disabled
# (i.e. use %dopar% via their backend). We can only restore "sequential", never
# an arbitrary third-party backend, so we never touch one we did not create.
#
# Used by evzinb(), evinb(), compare_models(), lr_test(), predict.evzinb(),
# predict.evinb() and the quantile helpers.

evinf_setup_backend <- function(multicore = FALSE, ncores = NULL) {
  had_backend <- isTRUE(try(foreach::getDoParRegistered(), silent = TRUE))
  backend_name <- if (had_backend) foreach::getDoParName() else NULL
  user_backend <- had_backend && !identical(backend_name, "doSEQ")

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
    return(list(operator = foreach::`%dopar%`, stop = stopper, ncores = ncores,
                cl = cl, multicore = TRUE, user_backend = user_backend))
  }

  # multicore = FALSE
  if (user_backend) {
    # Run on the user's own backend; leave it exactly as we found it.
    return(list(operator = foreach::`%dopar%`, stop = function() invisible(NULL),
                ncores = foreach::getDoParWorkers(), cl = NULL, multicore = FALSE,
                user_backend = TRUE))
  }

  # No (real) backend registered: make sure a sequential one is, so %dopar% /
  # %dorng% do not warn. Only restore sequential on exit if we registered it
  # here (which is a no-op anyway once doSEQ is active).
  if (!had_backend) {
    foreach::registerDoSEQ()
  }
  list(operator = foreach::`%do%`, stop = function() invisible(NULL),
       ncores = 1L, cl = NULL, multicore = FALSE, user_backend = FALSE)
}
