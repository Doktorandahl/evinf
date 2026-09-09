# Progress reporting for the bootstrap loops (audit 4.9), via progressr.
#
# evinf_progress_run(n, verbose, fun) evaluates fun(p) where `p` is a function to
# call once per completed iteration. Progress is shown when verbose = TRUE or the
# user has enabled a global progressr handler; otherwise `p` is a no-op and there
# is no overhead. Works with the foreach / %dorng% backend.

evinf_progress_active <- function(verbose) {
  if (isTRUE(verbose)) {
    return(TRUE)
  }
  isTRUE(tryCatch(progressr::handlers(global = NA), error = function(e) FALSE))
}

evinf_progress_run <- function(n, verbose, fun) {
  if (!evinf_progress_active(verbose)) {
    return(fun(function(...) invisible()))
  }
  progressr::with_progress({
    p <- progressr::progressor(steps = n)
    fun(function(...) p())
  })
}
