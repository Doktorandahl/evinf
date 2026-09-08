#' zinbboot and nbboot glance functions
#'
#' @param x An nbboot or zinbboot object
#' @param ... Further arguments to be passed to glance()
#'
#' @return A one-row tibble of goodness-of-fit statistics, including whether the
#'   full-sample model converged and the number of (failed) bootstraps.
#' @seealso \code{\link[generics]{glance}}
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' zinb_comp <- compare_models(model)
#' glance(zinb_comp$zinb)
#' }
glance.zinbboot <- function(x, ...) {
  boot <- evinf_boot_counts(x)

  tibble::tibble(
    nobs = x$full_run$n,
    npar = nrow(x$full_run$vcov) + 1,
    alpha = 1 / x$full_run$theta,
    aic = AIC(x$full_run),
    bic = BIC(x$full_run),
    logLik = x$full_run$loglik,
    converged = isTRUE(x$full_run$converged),
    n_bootstraps = boot$n_bootstraps,
    n_failed_bootstraps = boot$n_failed_bootstraps
  )
}

#' @rdname glance.zinbboot
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' zinb_comp <- compare_models(model)
#' glance(zinb_comp$nb)
#' }
glance.nbboot <- function(x, ...) {
  boot <- evinf_boot_counts(x)

  tibble::tibble(
    nobs = nrow(x$full_run$model),
    npar = x$full_run$rank,
    alpha = 1 / x$full_run$theta,
    aic = AIC(x$full_run),
    bic = BIC(x$full_run),
    logLik = as.numeric(logLik(x$full_run)),
    converged = isTRUE(x$full_run$converged),
    n_bootstraps = boot$n_bootstraps,
    n_failed_bootstraps = boot$n_failed_bootstraps
  )
}
