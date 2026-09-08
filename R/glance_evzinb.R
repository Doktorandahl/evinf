#' EVZINB and EVINB glance functions
#'
#' @param x An EVZINB or EVINB object
#' @param ... Further arguments to be passed to glance()
#'
#' @return A one-row tibble of goodness-of-fit statistics: number of observations
#'   and parameters, alpha_NB, C_EV, AIC, BIC, log-likelihood, whether the EM
#'   algorithm converged, and the number of (failed) bootstraps (\code{NA} when
#'   the model was fitted without bootstrapping).
#' @seealso \code{\link[generics]{glance}}
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' glance(model)
#' }
glance.evzinb <- function(x, ...) {
  boot <- evinf_boot_counts(x)

  tibble::tibble(
    nobs = nrow(x$data$x.nb),
    npar = length(x$par.all),
    alpha = x$coef$Alpha.NB,
    parameter = x$coef$C,
    aic = x$AIC,
    bic = x$BIC,
    logLik = x$log.lik,
    converged = isTRUE(x$converge),
    n_bootstraps = boot$n_bootstraps,
    n_failed_bootstraps = boot$n_failed_bootstraps
  )
}

#' EVZINB and EVINB glance functions
#'
#' @param x An EVZINB or EVINB object
#' @param ... Further arguments to be passed to glance()
#'
#' @return A one-row tibble of goodness-of-fit statistics: number of observations
#'   and parameters, alpha_NB, C_EV, AIC, BIC, log-likelihood, whether the EM
#'   algorithm converged, and the number of (failed) bootstraps (\code{NA} when
#'   the model was fitted without bootstrapping).
#' @seealso \code{\link[generics]{glance}}
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' glance(model)
#' }
glance.evinb <- function(x, ...) {
  boot <- evinf_boot_counts(x)

  tibble::tibble(
    nobs = nrow(x$data$x.nb),
    npar = length(x$par.all),
    alpha = x$coef$Alpha.NB,
    parameter = x$coef$C,
    aic = x$AIC,
    bic = x$BIC,
    logLik = x$log.lik,
    converged = isTRUE(x$converge),
    n_bootstraps = boot$n_bootstraps,
    n_failed_bootstraps = boot$n_failed_bootstraps
  )
}

# Number of successful / failed bootstraps, or NA when not bootstrapped.
evinf_boot_counts <- function(x) {
  if (is.null(x$bootstraps)) {
    return(list(n_bootstraps = NA_integer_, n_failed_bootstraps = NA_integer_))
  }
  n_total <- length(x$bootstraps)
  n_failed <- sum(vapply(x$bootstraps, function(b) inherits(b, 'try-error'), logical(1)))
  list(n_bootstraps = n_total - n_failed, n_failed_bootstraps = n_failed)
}
