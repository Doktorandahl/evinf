#' EVZINB and EVINB glance functions
#'
#' @param x An EVZINB or EVINB object
#' @param ... Further arguments to be passed to glance()
#'
#' @return A one-row tibble of goodness-of-fit statistics: number of observations
#'   and parameters, alpha_NB, C_EV, AIC, BIC, log-likelihood, whether the EM
#'   algorithm converged, the number of EM steps, the number of observations at
#'   or above C_EV, and, for a bootstrapped model, the bootstrap replicate
#'   counts (\code{NA} otherwise): \code{n_bootstraps} (usable),
#'   \code{n_failed_bootstraps}, \code{n_degenerate_bootstraps} --- these three
#'   partition the number of replicates requested (see \code{\link{evinf_control}}
#'   for "degenerate") --- and \code{n_c_on_boundary}, the number of replicates
#'   whose C_EV landed on a candidate-grid endpoint (informational, not counted
#'   as degenerate).
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
    n_above_c = if (is.null(x$n_above_c)) sum(x$data$y >= x$coef$C) else x$n_above_c,
    n_em_steps = if (is.null(x$n_em_steps)) NA_integer_ else x$n_em_steps,
    n_bootstraps = boot$n_bootstraps,
    n_failed_bootstraps = boot$n_failed_bootstraps,
    n_degenerate_bootstraps = boot$n_degenerate_bootstraps,
    n_c_on_boundary = if (is.null(x$n_c_on_boundary)) NA_integer_ else x$n_c_on_boundary
  )
}

#' EVZINB and EVINB glance functions
#'
#' @param x An EVZINB or EVINB object
#' @param ... Further arguments to be passed to glance()
#'
#' @return A one-row tibble of goodness-of-fit statistics: number of observations
#'   and parameters, alpha_NB, C_EV, AIC, BIC, log-likelihood, whether the EM
#'   algorithm converged, the number of EM steps, the number of observations at
#'   or above C_EV, and, for a bootstrapped model, the bootstrap replicate
#'   counts (\code{NA} otherwise): \code{n_bootstraps} (usable),
#'   \code{n_failed_bootstraps}, \code{n_degenerate_bootstraps} --- these three
#'   partition the number of replicates requested (see \code{\link{evinf_control}}
#'   for "degenerate") --- and \code{n_c_on_boundary}, the number of replicates
#'   whose C_EV landed on a candidate-grid endpoint (informational, not counted
#'   as degenerate).
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
    n_above_c = if (is.null(x$n_above_c)) sum(x$data$y >= x$coef$C) else x$n_above_c,
    n_em_steps = if (is.null(x$n_em_steps)) NA_integer_ else x$n_em_steps,
    n_bootstraps = boot$n_bootstraps,
    n_failed_bootstraps = boot$n_failed_bootstraps,
    n_degenerate_bootstraps = boot$n_degenerate_bootstraps,
    n_c_on_boundary = if (is.null(x$n_c_on_boundary)) NA_integer_ else x$n_c_on_boundary
  )
}

# Bootstrap replicate counts, or NA when the model was fitted without
# bootstrapping. `n_bootstraps` counts replicates that neither errored nor came
# out degenerate (see evinf_control(alpha_floor =)).
evinf_boot_counts <- function(x) {
  if (is.null(x$bootstraps)) {
    return(list(n_bootstraps = NA_integer_, n_failed_bootstraps = NA_integer_,
                n_degenerate_bootstraps = NA_integer_))
  }
  b <- x$bootstraps
  n_failed <- sum(vapply(b, function(z) inherits(z, 'try-error'), logical(1)))
  n_degen <- sum(vapply(b, function(z)
    !inherits(z, 'try-error') && isTRUE(z$degenerate), logical(1)))
  list(n_bootstraps = length(b) - n_failed - n_degen,
       n_failed_bootstraps = n_failed,
       n_degenerate_bootstraps = n_degen)
}
