#' EVZINB and EVINB glance functions
#'
#' @param x An EVZINB or EVINB object
#' @param ... Further arguments to be passed to glance()
#'
#' @return A one-row tibble of goodness-of-fit statistics: number of observations
#'   and parameters, family (round9 E.0, e.g. "nbinom/mixture" or "poisson/mixture"),
#'   alpha_NB (NA when the family has none, round9 E.1), C_EV, AIC, BIC,
#'   log-likelihood, whether the EM
#'   algorithm converged (\code{converged}) and whether the C_EV profile settled
#'   within \code{max.c.iter} (\code{c_converged}; \code{NA} for a model fitted
#'   before this field existed), the number of EM steps, the number of
#'   observations at or above C_EV, the smallest fitted Pareto shape
#'   (\code{min_alpha_pl}; see \code{alpha_pl_floor} in
#'   \code{\link{evinf_control}}), \code{sum_weights} (round9 D.2: the sum of
#'   \code{weights =}, equal to \code{nobs} for an unweighted fit -- this,
#'   not \code{nobs}, is what \code{aic}/\code{bic} and the model's degrees of
#'   freedom are computed from), and, for a bootstrapped model, the bootstrap
#'   replicate counts (\code{NA} otherwise): \code{n_bootstraps} (usable),
#'   \code{n_failed_bootstraps}, \code{n_degenerate_bootstraps} --- these three
#'   partition the number of replicates requested (see \code{\link{evinf_control}}
#'   for "degenerate") --- and \code{n_c_on_boundary}, the number of replicates
#'   whose C_EV landed on a candidate-grid endpoint (informational, not counted
#'   as degenerate); \code{oob_fraction_mean}/\code{oob_fraction_min}/
#'   \code{oob_fraction_max} (round10 0.9), the mean/min/max of each usable
#'   replicate's out-of-bag row fraction (\code{NA} without bootstraps) --
#'   most informative for the block schemes, where overlapping blocks change
#'   how much of the data a replicate leaves out.
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
  oob <- evinf_oob_fraction_summary(x)

  tibble::tibble(
    nobs = nrow(x$data$x.nb),
    sum_weights = evinf_nobs(x),
    npar = length(x$par.all),
    family = paste0(x$family$count %||% "nbinom", "/", x$family$zero %||% "mixture"),
    alpha = x$coef$Alpha.NB %||% NA_real_,
    parameter = x$coef$C,
    aic = x$AIC,
    bic = x$BIC,
    logLik = x$log.lik,
    converged = isTRUE(x$converge),
    c_converged = if (is.null(x$c_converged)) NA else isTRUE(x$c_converged),
    n_above_c = if (is.null(x$n_above_c)) sum(x$data$y >= x$coef$C) else x$n_above_c,
    n_em_steps = if (is.null(x$n_em_steps)) NA_integer_ else x$n_em_steps,
    min_alpha_pl = if (is.null(x$fitted$alpha.pl)) NA_real_ else min(x$fitted$alpha.pl),
    n_bootstraps = boot$n_bootstraps,
    n_failed_bootstraps = boot$n_failed_bootstraps,
    n_degenerate_bootstraps = boot$n_degenerate_bootstraps,
    n_c_on_boundary = if (is.null(x$n_c_on_boundary)) NA_integer_ else x$n_c_on_boundary,
    oob_fraction_mean = oob$oob_fraction_mean,
    oob_fraction_min = oob$oob_fraction_min,
    oob_fraction_max = oob$oob_fraction_max
  )
}

#' EVZINB and EVINB glance functions
#'
#' @param x An EVZINB or EVINB object
#' @param ... Further arguments to be passed to glance()
#'
#' @return A one-row tibble of goodness-of-fit statistics: number of observations
#'   and parameters, family (round9 E.0, e.g. "nbinom/mixture" or "poisson/mixture"),
#'   alpha_NB (NA when the family has none, round9 E.1), C_EV, AIC, BIC,
#'   log-likelihood, whether the EM
#'   algorithm converged (\code{converged}) and whether the C_EV profile settled
#'   within \code{max.c.iter} (\code{c_converged}; \code{NA} for a model fitted
#'   before this field existed), the number of EM steps, the number of
#'   observations at or above C_EV, the smallest fitted Pareto shape
#'   (\code{min_alpha_pl}; see \code{alpha_pl_floor} in
#'   \code{\link{evinf_control}}), \code{sum_weights} (round9 D.2: the sum of
#'   \code{weights =}, equal to \code{nobs} for an unweighted fit -- this,
#'   not \code{nobs}, is what \code{aic}/\code{bic} and the model's degrees of
#'   freedom are computed from), and, for a bootstrapped model, the bootstrap
#'   replicate counts (\code{NA} otherwise): \code{n_bootstraps} (usable),
#'   \code{n_failed_bootstraps}, \code{n_degenerate_bootstraps} --- these three
#'   partition the number of replicates requested (see \code{\link{evinf_control}}
#'   for "degenerate") --- and \code{n_c_on_boundary}, the number of replicates
#'   whose C_EV landed on a candidate-grid endpoint (informational, not counted
#'   as degenerate); \code{oob_fraction_mean}/\code{oob_fraction_min}/
#'   \code{oob_fraction_max} (round10 0.9), the mean/min/max of each usable
#'   replicate's out-of-bag row fraction (\code{NA} without bootstraps) --
#'   most informative for the block schemes, where overlapping blocks change
#'   how much of the data a replicate leaves out.
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
  oob <- evinf_oob_fraction_summary(x)

  tibble::tibble(
    nobs = nrow(x$data$x.nb),
    sum_weights = evinf_nobs(x),
    npar = length(x$par.all),
    family = paste0(x$family$count %||% "nbinom", "/", x$family$zero %||% "mixture"),
    alpha = x$coef$Alpha.NB %||% NA_real_,
    parameter = x$coef$C,
    aic = x$AIC,
    bic = x$BIC,
    logLik = x$log.lik,
    converged = isTRUE(x$converge),
    c_converged = if (is.null(x$c_converged)) NA else isTRUE(x$c_converged),
    n_above_c = if (is.null(x$n_above_c)) sum(x$data$y >= x$coef$C) else x$n_above_c,
    n_em_steps = if (is.null(x$n_em_steps)) NA_integer_ else x$n_em_steps,
    min_alpha_pl = if (is.null(x$fitted$alpha.pl)) NA_real_ else min(x$fitted$alpha.pl),
    n_bootstraps = boot$n_bootstraps,
    n_failed_bootstraps = boot$n_failed_bootstraps,
    n_degenerate_bootstraps = boot$n_degenerate_bootstraps,
    n_c_on_boundary = if (is.null(x$n_c_on_boundary)) NA_integer_ else x$n_c_on_boundary,
    oob_fraction_mean = oob$oob_fraction_mean,
    oob_fraction_min = oob$oob_fraction_min,
    oob_fraction_max = oob$oob_fraction_max
  )
}

# Fit-level oob_fraction summary across usable bootstrap replicates (round10
# 0.9, review §7): oob_fraction is stored on each replicate
# (evzinb.R/evinb.R's bootrun_*()), not on the fit -- this rolls it up for
# glance(). It matters most for the block schemes (moving_block/stationary),
# since overlapping blocks change the out-of-bag set's size relative to
# iid/cluster resampling, but is computed for any bootstrapped fit.
evinf_oob_fraction_summary <- function(x) {
  none <- list(oob_fraction_mean = NA_real_, oob_fraction_min = NA_real_,
              oob_fraction_max = NA_real_)
  if (is.null(x$bootstraps)) {
    return(none)
  }
  ok <- Filter(function(z) !inherits(z, "try-error") && !is.null(z$oob_fraction),
              x$bootstraps)
  if (length(ok) == 0) {
    return(none)
  }
  fr <- vapply(ok, function(z) z$oob_fraction, numeric(1))
  list(oob_fraction_mean = mean(fr), oob_fraction_min = min(fr),
       oob_fraction_max = max(fr))
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
