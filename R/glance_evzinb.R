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
#'   how much of the data a replicate leaves out; \code{n_starts}/
#'   \code{n_starts_at_best} (round10 G.1, \code{NA} for \code{n_starts <= 1}),
#'   how many perturbed starts were tried and how many of them reached the
#'   winning log-likelihood (within 1e-4); \code{median_boot_em_steps}/
#'   \code{n_boot_c_capped} (round10 G.3, \code{NA} without bootstraps), the
#'   median EM step count across usable replicates and how many had their
#'   C_EV profile hit \code{max.c.iter} without settling in either phase.
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
  starts <- evinf_n_starts_at_best(x)
  boot_conv <- evinf_boot_convergence_summary(x)

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
    oob_fraction_max = oob$oob_fraction_max,
    n_starts = starts$n_starts,
    n_starts_at_best = starts$n_starts_at_best,
    median_boot_em_steps = boot_conv$median_boot_em_steps,
    n_boot_c_capped = boot_conv$n_boot_c_capped
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
#'   how much of the data a replicate leaves out; \code{n_starts}/
#'   \code{n_starts_at_best} (round10 G.1, \code{NA} for \code{n_starts <= 1}),
#'   how many perturbed starts were tried and how many of them reached the
#'   winning log-likelihood (within 1e-4); \code{median_boot_em_steps}/
#'   \code{n_boot_c_capped} (round10 G.3, \code{NA} without bootstraps), the
#'   median EM step count across usable replicates and how many had their
#'   C_EV profile hit \code{max.c.iter} without settling in either phase.
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
  starts <- evinf_n_starts_at_best(x)
  boot_conv <- evinf_boot_convergence_summary(x)

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
    oob_fraction_max = oob$oob_fraction_max,
    n_starts = starts$n_starts,
    n_starts_at_best = starts$n_starts_at_best,
    median_boot_em_steps = boot_conv$median_boot_em_steps,
    n_boot_c_capped = boot_conv$n_boot_c_capped
  )
}

# Multi-start summary (round10 G.1): how many of object$starts's replicates
# reached the winning log-likelihood, within 1e-4 -- NA (not 1-of-1) for a
# single-start fit, where $starts is NULL.
evinf_n_starts_at_best <- function(x) {
  if (is.null(x$starts)) {
    return(list(n_starts = NA_integer_, n_starts_at_best = NA_integer_))
  }
  best <- max(x$starts$loglik)
  list(n_starts = nrow(x$starts),
       n_starts_at_best = sum(abs(x$starts$loglik - best) < 1e-4))
}

# Per-bootstrap convergence summary across usable replicates (round10 G.3,
# audit §5.10): median_boot_em_steps, and n_boot_c_capped (replicates whose
# C_EV profile hit max.c.iter without settling, in either phase --
# c_warmup_capped or a FALSE c_converged).
evinf_boot_convergence_summary <- function(x) {
  none <- list(median_boot_em_steps = NA_real_, n_boot_c_capped = NA_integer_)
  if (is.null(x$bootstraps)) {
    return(none)
  }
  ok <- evinf_usable_bootstraps(x)
  if (length(ok) == 0) {
    return(none)
  }
  steps <- vapply(ok, function(z) z$n_em_steps %||% NA_integer_, integer(1))
  capped <- vapply(ok, function(z)
    isTRUE(z$c_warmup_capped) || !isTRUE(z$c_converged), logical(1))
  list(median_boot_em_steps = stats::median(steps, na.rm = TRUE),
       n_boot_c_capped = sum(capped))
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
