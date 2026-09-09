# Estimation internals: fitted-value bookkeeping.
# See R/em_fit.R for the map of how the em_*.R files fit together.

#' Fitted quantities from a fitted EVZINB / EVINB parameter set
#'
#' Computes, for every observation, the count-component mean \eqn{\mu_{NB}}, the
#' Pareto shape \eqn{\alpha_{PL}}, several summaries of the Pareto tail (mean,
#' median, \eqn{\exp E[\log y]}, \eqn{1/E[1/y]}) and the corresponding
#' state-probability-weighted predictions \eqn{\hat y}. This is the vectorised
#' equivalent of the per-observation loop that used to sit at the end of the EM
#' driver; it is identical for the EVZINB and EVINB models.
#'
#' @param x_obj List of component design matrices (see
#'   \code{\link{em_extend_design}}); its \code{offset.nb} element (or a zero
#'   offset) enters \eqn{\mu_{NB}}.
#' @param par List of estimated parameters with elements \code{Beta.NB} and
#'   \code{Beta.PL}.
#' @param props Numeric matrix \eqn{n \times 3} of prior state probabilities
#'   (zero, count, extreme-value), as returned in \code{par.mat$Props}.
#' @param c_ev The estimated extreme-value threshold \eqn{C_{EV}}.
#' @param model Either \code{"evzinb"} or \code{"evinb"}; accepted for interface
#'   stability, but the computation is currently the same for both.
#'
#' @return A named list with \code{mu.nb.vec}, \code{alpha.pl.vec},
#'   \code{exp.E.log.y}, \code{E.inv.y}, \code{mean.pl.vec}, \code{median.pl.vec},
#'   \code{y.hat.plmedian}, \code{y.hat.plmean}, \code{y.hat.plexpElogy} and
#'   \code{y.hat.pl.E.inv.y}, each a length-\eqn{n} numeric vector.
#'
#' @seealso \code{\link{evzinb}()}, \code{\link{evinb}()}
#' @keywords internal
em_fitted_values <- function(x_obj, par, props, c_ev,
                             model = c("evzinb", "evinb")) {
  ext <- em_extend_design(x_obj, nrow(props))

  mu.nb.vec <- as.numeric(exp(ext$nb %*% par$Beta.NB + ext$offset))
  alpha.pl.vec <- as.numeric(exp(ext$pl %*% par$Beta.PL))
  exp.E.log.y <- c_ev * exp(1 / alpha.pl.vec)
  E.inv.y <- alpha.pl.vec / (c_ev * (alpha.pl.vec + 1))
  mean.pl.vec <- ifelse(
    alpha.pl.vec > 1,
    alpha.pl.vec * c_ev / (alpha.pl.vec - 1),
    NA_real_
  )
  median.pl.vec <- c_ev * 2^(1 / alpha.pl.vec)

  list(
    mu.nb.vec        = mu.nb.vec,
    alpha.pl.vec     = alpha.pl.vec,
    exp.E.log.y      = exp.E.log.y,
    E.inv.y          = E.inv.y,
    mean.pl.vec      = mean.pl.vec,
    median.pl.vec    = median.pl.vec,
    y.hat.plmedian   = props[, 2] * mu.nb.vec + props[, 3] * median.pl.vec,
    y.hat.plmean     = props[, 2] * mu.nb.vec + props[, 3] * mean.pl.vec,
    y.hat.plexpElogy = props[, 2] * mu.nb.vec + props[, 3] * exp.E.log.y,
    y.hat.pl.E.inv.y = props[, 2] * mu.nb.vec + props[, 3] / E.inv.y
  )
}
