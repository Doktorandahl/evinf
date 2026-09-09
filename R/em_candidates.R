# Estimation internals: input preparation for the EM driver.
#
# Two helpers shared by em_fit() (R/em_fit.R) and its sub-routines:
#   em_c_candidates()  - the candidate set for the extreme-value threshold C_EV
#   em_extend_design()  - prepend the intercept column the C++ routines expect
#
# See R/em_fit.R for the map of how the em_*.R files fit together.

#' Candidate set for the extreme-value threshold \eqn{C_{EV}}
#'
#' Builds the set of values over which the EM driver profiles the
#' log-likelihood for the extreme-value threshold \eqn{C_{EV}}: the unique
#' observed response values inside \code{c.lim}, optionally thinned.
#'
#' @param y Numeric response vector.
#' @param c.lim Length-2 numeric giving the lower and upper bound of the search.
#' @param prune.c.range \code{FALSE} for no thinning, or a number in
#'   \code{[0, 1]} giving the proportion of interior candidates to drop (sampled
#'   with probability proportional to the gaps between consecutive candidates, so
#'   the endpoints are always kept).
#'
#' @return A sorted numeric vector of candidate \eqn{C_{EV}} values.
#'
#' @details This is the ECME step of the algorithm in Appendix A1 of Randahl and
#'   Vegelius (2024): \eqn{C_{EV}} is updated by a grid search over the observed
#'   support rather than by a smooth optimiser. A warning is emitted when the set
#'   has more than 100 elements and no pruning was requested, because each extra
#'   candidate is one extra full-data log-likelihood evaluation per EM iteration.
#'
#' @seealso \code{\link{evzinb}()}, \code{\link{evinb}()}
#' @keywords internal
em_c_candidates <- function(y, c.lim, prune.c.range) {
  uy <- unique(sort(y))
  c.range <- uy[uy >= c.lim[1] & uy <= c.lim[2]]

  if (length(c.range) > 100 && isFALSE(prune.c.range)) {
    warning(
      "The c.range contains more than 100 values. If the estimation is slow ",
      "consider reducing the c-range with the c.lim argument of the function, ",
      "or pruning the c-range using the prune.c.range argument."
    )
  }

  if (is.numeric(prune.c.range)) {
    if (prune.c.range < 0 || prune.c.range > 1) {
      stop("The prune.c.range argument must be FALSE or be between 0 and 1")
    }
    gaps <- diff(c.range)
    sample.size <- ceiling(length(c.range) * (1 - prune.c.range)) - 2
    keep.indicies <- c(
      1,
      sample(
        seq(2, length(c.range) - 1),
        size = sample.size,
        prob = gaps[-c(length(c.range) - 1)] / sum(gaps)
      ),
      length(c.range)
    )
    c.range <- c.range[sort(keep.indicies)]
    if (length(c.range) > 100) {
      warning(
        "The pruned c.range still contains more than 100 values. If the ",
        "estimation is slow consider reducing the c-range with the c.lim ",
        "argument of the function, or increasing the prune.c.range argument."
      )
    }
  }

  c.range
}

#' Prepend the intercept column to each component design matrix
#'
#' The C++ estimation routines (\code{log_lik_fun()}, \code{update_bfgs_fun()})
#' expect each component's design matrix to include a leading column of ones.
#' This helper turns the raw \code{x_obj} (design matrices without intercept, as
#' built by \code{evinf_design()}) into the extended matrices, supplying a
#' ones-only matrix for any component with no covariates and a zero offset when
#' none was given.
#'
#' @param x_obj List with elements \code{X.multinom.ZC}, \code{X.multinom.PL},
#'   \code{X.NB}, \code{X.PL} (each a numeric matrix \eqn{n \times p} without an
#'   intercept column, or \code{NULL}) and optionally \code{offset.nb}.
#' @param n Number of observations.
#'
#' @return A list with elements \code{zc}, \code{pl_mult}, \code{nb}, \code{pl}
#'   (the extended numeric matrices) and \code{offset} (length-\code{n} numeric).
#'
#' @seealso \code{\link{evzinb}()}, \code{\link{evinb}()}
#' @keywords internal
em_extend_design <- function(x_obj, n) {
  extend <- function(x) {
    if (is.null(dim(x))) matrix(1, nrow = n, ncol = 1) else cbind(1, x)
  }
  offset <- x_obj$offset.nb
  if (is.null(offset)) offset <- rep(0, n)
  list(
    zc      = extend(x_obj$X.multinom.ZC),
    pl_mult = extend(x_obj$X.multinom.PL),
    nb      = extend(x_obj$X.NB),
    pl      = extend(x_obj$X.PL),
    offset  = offset
  )
}
