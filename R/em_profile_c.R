# Estimation internals: the log-likelihood profile over the C_EV candidate set.
# See R/em_fit.R for the map of how the em_*.R files fit together.

#' Profile the log-likelihood over the extreme-value threshold
#'
#' Evaluates the full-data log-likelihood at each candidate value of the
#' extreme-value threshold \eqn{C_{EV}}, holding all other parameters fixed, and
#' returns the profile together with its maximiser. This is the ECME update of
#' \eqn{C_{EV}} in Appendix A1 of Randahl and Vegelius (2024); it is run once per
#' EM iteration, in both the warm-up and the convergence phase.
#'
#' @param y Numeric response vector.
#' @param x_obj List of component design matrices (see
#'   \code{\link{em_extend_design}}).
#' @param par List of current parameter values with elements
#'   \code{Beta.multinom.ZC}, \code{Beta.multinom.PL}, \code{Beta.NB},
#'   \code{Alpha.NB}, \code{Beta.PL}.
#' @param c_candidates Numeric vector of candidate \eqn{C_{EV}} values, from
#'   \code{\link{em_c_candidates}}.
#'
#' @return A list with
#'   \describe{
#'     \item{profile}{a data frame with columns \code{c} and \code{loglik};}
#'     \item{c_hat}{the candidate(s) attaining the maximum log-likelihood;}
#'     \item{loglik_max}{that maximum.}
#'   }
#'
#' @seealso \code{\link{evzinb}()}, \code{\link{evinb}()}
#' @keywords internal
em_profile_c <- function(y, x_obj, par, c_candidates) {
  ext <- em_extend_design(x_obj, length(y))

  loglik <- vapply(c_candidates, function(cc) {
    log_lik_fun(
      par$Beta.multinom.ZC, par$Beta.multinom.PL, par$Beta.NB, par$Alpha.NB,
      par$Beta.PL, cc,
      ext$zc, ext$pl_mult, ext$nb, ext$pl, y, ext$offset
    )
  }, numeric(1))

  if (is.infinite(max(loglik))) {
    stop(
      "The log-likelihood is infinite for all tested values of c in the ",
      "c-range. Try expanding the c-range with the c.lim argument of the function"
    )
  }

  list(
    profile = data.frame(c = c_candidates, loglik = loglik),
    c_hat = c_candidates[which(loglik == max(loglik))],
    loglik_max = max(loglik)
  )
}
