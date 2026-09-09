# Estimation internals: the EM loop at a fixed threshold.
# See R/em_fit.R for the map of how the em_*.R files fit together.

#' EM iterations for the EVZINB / EVINB mixture at a fixed \eqn{C_{EV}}
#'
#' Runs \code{\link{em_step}} repeatedly, starting from \code{ini.val}, until the
#' largest absolute parameter change falls below \code{control$max.diff.par} or
#' \code{control$max.no.em.steps} iterations have been taken. The extreme-value
#' threshold is held at \code{ini.val$C} throughout; the driver
#' (\code{\link{em_fit}}) alternates calls to this routine with the ECME update
#' of \eqn{C_{EV}} (\code{\link{em_profile_c}}).
#'
#' @param y Numeric response vector.
#' @param x_obj List of raw component design matrices (see
#'   \code{\link{em_extend_design}}).
#' @param ini.val List of starting values: \code{Beta.multinom.ZC},
#'   \code{Beta.multinom.PL}, \code{Beta.NB}, \code{Alpha.NB}, \code{Beta.PL},
#'   \code{C}.
#' @param control An \code{\link{evinf_control}} object.
#' @param fixed_zc When \code{TRUE} (EVINB) the zero-inflation multinomial block
#'   is held at \code{ini.val$Beta.multinom.ZC}.
#'
#' @return A list with
#'   \describe{
#'     \item{par.mat}{the estimated parameters plus \code{Props} (the E-step
#'       prior state probabilities) and \code{C};}
#'     \item{par.all}{the free parameters as a flat numeric vector;}
#'     \item{resp}{the E-step responsibilities (posterior state probabilities);}
#'     \item{log.lik.vec}{the log-likelihood after each EM step;}
#'     \item{log.lik}{the final log-likelihood;}
#'     \item{converge}{\code{TRUE} if the loop stopped on the tolerance, not the
#'       iteration cap;}
#'     \item{n_em_steps}{the number of EM steps taken.}
#'   }
#'
#' @seealso \code{\link{em_fit}}, \code{\link{evzinb}()}, \code{\link{evinb}()}
#' @keywords internal
em_fit_fixed_c <- function(y, x_obj, ini.val, control, fixed_zc = FALSE) {
  ext <- em_extend_design(x_obj, length(y))

  par <- list(
    Beta.multinom.ZC = ini.val$Beta.multinom.ZC,
    Beta.multinom.PL = ini.val$Beta.multinom.PL,
    Beta.NB          = ini.val$Beta.NB,
    Alpha.NB         = ini.val$Alpha.NB,
    Beta.PL          = ini.val$Beta.PL,
    C                = ini.val$C
  )

  max_abs_par_diff <- 100
  i_em <- 1L
  log_lik_vec <- c()
  step <- NULL

  while (i_em < control$max.no.em.steps &&
         max_abs_par_diff > control$max.diff.par) {
    step <- em_step(y, ext, par, control, fixed_zc = fixed_zc)
    par <- step$par
    max_abs_par_diff <- step$max_abs_par_diff
    log_lik_vec[i_em] <- step$log_lik

    cat(
      "Iteration ", i_em,
      ": The max abs diff in parameters is ", max_abs_par_diff,
      ". The function value is ", round(step$log_lik, 4), "\n",
      sep = ""
    )
    i_em <- i_em + 1L
  }

  list(
    par.mat = list(
      Props            = step$upd_obj$prop,
      Beta.multinom.ZC = par$Beta.multinom.ZC,
      Beta.multinom.PL = par$Beta.multinom.PL,
      Beta.NB          = par$Beta.NB,
      Alpha.NB         = par$Alpha.NB,
      Beta.PL          = as.numeric(par$Beta.PL),
      C                = par$C
    ),
    par.all     = step$par_end,
    resp        = step$upd_obj$resp,
    log.lik.vec = log_lik_vec,
    log.lik     = step$log_lik,
    converge    = i_em < control$max.no.em.steps,
    n_em_steps  = i_em - 1L
  )
}
