# Estimation internals: the EM / ECME driver for evzinb() and evinb().
#
# How the R/em_*.R files fit together:
#
#   em_fit()                          <- the driver (this file)
#     em_c_candidates()  ............ R/em_candidates.R  (the C_EV search grid)
#     loop over two phases (warm-up, convergence):
#       em_fit_fixed_c()  .......... R/em_fixed_c.R      (EM iterations at fixed C_EV)
#         em_step()  ............... R/em_step.R         (one EM iteration)
#           em_extend_design()  ... R/em_candidates.R    (add the intercept column)
#       em_profile_c()  ............ R/em_profile_c.R    (ECME update of C_EV)
#     em_fitted_values()  ......... R/em_fitted.R        (fitted quantities)
#
# run_evzinb() / run_evinb() (R/evzinb.R, R/evinb.R) call em_fit() and reshape
# its return list into the fitted model object; bootrun_*() do the same on a
# resample. The list element names below are read by those callers, by
# summary()/glance()/c_profile() and by the test suite, so they must not change.

#' EM / ECME driver for the EVZINB and EVINB models
#'
#' Estimates an EVZINB or EVINB model by the generalised-EM / ECME algorithm of
#' Appendix A1 in Randahl and Vegelius (2024). It alternates between EM
#' iterations at a fixed extreme-value threshold (\code{\link{em_fit_fixed_c}})
#' and a grid update of that threshold (\code{\link{em_profile_c}}), first in a
#' short "warm-up" phase (\code{control$max.no.em.steps.warmup} EM steps) and
#' then in a "convergence" phase that runs to \code{control$max.diff.par}.
#'
#' @param y Numeric response vector.
#' @param x.obj List of raw component design matrices (elements
#'   \code{X.multinom.ZC}, \code{X.multinom.PL}, \code{X.NB}, \code{X.PL}, each a
#'   numeric matrix \eqn{n \times p} without an intercept column or \code{NULL};
#'   optional \code{offset.nb}).
#' @param ini.val List of starting values (\code{Beta.multinom.ZC},
#'   \code{Beta.multinom.PL}, \code{Beta.NB}, \code{Alpha.NB}, \code{Beta.PL},
#'   \code{C}).
#' @param control An \code{\link{evinf_control}} object.
#' @param model \code{"evzinb"} (all three latent states free) or \code{"evinb"}
#'   (the zero-inflation component is switched off: its coefficients are held at
#'   \code{ini.val$Beta.multinom.ZC} in the convergence phase and in every
#'   \eqn{C_{EV}} profile).
#'
#' @return A list with, among others, \code{par.mat} (estimated parameters),
#'   \code{log.lik}, \code{AIC}, \code{BIC}, \code{resp} (posterior state
#'   probabilities), \code{converge}, \code{c_profile}, \code{c_trace},
#'   \code{log.lik.vec.all}, \code{loglik_recomputed}, the fitted-value vectors
#'   (\code{mu.nb.vec}, \code{alpha.pl.vec}, \code{y.hat.pl*}, ...) and the design
#'   matrices. The exact set and names are consumed by \code{run_evzinb()} /
#'   \code{run_evinb()}.
#'
#' @details For \code{model = "evinb"} the zero-inflation multinomial block is
#'   held at its initial value (\code{ini.val$Beta.multinom.ZC}) throughout both
#'   phases and in every \eqn{C_{EV}} profile.
#'
#' @seealso \code{\link{evzinb}()}, \code{\link{evinb}()}
#' @keywords internal
em_fit <- function(y, x.obj, ini.val, control,
                   model = c("evzinb", "evinb")) {
  model <- match.arg(model)
  ext <- em_extend_design(x.obj, length(y))
  n <- length(y)

  c.range <- em_c_candidates(y, control$c.lim, control$prune.c.range)

  control.warmup <- control
  control.warmup$max.no.em.steps <- control$max.no.em.steps.warmup

  prel.val <- ini.val
  log.lik.vec.all <- NULL
  c_trace <- numeric(0)
  log.lik.vec <- NULL
  est.obj <- NULL

  fixed_zc <- model == "evinb"

  # --- warm-up phase -------------------------------------------------------
  cat("Begin warm-up", "\n", sep = "")
  c.abs.diff <- 100
  while (c.abs.diff > 0) {
    est.obj <- em_fit_fixed_c(y, x.obj, prel.val, control.warmup,
                              fixed_zc = fixed_zc)
    log.lik.vec.all <- c(log.lik.vec.all, est.obj$log.lik.vec)
    prel.val <- est.obj$par.mat

    prof <- em_profile_c(y, x.obj, prel.val, c.range)
    log.lik.vec <- prof$profile$loglik
    c.pl.new <- prof$c_hat
    c_trace <- c(c_trace, c.pl.new)

    c.abs.diff <- abs(c.pl.new - prel.val$C)
    prel.val$C <- c.pl.new
    func.val <- max(log.lik.vec)
    log.lik.vec.all <- c(log.lik.vec.all, func.val)
    cat("The new c is ", c.pl.new, ". The function value is ",
        round(func.val, 4), "\n", sep = "")
    prel.val$C <- c.pl.new
  }

  # --- convergence phase --------------------------------------------------
  cat("End warm-up. Run until convergence", "\n", sep = "")
  c.abs.diff <- 100
  while (c.abs.diff > 0) {
    est.obj <- em_fit_fixed_c(y, x.obj, prel.val, control, fixed_zc = fixed_zc)
    log.lik.vec.all <- c(log.lik.vec.all, est.obj$log.lik.vec)
    prel.val <- est.obj$par.mat

    prof <- em_profile_c(y, x.obj, prel.val, c.range)
    log.lik.vec <- prof$profile$loglik
    c.pl.new <- prof$c_hat
    c_trace <- c(c_trace, c.pl.new)

    c.abs.diff <- abs(c.pl.new - prel.val$C)
    prel.val$C <- c.pl.new
    func.val <- max(log.lik.vec.all)
    log.lik.vec.all <- c(log.lik.vec.all, func.val)
    cat("The new c is ", c.pl.new, ". The function value is ",
        round(func.val, 4), "\n", sep = "")
    prel.val$C <- c.pl.new
  }

  final.val <- prel.val
  props.old <- prel.val$Props

  fv <- em_fitted_values(x.obj, prel.val, props.old, c.pl.new, model = model)

  par.all <- c(
    final.val$Beta.multinom.ZC, final.val$Beta.multinom.PL, final.val$Beta.NB,
    final.val$Alpha.NB, final.val$Beta.PL, final.val$C
  )
  n.par <- length(par.all)

  # audit 2.13: the trace maximum is not guaranteed to be the log-likelihood at
  # the returned parameters. Recompute there; if they disagree, trust the value
  # at the returned parameters for log.lik / AIC / BIC.
  ll.at.par <- log_lik_fun(
    final.val$Beta.multinom.ZC, final.val$Beta.multinom.PL, final.val$Beta.NB,
    final.val$Alpha.NB, final.val$Beta.PL, final.val$C,
    ext$zc, ext$pl_mult, ext$nb, ext$pl, y, ext$offset
  )
  loglik_recomputed <- isTRUE(is.finite(ll.at.par) &&
                                abs(ll.at.par - func.val) > 1e-6)
  if (loglik_recomputed) {
    func.val <- ll.at.par
  }

  list(
    control          = control,
    par.mat          = final.val,
    log.lik.vec.all  = log.lik.vec.all,
    c_profile        = data.frame(c = c.range, loglik = log.lik.vec),
    c_trace          = c_trace,
    log.lik          = func.val,
    resp             = est.obj$resp,
    converge         = est.obj$converge,
    ini.val          = ini.val,
    x.nb             = x.obj$X.NB,
    x.pl             = x.obj$X.PL,
    x.multinom.zc    = x.obj$X.multinom.ZC,
    x.multinom.pl    = x.obj$X.multinom.PL,
    median.pl.vec    = fv$median.pl.vec,
    mean.pl.vec      = fv$mean.pl.vec,
    y                = y,
    y.hat.plmedian   = fv$y.hat.plmedian,
    y.hat.plmean     = fv$y.hat.plmean,
    y.hat.plexpElogy = fv$y.hat.plexpElogy,
    mu.nb.vec        = fv$mu.nb.vec,
    alpha.pl.vec     = fv$alpha.pl.vec,
    exp.E.log.y      = fv$exp.E.log.y,
    y.hat.pl.E.inv.y = fv$y.hat.pl.E.inv.y,
    par.all          = par.all,
    BIC              = log(n) * n.par - 2 * func.val,
    AIC              = 2 * n.par - 2 * func.val,
    loglik_recomputed = loglik_recomputed
  )
}
