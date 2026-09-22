# Estimation internals: one generalised EM step at a fixed threshold.
# See R/em_fit.R for the map of how the em_*.R files fit together.

# round10 0.7 (review §6): shrink the NB line search's eta interval so it
# never probes a negative alpha_nb. On near-Poisson data (alpha_nb already
# close to 0), the unconstrained eta.int let the Newton step's direction
# push alpha_new negative for some trial eta; log_lik_fun()'s NB branch is
# undefined there (returns NaN/Inf), so stats::optimise() silently
# substituted its own "maximum positive value" and warned on every probe --
# a stream of warnings on affected fits, not a correctness issue in the
# accepted step itself (which never lands outside eta.int).
#
# alpha_new(eta) = alpha_old + eta * d is linear in eta (d is the Newton
# step's alpha component), so the zero-crossing eta_boundary = -alpha_old/d
# is the one point that matters: admissible eta is on the same side of 0 as
# eta_boundary. Shrinking that boundary by `margin` toward 0 keeps every eta
# stats::optimise() actually evaluates (interior to the returned interval,
# never the interval's own endpoints) strictly on the positive-alpha side,
# since alpha_new() is monotonic in eta.
evinf_eta_interval_positive_alpha <- function(eta_int, alpha_old, d, margin = 0.01) {
  if (!is.finite(d) || d == 0 || !is.finite(alpha_old)) {
    return(eta_int)
  }
  eta_boundary <- -alpha_old / d
  if (d > 0) {
    c(max(eta_int[1], eta_boundary * (1 - margin)), eta_int[2])
  } else {
    c(eta_int[1], min(eta_int[2], eta_boundary * (1 - margin)))
  }
}

#' One EM step for the EVZINB / EVINB mixture at fixed \eqn{C_{EV}}
#'
#' Performs a single generalised-EM iteration: the C++ routine
#' \code{update_bfgs_fun()} takes a bounded quasi-Newton step in all free
#' parameters, then each component's step is accepted only if a one-dimensional
#' line search over the damping factor \eqn{\eta} does not decrease the
#' log-likelihood (Appendix A1 of Randahl and Vegelius, 2024). The extreme-value
#' threshold \eqn{C_{EV}} is held fixed.
#'
#' @param y Numeric response vector.
#' @param ext Extended design matrices from \code{\link{em_extend_design}}.
#' @param par List of current parameter values: \code{Beta.multinom.ZC},
#'   \code{Beta.multinom.PL}, \code{Beta.NB}, \code{Alpha.NB}, \code{Beta.PL},
#'   \code{C}.
#' @param control An \code{\link{evinf_control}} object (for
#'   \code{max.upd.par.nb}, \code{no.m.bfgs.steps.nb}, \code{eta.int}).
#' @param fixed_zc When \code{TRUE} (the EVINB case) the zero-inflation
#'   multinomial block is not updated and its line search is skipped.
#' @param family An \code{\link{evinf_family}()} object (round9 E.1); the
#'   default reproduces today's NB/mixture model exactly. For
#'   \code{count = "poisson"}, \code{Alpha.NB} never moves (the M-step's
#'   Newton system is block-diagonal with an exactly-zero step in that
#'   coordinate) -- it is simply dropped from \code{par.all} downstream.
#'
#' @return A list with
#'   \describe{
#'     \item{par}{the updated parameter list (same names as the input);}
#'     \item{par_end}{the updated free parameters as a flat numeric vector;}
#'     \item{max_abs_par_diff}{the largest absolute parameter change this step;}
#'     \item{log_lik}{the log-likelihood at the updated parameters;}
#'     \item{upd_obj}{the raw return of \code{update_bfgs_fun()} (its
#'       \code{prop} and \code{resp} elements are the E-step quantities).}
#'   }
#'
#' @seealso \code{\link{em_fit_fixed_c}}, \code{\link{evzinb}()}, \code{\link{evinb}()}
#' @keywords internal
em_step <- function(y, ext, par, control, fixed_zc = FALSE, family = evinf_family()) {
  n_beta_nb <- ncol(ext$nb)
  family_count_code <- evinf_family_count_code(family)
  family_zero_code <- evinf_family_zero_code(family)

  zc_old      <- par$Beta.multinom.ZC
  pl_mult_old <- par$Beta.multinom.PL
  nb_old      <- par$Beta.NB
  alpha_old   <- par$Alpha.NB
  pl_old      <- par$Beta.PL
  c_pl        <- par$C

  par_start <- c(zc_old, pl_mult_old, nb_old, abs(alpha_old), pl_old)

  ll <- function(zc, plm, nb, al, pl) {
    log_lik_fun(zc, plm, nb, al, pl, c_pl,
                ext$zc, ext$pl_mult, ext$nb, ext$pl, y, ext$offset,
                ext$offset_zc, ext$offset_pl_mult, ext$weights, family_count_code,
                family_zero_code, ext$has_weights)
  }

  # audit0.10 §1.8: pdf.pl.type = "exact" uses the discretised-Pareto
  # derivatives in the M-step instead of the continuous-Pareto approximation;
  # the likelihood itself (ll() above) is always the discretised pmf either way.
  upd <- update_bfgs_fun(
    zc_old, pl_mult_old, nb_old, alpha_old, pl_old, c_pl,
    ext$zc, ext$pl_mult, ext$nb, ext$pl, y,
    control$max.upd.par.nb, control$no.m.bfgs.steps.nb, ext$offset,
    ext$offset_zc, ext$offset_pl_mult, ext$weights,
    family_count_code, identical(control$pdf.pl.type, "exact"), family_zero_code,
    ext$has_weights
  )

  # --- take the BFGS values where they are finite, otherwise keep the old ----
  nb_new      <- if (sum(is.na(upd$beta_nb_old))  == 0) upd$beta_nb_old  else nb_old
  alpha_new   <- if (sum(is.na(upd$alpha_nb_old)) == 0) upd$alpha_nb_old else alpha_old
  pl_new      <- if (sum(is.na(upd$beta_pl_old))  == 0) upd$beta_pl_old  else pl_old
  pl_mult_new <- if (sum(is.na(upd$gamma_pl_old)) == 0) upd$gamma_pl_old else pl_mult_old
  zc_new <- if (fixed_zc) {
    zc_old
  } else if (sum(is.na(upd$gamma_z_old)) == 0) {
    upd$gamma_z_old
  } else {
    zc_old
  }

  # --- line search over eta for each component that stepped downhill ---------
  if (!is.na(upd$func_val_after_nb) &&
      upd$func_val_after_nb < upd$func_val_before_bfgs &&
      sum(is.na(upd$change_nb_bfgs)) == 0) {
    theta0 <- c(nb_old, alpha_old)
    obj_nb <- function(eta) {
      th <- theta0 + eta * upd$change_nb_bfgs
      -1.0 * ll(zc_old, pl_mult_old, th[seq_len(n_beta_nb)],
                th[n_beta_nb + 1], pl_old)
    }
    eta_int_nb <- evinf_eta_interval_positive_alpha(
      control$eta.int, alpha_old, upd$change_nb_bfgs[n_beta_nb + 1]
    )
    eta_nb <- stats::optimise(obj_nb, interval = eta_int_nb)$minimum
    th <- theta0 + eta_nb * upd$change_nb_bfgs
    nb_new    <- th[seq_len(n_beta_nb)]
    alpha_new <- th[n_beta_nb + 1]
  }

  if (!is.na(upd$func_val_after_pl) &&
      upd$func_val_after_pl < upd$func_val_before_bfgs) {
    obj_pl <- function(eta) {
      -1.0 * ll(zc_old, pl_mult_old, nb_old, alpha_old,
                pl_old + eta * upd$change_pl_bfgs)
    }
    eta_pl <- stats::optimise(obj_pl, interval = control$eta.int)$minimum
    pl_new <- pl_old + eta_pl * upd$change_pl_bfgs
  }

  if (!fixed_zc &&
      !is.na(upd$func_val_after_mult_z) &&
      upd$func_val_after_mult_z < upd$func_val_before_bfgs) {
    obj_zc <- function(eta) {
      -1.0 * ll(zc_old + eta * upd$change_mult_z_bfgs, pl_mult_old,
                nb_old, alpha_old, pl_old)
    }
    eta_zc <- stats::optimise(obj_zc, interval = control$eta.int)$minimum
    zc_new <- zc_old + eta_zc * upd$change_mult_z_bfgs
  }

  if (!is.na(upd$func_val_after_mult_pl) &&
      upd$func_val_after_mult_pl < upd$func_val_before_bfgs) {
    obj_plm <- function(eta) {
      -1.0 * ll(zc_old, pl_mult_old + eta * upd$change_mult_pl_bfgs,
                nb_old, alpha_old, pl_old)
    }
    eta_plm <- stats::optimise(obj_plm, interval = control$eta.int)$minimum
    pl_mult_new <- pl_mult_old + eta_plm * upd$change_mult_pl_bfgs
  }

  par_end <- c(zc_new, pl_mult_new, nb_new, alpha_new, pl_new)

  list(
    par = list(
      Beta.multinom.ZC = zc_new,
      Beta.multinom.PL = pl_mult_new,
      Beta.NB          = nb_new,
      Alpha.NB         = alpha_new,
      Beta.PL          = pl_new,
      C                = c_pl
    ),
    par_end          = par_end,
    max_abs_par_diff = max(abs(par_end - par_start)),
    log_lik          = ll(zc_new, pl_mult_new, nb_new, alpha_new, pl_new),
    upd_obj          = upd
  )
}
