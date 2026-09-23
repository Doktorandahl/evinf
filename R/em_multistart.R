# Estimation internals: multiple random starts (round10 G.1, audit §5.10).
# See R/em_fit.R for the map of how the em_*.R files fit together.

# em_fit(), always with its own console chatter captured -- run_evzinb()/
# run_evinb() already do this for verbose = FALSE; the multi-start path below
# always suppresses it (n_starts interleaved EM traces are unreadable), and
# `verbose` there instead controls evinf_pmap()'s progress bar.
evinf_em_fit_silent <- function(y, x.obj, ini.val, control, model, family) {
  out <- NULL
  utils::capture.output(
    out <- em_fit(y, x.obj, ini.val, control, model = model, family = family)
  )
  out
}

# Perturb one starting point (round10 G.1): each component's coefficients
# jittered by N(0, start_jitter^2); C_EV drawn uniformly from the candidate
# grid instead of jittered, since it must stay one of the discrete candidates.
evinf_jitter_ini_val <- function(ini.val, start_jitter, c_candidates) {
  jit <- function(v) v + stats::rnorm(length(v), mean = 0, sd = start_jitter)
  ini.val$Beta.multinom.ZC <- jit(ini.val$Beta.multinom.ZC)
  ini.val$Beta.multinom.PL <- jit(ini.val$Beta.multinom.PL)
  ini.val$Beta.NB <- jit(ini.val$Beta.NB)
  ini.val$Beta.PL <- jit(ini.val$Beta.PL)
  ini.val$C <- sample(c_candidates, 1L)
  ini.val
}

#' Run em_fit() from one or more starting points and keep the best (round10 G.1)
#'
#' With \code{control$n_starts <= 1} (the default), this is a single,
#' unperturbed call to \code{\link{em_fit}} -- bit-identical to before
#' \code{n_starts} existed. With \code{control$n_starts > 1}, the default
#' start (\code{ini.val} as given) plus \code{control$n_starts - 1} perturbed
#' starts (\code{Beta.*} jittered by \code{N(0, start_jitter^2)}, \code{C}
#' drawn uniformly from the candidate grid) each run through
#' \code{em_fit()}, in parallel via \code{\link{evinf_pmap}} seeded from
#' \code{start_seed}, and the replicate with the highest final
#' log-likelihood is kept.
#'
#' @param y,x.obj,ini.val,control,model,family As for \code{\link{em_fit}}.
#' @param verbose Controls \code{evinf_pmap()}'s progress bar for the
#'   multi-start pass (each start's own EM console output is always
#'   suppressed there, regardless of \code{verbose} -- see
#'   \code{evinf_em_fit_silent()}). For \code{n_starts <= 1}, controls the
#'   single \code{em_fit()} call's own output, exactly as before.
#' @param start_seed Integer seed for the perturbed starts; drawn if
#'   \code{NULL} and \code{n_starts > 1}.
#'
#' @return A list with \code{best} (the winning \code{em_fit()} return
#'   list), \code{starts} (a tibble with one row per start --
#'   \code{start}, \code{loglik}, \code{C}, \code{converged},
#'   \code{n_em_steps}, and list-columns \code{loglik_trace}/\code{c_trace};
#'   \code{NULL} when \code{n_starts <= 1}, to keep single-start objects
#'   small) and \code{start_seed} (\code{NULL} when \code{n_starts <= 1}).
#' @keywords internal
evinf_run_starts <- function(y, x.obj, ini.val, control, model, family,
                             verbose = FALSE, start_seed = NULL) {
  n_starts <- control$n_starts %||% 1L

  if (n_starts <= 1L) {
    best <- if (verbose) {
      em_fit(y, x.obj, ini.val, control, model = model, family = family)
    } else {
      evinf_em_fit_silent(y, x.obj, ini.val, control, model, family)
    }
    return(list(best = best, starts = NULL, start_seed = NULL))
  }

  if (is.null(start_seed)) {
    start_seed <- sample.int(.Machine$integer.max, 1L)
  }
  # audit0.10-style: the candidate grid is recomputed inside em_fit() anyway;
  # suppress its ">100 values" warning here so it fires (at most) once, from
  # the actual fit, not once per start on top of that.
  c_candidates <- suppressWarnings(
    em_c_candidates(y, control$c.lim, control$prune.c.range)
  )

  fits <- evinf_pmap(
    seq_len(n_starts),
    function(i, ini.val, y, x.obj, control, model, family, start_jitter, c_candidates) {
      iv <- if (i == 1L) {
        ini.val
      } else {
        evinf_jitter_ini_val(ini.val, start_jitter, c_candidates)
      }
      fit <- evinf_em_fit_silent(y, x.obj, iv, control, model, family)
      fit$start_ini_C <- iv$C
      fit
    },
    ini.val = ini.val, y = y, x.obj = x.obj, control = control, model = model,
    family = family, start_jitter = control$start_jitter %||% 0.5,
    c_candidates = c_candidates,
    seed = start_seed, label = "start", verbose = verbose,
    chunk_size = control$chunk_size
  )

  logliks <- vapply(fits, function(f) f$log.lik, numeric(1))
  best_idx <- which.max(logliks)

  starts <- tibble::tibble(
    start = seq_len(n_starts),
    loglik = logliks,
    C = vapply(fits, function(f) f$par.mat$C, numeric(1)),
    converged = vapply(fits, function(f) isTRUE(f$converge), logical(1)),
    n_em_steps = vapply(fits, function(f) length(f$log.lik.vec.all), integer(1)),
    loglik_trace = lapply(fits, function(f) f$log.lik.vec.all),
    c_trace = lapply(fits, function(f) f$c_trace)
  )

  list(best = fits[[best_idx]], starts = starts, start_seed = start_seed)
}
