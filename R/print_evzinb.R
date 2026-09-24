#' EVZINB print function
#'
#' @param x A fitted evzinb model
#' @param ... Not used
#' @return \code{x}, invisibly.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' print(model)
#' }
print.evzinb <- function(x, ...) {
  evinf_print_fit(x, "EVZINB", zi = TRUE)
}

# Shared body for print.evzinb() / print.evinb().
evinf_print_fit <- function(x, kind, zi) {
  fmt_formula <- function(f) paste(deparse(f), collapse = " ")
  n_boot <- evinf_boot_counts(x)
  c.lim <- x$control$c.lim
  on_boundary <- !is.null(x$c_profile) &&
    (x$coef$C <= min(x$c_profile$c) || x$coef$C >= max(x$c_profile$c))

  fam <- x$family %||% evinf_family()
  cat('\n', paste0('Fitted ', kind, ' model (', fam$count, ' count, ', fam$zero,
                   ' zero) with formulas:'),
      '\n NB:    ', fmt_formula(x$formulas$formula_nb), sep = ' ')
  if (zi) {
    cat('\n ZI:    ', fmt_formula(x$formulas$formula_zi), sep = ' ')
  }
  cat('\n EVI:   ', fmt_formula(x$formulas$formula_evi),
      '\n Pareto:', fmt_formula(x$formulas$formula_pareto),
      '\n ______',
      '\n Converged:                     ', isTRUE(x$converge),
      '\n C_EV:                          ', x$coef$C,
      '\n Candidate range for C_EV:      ',
      if (is.null(c.lim)) 'default' else paste0('[', c.lim[1], ', ', c.lim[2], ']'),
      if (isTRUE(x$c_lim_default)) ' (data-driven)' else '',
      '\n Observations at or above C_EV: ', sum(x$data$y >= x$coef$C),
      '\n Parameters:                    ', length(x$par.all),
      '\n Bootstraps (failed, degenerate):',
      if (is.na(n_boot$n_bootstraps)) 'none' else
        paste0(n_boot$n_bootstraps, ' (', n_boot$n_failed_bootstraps, ', ',
               n_boot$n_degenerate_bootstraps, ')'),
      sep = ' ')
  if (isTRUE(n_boot$n_failed_bootstraps + n_boot$n_degenerate_bootstraps > 0)) {
    cat('\n   (call failed_bootstraps() for the details)')
  }
  n_c_bnd <- x$n_c_on_boundary
  n_boot_total <- if (is.null(x$bootstraps)) NA_integer_ else length(x$bootstraps)
  if (on_boundary) {
    cat('\n Note: C_EV lies on the boundary of the candidate range')
    if (isTRUE(n_c_bnd > 0)) {
      cat(' (and in ', n_c_bnd, ' of ', n_boot_total,
          ' bootstrap replicates)', sep = '')
    }
    cat('.')
  } else if (isTRUE(n_c_bnd > 0)) {
    cat('\n Note: C_EV reached the boundary of the candidate range in ',
        n_c_bnd, ' of ', n_boot_total, ' bootstrap replicates.', sep = '')
  }
  if (isTRUE(x$loglik_recomputed)) {
    cat('\n Note: log-likelihood / AIC / BIC recomputed from the returned parameters.')
  }
  min_alpha_pl <- if (is.null(x$fitted$alpha.pl)) NA_real_ else min(x$fitted$alpha.pl)
  alpha_pl_floor <- x$control$alpha_pl_floor %||% 0.01
  if (isTRUE(is.finite(min_alpha_pl) && min_alpha_pl < alpha_pl_floor)) {
    cat('\n Note: fitted Pareto shape (alpha_pl) has collapsed for at least one',
        'observation (min = ', signif(min_alpha_pl, 3),
        '); predictions involving it are floored at alpha_pl_floor = ',
        alpha_pl_floor, '. See glance()$min_alpha_pl.', sep = '')
  }
  if (!is.null(x$starts)) {
    # round10 G.1 (audit §5.10)
    st <- evinf_n_starts_at_best(x)
    cat('\n Note: ', st$n_starts_at_best, ' of ', st$n_starts,
        ' starts reached the best log-likelihood (within 1e-4).', sep = '')
  }
  cat('\n')
  invisible(x)
}

#' EVINB print function
#'
#' @param x A fitted evinb model
#' @param ... Not used
#' @return \code{x}, invisibly.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' print(model)
#' }
print.evinb <- function(x, ...) {
  evinf_print_fit(x, "EVINB", zi = FALSE)
}

#' Print method for compare_models() output
#'
#' Shows the compared models, a compact \code{\link{glance}()} table (one row
#' per slot, key fit statistics only), and the paired bootstrap comparison
#' from \code{\link{compare_fit}()} (round10 I.4, audit §5.9) -- previously
#' only the compared-model names and bootstrap count.
#'
#' @param x An \code{evzinbcomp} object returned by \code{\link{compare_models}}.
#' @param metrics Metrics for the \code{compare_fit()} table; see there.
#' @param ... Not used
#' @return \code{x}, invisibly.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' print(compare_models(model))
#' }
print.evzinbcomp <- function(x, metrics = c("aic", "bic", "rmse", "rmsle"), ...) {
  comp_slots <- setdiff(names(x), c('model', 'evzinb'))
  comp_class <- class(x$model)

  cat('\n', 'Model comparison of ', comp_class,
      '\n ', 'Compared models: ', paste(comp_slots, collapse = ', '),
      '\n ', 'Number of compared models: ', length(comp_slots),
      '\n Number of bootstraps:', length(x$model$bootstraps), '\n\n', sep = '')

  gl <- tryCatch(glance.evzinbcomp(x), error = function(e) NULL)
  if (!is.null(gl)) {
    cols <- intersect(c("model", "nobs", "npar", "logLik", "aic", "bic"), names(gl))
    df <- as.data.frame(gl[, cols, drop = FALSE])
    round_cols <- intersect(c("logLik", "aic", "bic"), names(df))
    df[round_cols] <- lapply(df[round_cols], round, digits = 1)
    cat("Fit summary\n")
    print(df, row.names = FALSE)
    cat("\n")
  }

  cf <- tryCatch(compare_fit(x, metrics = metrics), error = function(e) NULL)
  if (!is.null(cf)) {
    print(cf)
  }

  invisible(x)
}


#' Print methods for evzinb / evinb summaries
#'
#' @param x A \code{summary.evzinb} or \code{summary.evinb} object.
#' @param digits Number of significant digits for the coefficient tables.
#' @param signif.stars Logical; show significance stars.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @name print.summary.evzinb
#' @export
print.summary.evzinb <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"), ...) {
  evinf_print_summary(x, 'EVZINB', digits = digits, signif.stars = signif.stars)
}

#' @rdname print.summary.evzinb
#' @export
print.summary.evinb <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"), ...) {
  evinf_print_summary(x, 'EVINB', digits = digits, signif.stars = signif.stars)
}

# Shared formatter for print.summary.evzinb() / print.summary.evinb().
evinf_print_summary <- function(x, model_type, digits, signif.stars) {

  fam <- x$family %||% evinf_family()
  # round9 E.1: the label names whichever count distribution was actually fit.
  count_label <- if (fam$count == "poisson") "Count component (Poisson)" else
    "Count component (negative binomial)"
  comp_labels <- c(count = count_label,
                   zero = 'Zero-inflation component',
                   evi = 'Extreme-value inflation component',
                   pareto = 'Pareto (extreme value) component')

  cat(model_type, 'model summary\n')
  cat(strrep('=', nchar(model_type) + 14), '\n', sep = '')

  n_boot_used <- x$n_bootstraps_used %||% NA_integer_

  for (cn in names(x$coefficients)) {
    tab <- x$coefficients[[cn]]
    cat('\n', comp_labels[[cn]], '\n', sep = '')

    # audit0.10 §1.2: standard_error = FALSE / approx_t_value = FALSE drop the
    # `se` / `approx_t` columns, so cs.ind/tst.ind must be derived from which
    # columns are actually present rather than assumed from ncol(m).
    orig_names <- setdiff(names(tab), 'Variable')
    m <- as.matrix(tab[, orig_names, drop = FALSE])
    rownames(m) <- tab$Variable
    colnames(m) <- vapply(orig_names, function(nm) switch(nm,
      Estimate = 'Estimate', se = 'Std. Error', approx_t = 'approx t',
      bootstrap_p = 'Pr(boot)', approx_p = 'Pr(>|t|)', nm), character(1))

    cs_ind <- which(orig_names %in% c('Estimate', 'se'))
    tst_ind <- which(orig_names == 'approx_t')
    p_col <- which(orig_names %in% c('bootstrap_p', 'approx_p'))
    # audit0.10 §1.11: a bootstrap p-value from n_boot_used draws can never
    # resolve below ~1/n_boot_used; show "< 1/B" (e.g. "<0.01" for 100 usable
    # bootstraps) there instead of the default eps.Pvalue's "<2e-16", which
    # implies a precision the bootstrap never had. The tiny relative nudge
    # makes the floored value itself (exactly 1/B) also print with "<",
    # since format.pval()/printCoefmat() only do that for values strictly
    # below eps.Pvalue.
    eps_pvalue <- if ('bootstrap_p' %in% orig_names && isTRUE(n_boot_used > 0)) {
      (1 / n_boot_used) * (1 + 1e-8)
    } else {
      .Machine$double.eps
    }
    stats::printCoefmat(m, digits = digits,
                        signif.stars = signif.stars && length(p_col) > 0,
                        has.Pvalue = length(p_col) > 0,
                        P.values = length(p_col) > 0,
                        cs.ind = cs_ind,
                        tst.ind = tst_ind,
                        eps.Pvalue = eps_pvalue,
                        na.print = '')
  }

  ms <- x$model_statistics
  cat('\n', strrep('-', 40), '\n', sep = '')
  # round9 E.1: a Poisson count state has no Alpha.NB (Alpha_nb is a length-1
  # NA_real_ vector, not empty, so this indexes safely and prints "NA").
  alpha_nb_est <- if ("Alpha_NB" %in% names(ms$Alpha_nb)) ms$Alpha_nb[['Alpha_NB']] else NA_real_
  cat('alpha_NB: ', signif(alpha_nb_est, digits),
      '   C_EV: ', ms$C[['C']], '\n', sep = '')
  if (!is.null(ms$n_above_c)) {
    cat('Observations at or above C_EV: ', ms$n_above_c, '\n', sep = '')
  }

  props <- x$component_proportions
  cat('Mean state proportions:  ',
      paste(sprintf('%s = %.3f', props$state, props$mean_prop), collapse = '   '),
      '\n', sep = '')

  cat('Observations: ', ms$Obs[['Obs']],
      '   Parameters: ', ms$Obs[['pars']],
      '   df: ', ms$Obs[['df']], '\n', sep = '')
  cat('logLik: ', signif(ms$fit[['logLik']], digits),
      '   AIC: ', signif(ms$fit[['AIC']], digits),
      '   BIC: ', signif(ms$fit[['BIC']], digits),
      '   Converged: ', isTRUE(ms$converged), '\n', sep = '')
  if (isTRUE(ms$loglik_recomputed)) {
    cat('Note: log-likelihood / AIC / BIC recomputed from the returned ',
        'parameters (differed from the EM trace maximum).\n', sep = '')
  }

  n_failed <- x$n_failed_bootstraps
  n_degen <- x$n_degenerate_bootstraps %||% NA_integer_
  cat('Bootstraps: ', if (is.na(n_failed)) 'none (bootstrap = FALSE)' else
      paste0('failed = ', n_failed, ', degenerate = ',
             if (is.na(n_degen)) 0L else n_degen), '\n', sep = '')

  invisible(x)
}
