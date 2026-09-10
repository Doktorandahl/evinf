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

  cat('\n', paste0('Fitted ', kind, ' model with formulas:'),
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
  if (on_boundary) {
    cat('\n Note: C_EV lies on the boundary of the candidate range.')
  }
  if (isTRUE(x$loglik_recomputed)) {
    cat('\n Note: log-likelihood / AIC / BIC recomputed from the returned parameters.')
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
#' @param x An \code{evzinbcomp} object returned by \code{\link{compare_models}}.
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
print.evzinbcomp <- function(x, ...) {
  comp_slots <- setdiff(names(x), c('model', 'evzinb'))
  comp_class <- class(x$model)

  cat('\n', 'Model comparison of ', comp_class,
      '\n ', 'Compared models: ', paste(comp_slots, collapse = ', '),
      '\n ', 'Number of compared models: ', length(comp_slots),
      '\n Number of bootstraps:', length(x$model$bootstraps), '\n')
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

  comp_labels <- c(count = 'Count component (negative binomial)',
                   zero = 'Zero-inflation component',
                   evi = 'Extreme-value inflation component',
                   pareto = 'Pareto (extreme value) component')

  cat(model_type, 'model summary\n')
  cat(strrep('=', nchar(model_type) + 14), '\n', sep = '')

  for (cn in names(x$coefficients)) {
    tab <- x$coefficients[[cn]]
    cat('\n', comp_labels[[cn]], '\n', sep = '')

    m <- as.matrix(tab[, setdiff(names(tab), 'Variable'), drop = FALSE])
    rownames(m) <- tab$Variable
    colnames(m) <- vapply(colnames(m), function(nm) switch(nm,
      Estimate = 'Estimate', se = 'Std. Error', approx_t = 'approx t',
      bootstrap_p = 'Pr(boot)', approx_p = 'Pr(>|t|)', nm), character(1))

    p_col <- which(colnames(m) %in% c('Pr(boot)', 'Pr(>|t|)'))
    stats::printCoefmat(m, digits = digits,
                        signif.stars = signif.stars && length(p_col) > 0,
                        has.Pvalue = length(p_col) > 0,
                        P.values = length(p_col) > 0,
                        cs.ind = if (ncol(m) >= 2) 1:2 else 1,
                        tst.ind = integer(0),
                        na.print = '')
  }

  ms <- x$model_statistics
  cat('\n', strrep('-', 40), '\n', sep = '')
  cat('alpha_NB: ', signif(ms$Alpha_nb[['Alpha_NB']], digits),
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
