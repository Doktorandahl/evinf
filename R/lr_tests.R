



# Refit one restricted model (audit R0.2): the full model's control settings
# and warm starts, so the fits are genuinely nested and the LR statistic cannot
# come out negative because of a different c-grid. Top-level (not a closure) so
# parallel workers do not receive the calling frame.
#   reduced  the $formulas list from formula_var_remover()
#   data     the (possibly resampled) data to fit on
#   obj      a stand-in carrying only $control, $coef, $family and $weights,
#            with the model class
#   md       the full model data (for resolving reduced design column names)
#   weights  the weights vector aligned to `data`'s rows (NULL for an
#            unweighted fit); round10 0.1
#
# audit0.10 §1.7: passed as control = restricted_control(...) rather than the
# ~20 individual tuning arguments, so nothing reaches resolve_evinf_control()'s
# deprecation warning (which fires on exactly those arguments).
#
# round10 0.1: the restricted refit must reproduce the full model exactly
# except for the dropped variables -- so it needs the same `family` and
# `weights`, not just the same control/warm start. Without these, a weighted
# fit's restricted log-likelihood came from an *unweighted* refit (statistic
# silently wrong, occasionally negative), and a Poisson or hurdle fit's
# restricted refit errored or compared against the wrong nested model.
lr_refit_restricted <- function(reduced, data, obj, md, verbose = FALSE,
                                 weights = NULL) {
  ctrl <- restricted_control(obj, reduced, md)
  fam <- obj$family %||% evinf_family()
  if (inherits(obj, "evzinb")) {
    evinf::evzinb(reduced$nb, reduced$zi, reduced$evinf, reduced$pareto,
                 data = data, control = ctrl, bootstrap = FALSE, verbose = verbose,
                 weights = weights, family = fam)
  } else {
    evinf::evinb(reduced$nb, reduced$evinf, reduced$pareto,
                data = data, control = ctrl, bootstrap = FALSE, verbose = verbose,
                weights = weights, family = fam)
  }
}

#' Likelihood ratio test for individual variables of evzinb
#'
#' @param object EVZINB or EVINB object to perform likelihood ratio test on
#' @param vars Either a list of character vectors with variable names which to be restricted in the LR test or a character vector of variable names. If a list, each character vector of the list will be run separately, allowing for multiple variables to be restricted as once. If a character vector, parameter 'single' can be used to determine whether all variables in the vector should be restricted at once (single = FALSE) or if the variables should be restricted one by one (single = TRUE)
#' @param single Logical. Determining whether variables in 'vars' should be restricted individually (single = TRUE) or all at once (single = FALSE)
#' @param bootstrap Should LR tests be conducted on each bootstrapped sample or only on the original sample.
#' @inheritParams evzinb
#' @param exclude_degenerate Logical. When \code{bootstrap = TRUE}, also drop
#'   bootstrap replicates flagged degenerate (see \code{\link{failed_bootstraps}}),
#'   not just those that errored. Default \code{TRUE}.
#' @param verbose Logical. Should the function be verbose?
#'
#' @details The likelihood ratio statistic is \eqn{2(\ell_{full} - \ell_{restricted})}
#'   and is compared to a chi-square distribution with degrees of freedom equal to
#'   the number of design-matrix columns removed across all model components (so a
#'   restricted four-level factor contributes three degrees of freedom).
#'
#' @inheritSection evzinb Parallel processing
#' @inheritSection evzinb Reproducibility
#'
#' @return A tibble with one row per performed LR test, or, when \code{bootstrap = TRUE},
#'   a list with the summary tibble (\code{results}) and the per-bootstrap statistics
#'   (\code{boot_results}).
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' lr_test(model,'x1')
#' }
lr_test <- function(object, vars, single = TRUE, bootstrap = FALSE, multicore = NULL, ncores = NULL,
                     exclude_degenerate = TRUE, verbose = FALSE){

  model_data <- object$data$data

  if(!is(vars, 'list')){
    if(single){
      formulas_dfs <- purrr::map(seq_along(vars), function(i)
        formula_var_remover(object$formulas, vars[i], data = model_data))
    }else{
      formulas_dfs <- list(
        formula_var_remover(object$formulas, vars, data = model_data))
      vars <- paste(vars,collapse ='_')
    }
  }else{
    formulas_dfs <- purrr::map(vars, function(v)
      formula_var_remover(object$formulas, v, data = model_data))
    vars <- vars  %>% purrr::map(~paste(.x,collapse = '_')) %>% purrr::reduce(c)
  }

  # A lightweight stand-in carrying only what lr_refit_restricted() needs, so
  # the bootstrap workers never receive the (large) $bootstraps list.
  # round10 0.1: $family and $weights travel with it so the refit is the same
  # model as the full fit, minus the restricted variables.
  lr_obj <- structure(list(control = object$control, coef = object$coef,
                           family = object$family, weights = object$weights),
                      class = class(object))

  reruns <- purrr::map(formulas_dfs, function(fd)
    lr_refit_restricted(fd$formulas, model_data, lr_obj, model_data, verbose,
                        weights = lr_obj$weights))

    logliks <- reruns %>% purrr::map('log.lik') %>% purrr::reduce(c)
    dfs <- formulas_dfs %>% purrr::map('df') %>% purrr::reduce(c)

    res_full <- tibble::tibble(vars = vars,
                          loglik_full = object$log.lik,
                          loglik_restricted = logliks,
                          df = dfs) %>%
      dplyr::mutate(statistic = 2 * (.data$loglik_full - .data$loglik_restricted),
             prob = pchisq(.data$statistic,.data$df, lower.tail = FALSE))

    if(bootstrap){

    # Filter to usable replicates first (audit0.10 §1.1): boot_ids and
    # logliks_boot are then both derived from that same list, with vapply(),
    # so they always have the same length and order as each other and as the
    # per-formula reduced-model refits below.
    boots <- evinf_usable_bootstraps(object, exclude_degenerate)
    n_f <- length(formulas_dfs)
    n_b <- length(boots)
    grid <- expand.grid(fi = seq_len(n_f), bj = seq_len(n_b))
    boot_ids <- purrr::map(boots, "boot_id")
    logliks_boot <- vapply(boots, function(b) b$log.lik, numeric(1))

    flat <- evinf_with_plan(multicore, ncores, {
      evinf_pmap(
        seq_len(nrow(grid)),
        function(k, grid, formulas_dfs, boot_ids, model_data, lr_obj) {
          fi <- grid$fi[k]; bj <- grid$bj[k]
          # round10 0.1: subset weights the same way as the resampled rows,
          # so a weighted bootstrap replicate's restricted refit is weighted
          # on exactly the rows it was resampled onto.
          w <- if (is.null(lr_obj$weights)) NULL else lr_obj$weights[boot_ids[[bj]]]
          try(lr_refit_restricted(formulas_dfs[[fi]]$formulas,
                                  model_data[boot_ids[[bj]], ],
                                  lr_obj, model_data, FALSE, weights = w))
        },
        grid = grid, formulas_dfs = formulas_dfs, boot_ids = boot_ids,
        model_data = model_data, lr_obj = lr_obj,
        seed = object$boot_seeds[[1]], label = "LR refit", verbose = verbose
      )
    })
    # split the flat results back into reruns[[fi]][[bj]]
    reruns <- lapply(seq_len(n_f), function(fi)
      flat[grid$fi == fi])

    logliks_boot_reduced <- purrr::map(reruns, function(fr)
      vapply(fr, function(r)
        if (inherits(r, "try-error")) NA_real_ else r$log.lik, numeric(1)))

   statistics_boot <- purrr::map(logliks_boot_reduced, function(llr)
     2 * (logliks_boot - llr))

   p_vals <- tibble::tibble(
     vars = vars,
     ks_p = vapply(seq_along(statistics_boot), function(i)
       stats::ks.test(stats::na.omit(statistics_boot[[i]]), pchisq, dfs[i])$p.value,
       numeric(1)),
     chisq_mean = vapply(seq_along(statistics_boot), function(i)
       pchisq(mean(stats::na.omit(statistics_boot[[i]])), dfs[i], lower.tail = FALSE),
       numeric(1)),
     chisq_median = vapply(seq_along(statistics_boot), function(i)
       pchisq(median(stats::na.omit(statistics_boot[[i]])), dfs[i], lower.tail = FALSE),
       numeric(1)),
     prop_sig = vapply(seq_along(statistics_boot), function(i)
       mean(stats::na.omit(statistics_boot[[i]]) > qchisq(0.95, dfs[i])),
       numeric(1)),
     n_failed_bootstraps = vapply(statistics_boot, function(s)
       sum(is.na(s)), integer(1)),
     # audit0.10 §1.1: how many replicates fed into this row at all, i.e. those
     # that survived evinf_usable_bootstraps() and had a successful reduced-
     # model refit (n_failed_bootstraps counts the latter kind of failure).
     n_bootstraps_used = vapply(statistics_boot, function(s)
       n_b - sum(is.na(s)), integer(1)))

   res_boot <- purrr::map(seq_along(statistics_boot), function(i)
     tibble::tibble(ll_reduced = logliks_boot_reduced[[i]],
            ll_full = logliks_boot,
            statistic = statistics_boot[[i]]))
   names(res_boot) <- vars

   res_full <- res_full %>% dplyr::left_join(p_vals,by = 'vars')

   out <- list(results = res_full,
               boot_results = res_boot)

   return(out)

  }else{
  return(res_full)
  }
}


#' Remove every term mentioning a set of variables from each model component
#'
#' @param formulas The \code{$formulas} list from a fitted evzinb/evinb object.
#' @param vars Character vector of variable names to restrict.
#' @param data The model data, used to count how many design-matrix columns are
#'   dropped (this is the degrees of freedom of the LR test).
#'
#' @return A list with \code{df} (columns dropped across components) and
#'   \code{formulas} (the restricted formulas per component).
#' @noRd
formula_var_remover <- function(formulas, vars, data){

  drop_terms <- function(f){
    if(is.null(f)){
      return(NULL)
    }
    tt <- stats::terms(f)
    tl <- attr(tt, 'term.labels')
    if(length(tl) == 0){
      return(f)
    }
    drop_idx <- which(vapply(tl, function(term){
      term_vars <- all.vars(stats::as.formula(paste('~', term)))
      any(vars %in% term_vars)
    }, logical(1)))
    if(length(drop_idx) == 0){
      return(f)
    }
    # round10 0.1 (review §1): rebuild from the terms object with
    # stats::drop.terms() rather than pasting the surviving term.labels back
    # together. offset() terms never appear in term.labels, so the old
    # paste-based reconstruction silently dropped every offset (and, since
    # dropping to an empty term.labels forced '1' onto the RHS, silently
    # dropped a `- 1` no-intercept specification too). drop.terms() carries
    # both through because it operates on the terms object itself.
    #
    # round12 A1: on R < 4.4, stats::drop.terms() calls reformulate() with a
    # zero-length termlabels when every term is dropped, which errors
    # ("'termlabels' must be a character vector of length at least one").
    # Newer R tolerates it. Since DESCRIPTION claims R >= 4.1.0, build the
    # emptied formula by hand instead of relying on that version-dependent
    # tolerance -- same result (response, offset()s, intercept setting kept),
    # every R version.
    if (length(drop_idx) == length(tl)) {
      variables <- attr(tt, 'variables')
      resp <- if (attr(tt, 'response') == 1) deparse(variables[[2]]) else NULL
      offset_idx <- attr(tt, 'offset')
      offset_terms <- if (is.null(offset_idx)) character(0) else
        vapply(offset_idx, function(i) deparse(variables[[i + 1]]), character(1))
      rhs <- paste(c(if (attr(tt, 'intercept') == 1) '1' else '0', offset_terms),
                  collapse = ' + ')
      return(stats::as.formula(paste(resp, '~', rhs), env = environment(f)))
    }
    new_tt <- stats::drop.terms(tt, dropx = drop_idx, keep.response = TRUE)
    stats::formula(new_tt)
  }

  ncol_mm <- function(f){
    if(is.null(f)){
      return(0L)
    }
    mf <- stats::model.frame(f, data, na.action = stats::na.omit)
    ncol(stats::model.matrix(f, mf))
  }

  new_nb <- drop_terms(formulas$formula_nb)
  new_zi <- drop_terms(formulas$formula_zi)
  new_evi <- drop_terms(formulas$formula_evi)
  new_pareto <- drop_terms(formulas$formula_pareto)

  df <- (ncol_mm(formulas$formula_nb) - ncol_mm(new_nb)) +
    (ncol_mm(formulas$formula_zi) - ncol_mm(new_zi)) +
    (ncol_mm(formulas$formula_evi) - ncol_mm(new_evi)) +
    (ncol_mm(formulas$formula_pareto) - ncol_mm(new_pareto))

  out <- list()
  out$df <- df

  if(is.null(formulas$formula_zi)){
    out$formulas <- list(nb = new_nb,
                         evinf = new_evi,
                         pareto = new_pareto)
  }else{
    out$formulas <- list(nb = new_nb,
                         zi = new_zi,
                         evinf = new_evi,
                         pareto = new_pareto)
  }
  return(out)
}


#' A copy of the full model's control, warm-started at its estimates (audit0.10 §1.7)
#'
#' Returns a modified copy of \code{object$control} with the \code{init.*}
#' fields (and \code{init.C}) set from the full-model estimate, so the
#' restricted refit is genuinely nested (same tuning settings, same C_EV
#' candidate grid) and starts near the constrained optimum. Every other
#' tuning setting -- including ones added after this object was fitted, e.g.
#' \code{max.c.iter} -- is carried over unchanged because it is a copy of the
#' same \code{evinf_control} object, not a hand-picked subset of arguments.
#'
#' @param object The fitted evzinb/evinb object.
#' @param reduced The \code{$formulas} list from \code{formula_var_remover()}.
#' @param data The full model data (for resolving the reduced design column names).
#' @return A modified \code{evinf_control} object.
#' @noRd
restricted_control <- function(object, reduced, data) {
  ctrl <- object$control

  # Retained coefficients start from the full-model estimate; anything the
  # reduced design somehow adds (it never should) starts at zero.
  beta_start <- function(full_named, reduced_formula) {
    if (is.null(reduced_formula)) {
      return(NULL)
    }
    cn <- c('(Intercept)', colnames(evinf_design(reduced_formula, data)$X))
    v <- full_named[cn]
    v[is.na(v)] <- 0
    as.numeric(v)
  }

  # round10 0.1: a Poisson count state has no Alpha.NB (audit's E.1 drop) --
  # object$coef$Alpha.NB is NULL there, so leave ctrl$init.Alpha.NB at its
  # inherited (unused) value instead of overwriting it with numeric(0), which
  # evinf_control() rejects as "must be a single positive number".
  if (!isTRUE((object$family %||% evinf_family())$count == "poisson")) {
    ctrl$init.Alpha.NB <- as.numeric(object$coef$Alpha.NB)
  }
  ctrl$init.C <- as.numeric(object$coef$C)
  ctrl$init.Beta.NB <- beta_start(object$coef$Beta.NB, reduced$nb)
  ctrl$init.Beta.multinom.PL <- beta_start(object$coef$Beta.multinom.PL, reduced$evinf)
  ctrl$init.Beta.PL <- beta_start(object$coef$Beta.PL, reduced$pareto)
  if (inherits(object, 'evzinb')) {
    ctrl$init.Beta.multinom.ZC <- beta_start(object$coef$Beta.multinom.ZC, reduced$zi)
  }
  ctrl
}
