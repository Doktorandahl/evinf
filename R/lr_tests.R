



# Refit one restricted model (audit R0.2): the full model's control settings
# and warm starts, so the fits are genuinely nested and the LR statistic cannot
# come out negative because of a different c-grid. Top-level (not a closure) so
# parallel workers do not receive the calling frame.
#   reduced  the $formulas list from formula_var_remover()
#   data     the (possibly resampled) data to fit on
#   obj      a stand-in carrying only $control and $coef, with the model class
#   md       the full model data (for resolving reduced design column names)
lr_refit_restricted <- function(reduced, data, obj, md, verbose = FALSE) {
  if (inherits(obj, "evzinb")) {
    do.call(evinf::evzinb, c(
      list(reduced$nb, reduced$zi, reduced$evinf, reduced$pareto,
           data = data, bootstrap = FALSE, verbose = verbose),
      restricted_fit_args(obj, reduced, md)
    ))
  } else {
    do.call(evinf::evinb, c(
      list(reduced$nb, reduced$evinf, reduced$pareto,
           data = data, bootstrap = FALSE, verbose = verbose),
      restricted_fit_args(obj, reduced, md)
    ))
  }
}

#' Likelihood ratio test for individual variables of evzinb
#'
#' @param object EVZINB or EVINB object to perform likelihood ratio test on
#' @param vars Either a list of character vectors with variable names which to be restricted in the LR test or a character vector of variable names. If a list, each character vector of the list will be run separately, allowing for multiple variables to be restricted as once. If a character vector, parameter 'single' can be used to determine whether all variables in the vector should be restricted at once (single = FALSE) or if the variables should be restricted one by one (single = TRUE)
#' @param single Logical. Determining whether variables in 'vars' should be restricted individually (single = TRUE) or all at once (single = FALSE)
#' @param bootstrap Should LR tests be conducted on each bootstrapped sample or only on the original sample.
#' @inheritParams evzinb
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
lr_test <- function(object, vars, single = TRUE, bootstrap = FALSE, multicore = NULL, ncores = NULL,verbose = FALSE){

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
  lr_obj <- structure(list(control = object$control, coef = object$coef),
                      class = class(object))

  reruns <- purrr::map(formulas_dfs, function(fd)
    lr_refit_restricted(fd$formulas, model_data, lr_obj, model_data, verbose))

    logliks <- reruns %>% purrr::map('log.lik') %>% purrr::reduce(c)
    dfs <- formulas_dfs %>% purrr::map('df') %>% purrr::reduce(c)

    res_full <- tibble::tibble(vars = vars,
                          loglik_full = object$log.lik,
                          loglik_restricted = logliks,
                          df = dfs) %>%
      dplyr::mutate(statistic = 2 * (.data$loglik_full - .data$loglik_restricted),
             prob = pchisq(.data$statistic,.data$df, lower.tail = FALSE))

    if(bootstrap){

    n_f <- length(formulas_dfs)
    n_b <- length(object$bootstraps)
    grid <- expand.grid(fi = seq_len(n_f), bj = seq_len(n_b))
    boot_ids <- purrr::map(object$bootstraps, "boot_id")

    flat <- evinf_with_plan(multicore, ncores, {
      evinf_pmap(
        seq_len(nrow(grid)),
        function(k, grid, formulas_dfs, boot_ids, model_data, lr_obj) {
          fi <- grid$fi[k]; bj <- grid$bj[k]
          try(lr_refit_restricted(formulas_dfs[[fi]]$formulas,
                                  model_data[boot_ids[[bj]], ],
                                  lr_obj, model_data, FALSE))
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

   logliks_boot <- object$bootstraps %>% purrr::map('log.lik') %>% purrr::reduce(c)

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
       sum(is.na(s)), integer(1)))

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
    tl <- attr(terms.formula(f), 'term.labels')
    if(length(tl) == 0){
      return(f)
    }
    keep <- vapply(tl, function(term){
      term_vars <- all.vars(stats::as.formula(paste('~', term)))
      !any(vars %in% term_vars)
    }, logical(1))
    lhs <- deparse(f[[2]])
    rhs <- if(any(keep)) paste(tl[keep], collapse = ' + ') else '1'
    stats::as.formula(paste(lhs, '~', rhs), env = environment(f))
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


#' Control + warm-start arguments for the restricted LR-test refits (audit R0.2)
#'
#' @param object The fitted evzinb/evinb object.
#' @param reduced The \code{$formulas} list from \code{formula_var_remover()}.
#' @param data The full model data (for resolving the reduced design column names).
#' @return A named list of arguments to pass to \code{evzinb()} / \code{evinb()}.
#' @noRd
restricted_fit_args <- function(object, reduced, data) {
  ctrl <- object$control
  is_zinb <- inherits(object, 'evzinb')

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

  args <- list(
    max.diff.par = ctrl$max.diff.par,
    max.no.em.steps = ctrl$max.no.em.steps,
    max.no.em.steps.warmup = ctrl$max.no.em.steps.warmup,
    c.lim = ctrl$c.lim,
    prune.c.range = ctrl$prune.c.range,
    max.upd.par.pl.multinomial = ctrl$max.upd.par.pl.multinomial,
    max.upd.par.nb = ctrl$max.upd.par.nb,
    max.upd.par.pl = ctrl$max.upd.par.pl,
    no.m.bfgs.steps.multinomial = ctrl$no.m.bfgs.steps.multinomial,
    no.m.bfgs.steps.nb = ctrl$no.m.bfgs.steps.nb,
    no.m.bfgs.steps.pl = ctrl$no.m.bfgs.steps.pl,
    pdf.pl.type = ctrl$pdf.pl.type,
    eta.int = ctrl$eta.int,
    init.Alpha.NB = as.numeric(object$coef$Alpha.NB),
    init.C = as.numeric(object$coef$C),
    init.Beta.NB = beta_start(object$coef$Beta.NB, reduced$nb),
    init.Beta.multinom.PL = beta_start(object$coef$Beta.multinom.PL, reduced$evinf),
    init.Beta.PL = beta_start(object$coef$Beta.PL, reduced$pareto)
  )
  if (is_zinb) {
    args$max.upd.par.zc.multinomial <- ctrl$max.upd.par.zc.multinomial
    args$init.Beta.multinom.ZC <- beta_start(object$coef$Beta.multinom.ZC, reduced$zi)
  }
  args
}
