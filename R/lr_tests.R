



#' Likelihood ratio test for individual variables of evzinb
#'
#' @param object EVZINB or EVINB object to perform likelihood ratio test on
#' @param vars Either a list of character vectors with variable names which to be restricted in the LR test or a character vector of variable names. If a list, each character vector of the list will be run separately, allowing for multiple variables to be restricted as once. If a character vector, parameter 'single' can be used to determine whether all variables in the vector should be restricted at once (single = FALSE) or if the variables should be restricted one by one (single = TRUE)
#' @param single Logical. Determining whether variables in 'vars' should be restricted individually (single = TRUE) or all at once (single = FALSE)
#' @param bootstrap Should LR tests be conducted on each bootstrapped sample or only on the original sample.
#' @param multicore Logical. Should the function be run in parallel?
#' @param ncores Number of cores to use if multicore = TRUE
#' @param verbose Logical. Should the function be verbose?
#'
#' @details The likelihood ratio statistic is \eqn{2(\ell_{full} - \ell_{restricted})}
#'   and is compared to a chi-square distribution with degrees of freedom equal to
#'   the number of design-matrix columns removed across all model components (so a
#'   restricted four-level factor contributes three degrees of freedom).
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
lr_test <- function(object, vars, single = TRUE, bootstrap = FALSE, multicore = FALSE, ncores = NULL,verbose = FALSE){
  i <- 'temp_iter'
  j <- 'temp_iter'

  model_data <- object$data$data

  if(!is(vars, 'list')){
    if(single){
      formulas_dfs <- foreach::foreach(i = 1:length(vars)) %do%
        formula_var_remover(object$formulas,vars[i], data = model_data)
    }else{
      formulas_dfs <- foreach::foreach(i = 1:1) %do%
        formula_var_remover(object$formulas,vars, data = model_data)
      vars <- paste(vars,collapse ='_')
    }
  }else{
    formulas_dfs <- foreach::foreach(i = 1:length(vars)) %do%
      formula_var_remover(object$formulas,vars[[i]], data = model_data)
    vars <- vars  %>% purrr::map(~paste(.x,collapse = '_')) %>% purrr::reduce(c)
  }

  if(inherits(object,'evzinb')){
    reruns <- foreach::foreach(i = 1:length(formulas_dfs)) %do%
      evzinb(formulas_dfs[[i]]$formulas$nb,
             formulas_dfs[[i]]$formulas$zi,
             formulas_dfs[[i]]$formulas$evinf,
             formulas_dfs[[i]]$formulas$pareto,
             data = model_data,
             bootstrap = FALSE,
             verbose = verbose)
  }else{
    reruns <- foreach::foreach(i = 1:length(formulas_dfs)) %do%
      evinb(formulas_dfs[[i]]$formulas$nb,
             formulas_dfs[[i]]$formulas$evinf,
             formulas_dfs[[i]]$formulas$pareto,
             data = model_data,
             bootstrap = FALSE,
            verbose = verbose)
  }

    logliks <- reruns %>% purrr::map('log.lik') %>% purrr::reduce(c)
    dfs <- formulas_dfs %>% purrr::map('df') %>% purrr::reduce(c)

    res_full <- tibble::tibble(vars = vars,
                          loglik_full = object$log.lik,
                          loglik_restricted = logliks,
                          df = dfs) %>%
      dplyr::mutate(statistic = 2 * (.data$loglik_full - .data$loglik_restricted),
             prob = pchisq(.data$statistic,.data$df, lower.tail = FALSE))

    if(bootstrap){

    be <- evinf_setup_backend(multicore, ncores)
    on.exit(be$stop(), add = TRUE)
    `%op%` <- be$operator

    if(inherits(object,'evzinb')){
    reruns <- foreach::foreach(i = 1:length(formulas_dfs)) %:%
      foreach::foreach(j = 1:length(object$bootstraps)) %op% {
        if(verbose){
          message(paste('Running bootstrap',j,'for formula',i))
        }
        try(evzinb(formulas_dfs[[i]]$formulas$nb,
                   formulas_dfs[[i]]$formulas$zi,
                   formulas_dfs[[i]]$formulas$evinf,
                   formulas_dfs[[i]]$formulas$pareto,
                   data = model_data[object$bootstraps[[j]]$boot_id,],
                   bootstrap = FALSE))
      }

    }else{
      reruns <- foreach::foreach(i = 1:length(formulas_dfs)) %:%
        foreach::foreach(j = 1:length(object$bootstraps)) %op% {
          if(verbose){
            message(paste('Running bootstrap',j,'for formula',i))
          }
        try(evinb(formulas_dfs[[i]]$formulas$nb,
               formulas_dfs[[i]]$formulas$evinf,
               formulas_dfs[[i]]$formulas$pareto,
               data = model_data[object$bootstraps[[j]]$boot_id,],
               bootstrap = FALSE))
        }
    }


   logliks_boot_reduced <-  foreach::foreach(i = 1:length(formulas_dfs)) %:%
        foreach::foreach(j = 1:length(object$bootstraps),.final = unlist) %do%{
          if('try-error' %in% class(reruns[[i]][[j]])){
            NA
          }else{
            reruns[[i]][[j]]$log.lik
          }
        }

   logliks_boot <- object$bootstraps %>% purrr::map('log.lik') %>% purrr::reduce(c)

   statistics_boot <- foreach::foreach(i = 1:length(formulas_dfs)) %do%
     (2 * (logliks_boot - logliks_boot_reduced[[i]]))


   p_vals <- tibble::tibble(vars = vars,
                    ks_p = foreach::foreach(i = 1:length(statistics_boot),.final = unlist) %do%
                      stats::ks.test(stats::na.omit(statistics_boot[[i]]),pchisq, dfs[i])$p.value,
                    chisq_mean = foreach::foreach(i = 1:length(statistics_boot),.final = unlist) %do%
                      pchisq(mean(stats::na.omit(statistics_boot[[i]])),dfs[i], lower.tail = FALSE),
                    chisq_median = foreach::foreach(i = 1:length(statistics_boot),.final = unlist) %do%
                      pchisq(median(stats::na.omit(statistics_boot[[i]])),dfs[i], lower.tail = FALSE),
                    prop_sig = foreach::foreach(i = 1:length(statistics_boot),.final = unlist) %do%
                      mean(stats::na.omit(statistics_boot[[i]]) > qchisq(0.95,dfs[i])),
                    n_failed_bootstraps = foreach::foreach(i = 1:length(statistics_boot),.final = unlist) %do%
                      sum(is.na(statistics_boot[[i]])))

   res_boot <- foreach::foreach(i = 1:length(statistics_boot)) %do%
     tibble::tibble(ll_reduced = logliks_boot_reduced[[i]],
            ll_full = logliks_boot,
            statistic = statistics_boot[[i]])
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
