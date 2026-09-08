#' EVZINB summary function
#'
#' @param object an EVZINB object
#' @param coef Type of coefficients. Original are the coefficient estimates from the non-bootstrapped version of the model. 'bootstrapped_mean' are the mean coefficients across bootstraps, and 'bootstrapped_median' are the median coefficients across bootstraps
#' @param standard_error Should standard errors be computed?
#' @param p_value What type of p_values should be computed? 'bootstrapped' are bootstrapped p_values through confidence interval inversion. 'approx' are p-values based on the t-value produced by dividing the coefficient with the standard error. 'both' returns both.
#' @param bootstrapped_props Type of bootstrapped proportions of component proportions to be returned
#' @param approx_t_value Should approximate t-values be returned
#' @param symmetric_bootstrap_p Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.
#' @param ... Additional arguments passed to the summary function
#'
#' @details When \code{object} was fitted with \code{bootstrap = FALSE}, the
#'   bootstrap-based quantities (standard errors, p-values, approximate t-values
#'   and bootstrapped proportions) are unavailable; the summary then reports the
#'   point estimates only, \code{n_failed_bootstraps} is \code{NA}, and a message
#'   is emitted. Requesting bootstrapped coefficients for such a model is an error.
#'
#' @return An EVZINB summary object
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' summary(model)
#' }
summary.evzinb <- function(object, coef = c('original', 'bootstrapped_mean', 'bootstrapped_median'),
                           standard_error = TRUE, p_value = c('bootstrapped', 'approx', 'both', 'none'),
                           bootstrapped_props = c('none', 'mean', 'median'), approx_t_value = TRUE,
                           symmetric_bootstrap_p = TRUE, ...) {

  coef <- match.arg(coef, c('original', 'bootstrapped_mean', 'bootstrapped_median'))
  p_value <- match.arg(p_value, c('bootstrapped', 'approx', 'both', 'none'))
  bootstrapped_props <- match.arg(bootstrapped_props, c('none', 'mean', 'median'))

  parts <- evinf_summary_components(
    object,
    components = c('count', 'zero', 'evi', 'pareto'),
    coef = coef, standard_error = standard_error, p_value = p_value,
    bootstrapped_props = bootstrapped_props, approx_t_value = approx_t_value,
    symmetric_bootstrap_p = symmetric_bootstrap_p
  )

  res <- list(coefficients = parts$coefficients,
              model_statistics = parts$model_statistics,
              component_proportions = parts$component_proportions,
              n_failed_bootstraps = parts$n_failed_bootstraps)
  class(res) <- 'summary.evzinb'
  res
}

#' EVINB summary function
#'
#' @param object an EVINB object
#' @param coef Type of coefficients. Original are the coefficient estimates from the non-bootstrapped version of the model. 'bootstrapped_mean' are the mean coefficients across bootstraps, and 'bootstrapped_median' are the median coefficients across bootstraps
#' @param standard_error Should standard errors be computed?
#' @param p_value What type of p_values should be computed? 'bootstrapped' are bootstrapped p_values through confidence interval inversion. 'approx' are p-values based on the t-value produced by dividing the coefficient with the standard error. 'both' returns both.
#' @param bootstrapped_props Type of bootstrapped proportions of component proportions to be returned
#' @param approx_t_value Should approximate t-values be returned
#' @param symmetric_bootstrap_p Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.
#' @param ... Additional arguments passed to the summary function
#'
#' @inherit summary.evzinb details
#'
#' @return An EVINB summary object
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' summary(model)
#' }
summary.evinb <- function(object, coef = c('original', 'bootstrapped_mean', 'bootstrapped_median'),
                          standard_error = TRUE, p_value = c('bootstrapped', 'approx', 'both', 'none'),
                          bootstrapped_props = c('none', 'mean', 'median'), approx_t_value = TRUE,
                          symmetric_bootstrap_p = TRUE, ...) {

  coef <- match.arg(coef, c('original', 'bootstrapped_mean', 'bootstrapped_median'))
  p_value <- match.arg(p_value, c('bootstrapped', 'approx', 'both', 'none'))
  bootstrapped_props <- match.arg(bootstrapped_props, c('none', 'mean', 'median'))

  parts <- evinf_summary_components(
    object,
    components = c('count', 'evi', 'pareto'),
    coef = coef, standard_error = standard_error, p_value = p_value,
    bootstrapped_props = bootstrapped_props, approx_t_value = approx_t_value,
    symmetric_bootstrap_p = symmetric_bootstrap_p
  )

  res <- list(coefficients = parts$coefficients,
              model_statistics = parts$model_statistics,
              component_proportions = parts$component_proportions,
              n_failed_bootstraps = parts$n_failed_bootstraps)
  class(res) <- 'summary.evinb'
  res
}


# Shared engine for summary.evzinb() / summary.evinb().
evinf_summary_components <- function(object, components, coef, standard_error, p_value,
                                     bootstrapped_props, approx_t_value, symmetric_bootstrap_p) {

  has_zi <- 'zero' %in% components
  has_boot <- !is.null(object$bootstraps)
  n_failed_bootstraps <- NA_integer_

  if (!has_boot) {
    if (coef != 'original') {
      stop("coef = '", coef, "' requires a model fitted with bootstrap = TRUE.", call. = FALSE)
    }
    if (standard_error || approx_t_value || p_value != 'none' || bootstrapped_props != 'none') {
      message('summary(): bootstrap-based quantities (standard errors, p-values, ',
              't-values, bootstrapped proportions) are unavailable for a model ',
              'fitted with bootstrap = FALSE; returning estimates only.')
    }
    standard_error <- FALSE
    approx_t_value <- FALSE
    p_value <- 'none'
    bootstrapped_props <- 'none'
  }
  if (has_boot && p_value %in% c('approx', 'both')) {
    standard_error <- TRUE
  }

  nobs <- nrow(object$data$x.nb)
  npar <- length(object$par.all)

  prop_names <- if (has_zi) c('zero', 'count', 'evi') else c('count', 'evi')
  props <- object$props %>%
    dplyr::as_tibble(.name_repair = ~prop_names) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = 'state') %>%
    dplyr::group_by(.data$state) %>%
    dplyr::summarize(mean_prop = mean(.data$value))
  if (has_zi) {
    props <- props %>% dplyr::slice(c(3, 1, 2))
  }

  alpha_nb <- c(Alpha_NB = as.numeric(object$coef$Alpha.NB))
  C_est <- c(C = as.numeric(object$coef$C))

  boot_tabs <- list()
  if (has_boot) {
    n_bootstraps_org <- length(object$bootstraps)
    object$bootstraps <- object$bootstraps %>% purrr::discard(~ 'try-error' %in% class(.x))
    n_failed_bootstraps <- n_bootstraps_org - length(object$bootstraps)

    boot_long <- function(slot) {
      object$bootstraps %>%
        purrr::map('coef') %>%
        purrr::map(slot) %>%
        dplyr::bind_rows() %>%
        tidyr::pivot_longer(dplyr::everything(), names_to = 'Variable')
    }
    boot_tabs$count <- boot_long('Beta.NB')
    boot_tabs$evi <- boot_long('Beta.multinom.PL')
    boot_tabs$pareto <- boot_long('Beta.PL')
    if (has_zi) {
      boot_tabs$zero <- boot_long('Beta.multinom.ZC')
    }

    props_boot <- object$bootstraps %>%
      purrr::map('props') %>%
      purrr::map(colMeans) %>%
      purrr::reduce(rbind) %>%
      dplyr::as_tibble(.name_repair = ~prop_names) %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = 'state') %>%
      dplyr::group_by(.data$state) %>%
      dplyr::summarize(bootstrap_mean = mean(.data$value),
                       bootstrap_median = median(.data$value),
                       standard_error = sd(.data$value))
    props_boot <- switch(bootstrapped_props,
      median = dplyr::select(props_boot, .data$state, .data$bootstrap_median, .data$standard_error),
      mean = dplyr::select(props_boot, .data$state, .data$bootstrap_mean, .data$standard_error),
      dplyr::select(props_boot, .data$state, .data$standard_error))
    props <- dplyr::left_join(props, props_boot, by = 'state')

    alpha_nb_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Alpha.NB') %>% purrr::reduce(c)
    alpha_nb <- c(alpha_nb, bootstrap_mean = mean(alpha_nb_boot),
                  bootstrap_median = median(alpha_nb_boot), standard_error = sd(alpha_nb_boot))
    C_est_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('C') %>% purrr::reduce(c)
    C_est <- c(C_est, bootstrap_mean = mean(C_est_boot),
               bootstrap_median = median(C_est_boot), standard_error = sd(C_est_boot))
  }

  build_component <- function(orig) {
    boot <- boot_tabs[[orig$slot]]
    tab <- switch(coef,
      original = dplyr::tibble(Variable = orig$names, Estimate = orig$est),
      bootstrapped_mean = dplyr::tibble(Variable = orig$names) %>%
        dplyr::left_join(boot %>% dplyr::group_by(.data$Variable) %>%
                           dplyr::summarize(Estimate = mean(.data$value)), by = 'Variable'),
      bootstrapped_median = dplyr::tibble(Variable = orig$names) %>%
        dplyr::left_join(boot %>% dplyr::group_by(.data$Variable) %>%
                           dplyr::summarize(Estimate = median(.data$value)), by = 'Variable'))
    if (standard_error) {
      tab <- tab %>% dplyr::left_join(
        boot %>% dplyr::group_by(.data$Variable) %>% dplyr::summarize(se = sd(.data$value)),
        by = 'Variable')
    }
    if (approx_t_value) {
      tab <- tab %>% dplyr::mutate(approx_t = .data$Estimate / .data$se)
    }
    if (p_value %in% c('bootstrapped', 'both')) {
      tab <- tab %>% dplyr::left_join(
        dplyr::left_join(boot, tab, by = 'Variable') %>%
          dplyr::group_by(.data$Variable) %>%
          dplyr::summarize(bootstrap_p = bootstrap_p_value_calculator(
            .data$value, .data$Estimate[1], symmetric = symmetric_bootstrap_p)),
        by = 'Variable')
    }
    if (p_value %in% c('approx', 'both')) {
      tab <- tab %>% dplyr::mutate(
        approx_p = 2 * pt(abs(.data$Estimate / .data$se), df = nobs - npar, lower.tail = FALSE))
    }
    tab
  }

  specs <- list(
    count = list(slot = 'count',
                 names = names(object$coef$Beta.NB), est = object$coef$Beta.NB),
    evi = list(slot = 'evi',
               names = names(object$coef$Beta.multinom.PL), est = object$coef$Beta.multinom.PL),
    pareto = list(slot = 'pareto',
                  names = names(object$coef$Beta.PL), est = object$coef$Beta.PL)
  )
  if (has_zi) {
    specs$zero <- list(slot = 'zero',
                       names = names(object$coef$Beta.multinom.ZC), est = object$coef$Beta.multinom.ZC)
  }

  coefficients <- stats::setNames(lapply(components, function(cn) build_component(specs[[cn]])), components)

  list(
    coefficients = coefficients,
    model_statistics = list(Alpha_nb = alpha_nb, C = C_est,
                            Obs = c(Obs = nobs, pars = npar, df = nobs - npar),
                            fit = c(logLik = object$log.lik, AIC = object$AIC, BIC = object$BIC),
                            converged = isTRUE(object$converge)),
    component_proportions = props,
    n_failed_bootstraps = n_failed_bootstraps
  )
}


error_remover <- function(object){
  if('try-error' %in% class(object)){
    return(NULL)
  }else{
    return(object)
  }
}

bootstrap_p_value_calculator <- function(x,estimate = NULL, estimate_fallback = c('median','mean'), symmetric = TRUE){
  if(is.null(estimate)){
    estimate_fallback <- match.arg(estimate_fallback,c('median','mean'))
    estimate <- do.call(estimate_fallback,list(x=x))
  }
  if(symmetric){
   if(estimate>=0){
     return(min(1,2*mean(x<0)))
   }else{
     return(min(1,2*mean(x>=0)))
   }
  }else{
    return(mean(abs(x-estimate)>=abs(estimate)))
  }
}
