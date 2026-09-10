


#' Bootstrap coefficient extractor
#'
#' @param object a fitted model with bootstraps of class evzinb, evinb, nbboot, or zinbboot
#' @param ... Arguments passed to methods, in particular \code{component}
#'   (not for nbboot): one of \code{"count"}, \code{"zero"}, \code{"evi"},
#'   \code{"pareto"} or \code{"all"}. The pre-0.9.4 names \code{"nb"}, \code{"zi"}
#'   and \code{"evinf"} are still accepted with a deprecation warning.
#'
#' @return A tibble with coefficient values, one row per bootstrap and component
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' coefficient_extractor(model, component = 'all')
#' }
coefficient_extractor <- function(object,...){
  UseMethod('coefficient_extractor')
}


#' @rdname coefficient_extractor
#' @param component Which component should be extracted
#' @param exclude_degenerate Drop bootstrap replicates flagged degenerate
#'   (default \code{TRUE}); see \code{\link{evinf_control}}.
#' @export
coefficient_extractor.evzinb <- function(object, component = c('all','count','zero','evi','pareto'), exclude_degenerate = TRUE, ...){

  component <- normalize_component(component, c('all','count','zero','evi','pareto'))

  object$bootstraps <- evinf_usable_bootstraps(object, exclude_degenerate)

  count_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows()
  zero_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.ZC') %>% dplyr::bind_rows()
  evi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows()
  pareto_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows()

  if(component == 'count'){
    return(count_boot)
  }else if(component == 'zero'){
    return(zero_boot)
  }else if(component == 'evi'){
    return(evi_boot)
  }else if(component == 'pareto'){
    return(pareto_boot)
  }else{
    dplyr::bind_rows(count_boot %>% dplyr::mutate(.component = 'count'),
                     zero_boot %>% dplyr::mutate(.component = 'zero'),
                     evi_boot %>% dplyr::mutate(.component = 'evi'),
                     pareto_boot %>% dplyr::mutate(.component = 'pareto'))
  }
}

#' @rdname coefficient_extractor
#' @export
coefficient_extractor.evinb <- function(object, component = c('all','count','evi','pareto'), exclude_degenerate = TRUE, ...){

  component <- normalize_component(component, c('all','count','evi','pareto'))

  object$bootstraps <- evinf_usable_bootstraps(object, exclude_degenerate)

  count_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.NB') %>% dplyr::bind_rows()
  evi_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.multinom.PL') %>% dplyr::bind_rows()
  pareto_boot <- object$bootstraps %>% purrr::map('coef') %>% purrr::map('Beta.PL') %>% dplyr::bind_rows()

  if(component == 'count'){
    return(count_boot)
  }else if(component == 'evi'){
    return(evi_boot)
  }else if(component == 'pareto'){
    return(pareto_boot)
  }else{
    dplyr::bind_rows(count_boot %>% dplyr::mutate(.component = 'count'),
                     evi_boot %>% dplyr::mutate(.component = 'evi'),
                     pareto_boot %>% dplyr::mutate(.component = 'pareto'))
  }
}

#' @rdname coefficient_extractor
#' @export
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' zinb_comp <- compare_models(model)
#' coefficient_extractor(zinb_comp$zinb)
#' }
coefficient_extractor.zinbboot <- function(object, component = c('all','count','zero'), ...){

  component <- normalize_component(component, c('all','count','zero'))

  object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))

  count_boot <- object$bootstraps %>% purrr::map('coefficients') %>% purrr::map('count') %>% dplyr::bind_rows()
  zero_boot <- object$bootstraps %>% purrr::map('coefficients') %>% purrr::map('zero') %>% dplyr::bind_rows()

  if(component == 'count'){
    return(count_boot)
  }else if(component == 'zero'){
    return(zero_boot)
  }else{
    dplyr::bind_rows(count_boot %>% dplyr::mutate(.component = 'count'),
                     zero_boot %>% dplyr::mutate(.component = 'zero'))
  }
}

#' @rdname coefficient_extractor
#' @export
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' zinb_comp <- compare_models(model)
#' coefficient_extractor(zinb_comp$nb)
#' }
coefficient_extractor.nbboot <- function(object, ...){

  object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))

  object$bootstraps %>% purrr::map('coefficients') %>% dplyr::bind_rows()
}
