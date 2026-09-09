#' EVZINB tidy function
#'
#' @param x An evzinb object
#' @param component Which component should be shown? One of \code{"count"},
#'   \code{"zi"}, \code{"evi"}, \code{"pareto"} or \code{"all"} (the default).
#' @param coef_type Type of coefficients. Original are the coefficient estimates from the non-bootstrapped version of the model. 'bootstrap_mean' are the mean coefficients across bootstraps, and 'bootstrap_median' are the median coefficients across bootstraps
#' @param standard_error Should bootstrapped standard errors be computed?
#' @param p_value What type of p_values should be computed? 'bootstrapped' are bootstrapped p_values through confidence interval inversion. 'approx' are p-values based on the t-value produced by dividing the coefficient with the standard error.
#' @param confint What type of confidence interval should be computed: \code{"none"},
#'   \code{"bootstrapped"} (percentile interval from the bootstrap distribution) or
#'   \code{"approx"} (\code{estimate +/- qt() * std.error}).
#' @param conf_level Confidence level for the confidence interval
#' @param approx_t_value Should approximate t-values be returned
#' @param symmetric_bootstrap_p Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.
#' @param ... Other arguments passed to the tidy function
#'
#' @details When \code{x} was fitted with \code{bootstrap = FALSE}, standard
#'   errors, p-values, t-values and confidence intervals are unavailable; the
#'   estimate column is returned on its own and a message is emitted. Requesting
#'   bootstrapped coefficients for such a model is an error.
#'
#' @return A tibble with one row per coefficient
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' tidy(model)
#' }
tidy.evzinb <- function(x,
                        component = c('all', 'count', 'zero', 'evi', 'pareto'),
                        coef_type = c('original', 'bootstrap_mean', 'bootstrap_median'),
                        standard_error = TRUE,
                        p_value = c('bootstrapped', 'approx', 'none'),
                        confint = c('none', 'bootstrapped', 'approx'),
                        conf_level = 0.95,
                        approx_t_value = TRUE,
                        symmetric_bootstrap_p = TRUE,
                        ...) {
  coef_type <- match.arg(coef_type, c('original', 'bootstrap_mean', 'bootstrap_median'))
  p_value <- match.arg(p_value, c('bootstrapped', 'approx', 'none'))
  confint <- match.arg(confint, c('none', 'bootstrapped', 'approx'))
  component <- normalize_component(component, c('all', 'count', 'zero', 'evi', 'pareto'))

  evinf_tidy_engine(
    x,
    component = component,
    y_levels = c('zero', 'evi', 'count', 'pareto'),
    coef_type = coef_type, standard_error = standard_error, p_value = p_value,
    confint = confint, conf_level = conf_level, approx_t_value = approx_t_value,
    symmetric_bootstrap_p = symmetric_bootstrap_p
  )
}

#' EVINB tidy function
#'
#' @param x An evinb object
#' @param component Which component should be shown? One of \code{"count"},
#'   \code{"evi"}, \code{"pareto"} or \code{"all"} (the default).
#' @param coef_type Type of coefficients. Original are the coefficient estimates from the non-bootstrapped version of the model. 'bootstrap_mean' are the mean coefficients across bootstraps, and 'bootstrap_median' are the median coefficients across bootstraps
#' @param standard_error Should bootstrapped standard errors be computed?
#' @param p_value What type of p_values should be computed? 'bootstrapped' are bootstrapped p_values through confidence interval inversion. 'approx' are p-values based on the t-value produced by dividing the coefficient with the standard error.
#' @param confint What type of confidence interval should be computed: \code{"none"},
#'   \code{"bootstrapped"} (percentile interval from the bootstrap distribution) or
#'   \code{"approx"} (\code{estimate +/- qt() * std.error}).
#' @param conf_level Confidence level for the confidence interval
#' @param approx_t_value Should approximate t-values be returned
#' @param symmetric_bootstrap_p Should bootstrap p-values be computed as symmetric (leaving alpha/2 percent in each tail)? FALSE gives non-symmetric, but narrower, intervals. TRUE corresponds most closely to conventional p-values.
#' @param ... Other arguments passed to the tidy function
#'
#' @inherit tidy.evzinb details
#'
#' @return A tibble with one row per coefficient
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' tidy(model)
#' }
tidy.evinb <- function(x,
                       component = c('all', 'count', 'evi', 'pareto'),
                       coef_type = c('original', 'bootstrap_mean', 'bootstrap_median'),
                       standard_error = TRUE,
                       p_value = c('bootstrapped', 'approx', 'none'),
                       confint = c('none', 'bootstrapped', 'approx'),
                       conf_level = 0.95,
                       approx_t_value = TRUE,
                       symmetric_bootstrap_p = TRUE,
                       ...) {
  coef_type <- match.arg(coef_type, c('original', 'bootstrap_mean', 'bootstrap_median'))
  p_value <- match.arg(p_value, c('bootstrapped', 'approx', 'none'))
  confint <- match.arg(confint, c('none', 'bootstrapped', 'approx'))
  component <- normalize_component(component, c('all', 'count', 'evi', 'pareto'))

  evinf_tidy_engine(
    x,
    component = component,
    y_levels = c('evi', 'count', 'pareto'),
    coef_type = coef_type, standard_error = standard_error, p_value = p_value,
    confint = confint, conf_level = conf_level, approx_t_value = approx_t_value,
    symmetric_bootstrap_p = symmetric_bootstrap_p
  )
}


# Shared engine for tidy.evzinb() / tidy.evinb().
evinf_tidy_engine <- function(x, component, y_levels, coef_type, standard_error,
                              p_value, confint, conf_level, approx_t_value,
                              symmetric_bootstrap_p) {

  has_zi <- 'zero' %in% y_levels
  has_boot <- !is.null(x$bootstraps)

  if (!has_boot) {
    if (coef_type != 'original') {
      stop("coef_type = '", coef_type, "' requires a model fitted with bootstrap = TRUE.", call. = FALSE)
    }
    if (standard_error || approx_t_value || p_value != 'none' || confint != 'none') {
      message('tidy(): standard errors, p-values, t-values and confidence intervals ',
              'are unavailable for a model fitted with bootstrap = FALSE; ',
              'returning estimates only.')
    }
    standard_error <- FALSE
    approx_t_value <- FALSE
    p_value <- 'none'
    confint <- 'none'
  }
  if (has_boot && (p_value == 'approx' || confint == 'approx')) {
    standard_error <- TRUE
  }

  nobs <- nrow(x$data$x.nb)
  npar <- length(x$par.all)
  err_df <- nobs - npar

  boot_tabs <- list()
  if (has_boot) {
    x$bootstraps <- x$bootstraps %>% purrr::discard(~ 'try-error' %in% class(.x))
    boot_long <- function(slot) {
      x$bootstraps %>%
        purrr::map('coef') %>%
        purrr::map(slot) %>%
        dplyr::bind_rows() %>%
        tidyr::pivot_longer(dplyr::everything(), names_to = 'term')
    }
    boot_tabs$count <- boot_long('Beta.NB')
    boot_tabs$evi <- boot_long('Beta.multinom.PL')
    boot_tabs$pareto <- boot_long('Beta.PL')
    if (has_zi) {
      boot_tabs$zero <- boot_long('Beta.multinom.ZC')
    }
  }

  specs <- list(
    count = list(names = names(x$coef$Beta.NB), est = x$coef$Beta.NB),
    evi = list(names = names(x$coef$Beta.multinom.PL), est = x$coef$Beta.multinom.PL),
    pareto = list(names = names(x$coef$Beta.PL), est = x$coef$Beta.PL)
  )
  if (has_zi) {
    specs$zero <- list(names = names(x$coef$Beta.multinom.ZC), est = x$coef$Beta.multinom.ZC)
  }

  build_component <- function(cn) {
    spec <- specs[[cn]]
    boot <- boot_tabs[[cn]]
    tab <- switch(coef_type,
      original = dplyr::tibble(term = spec$names, estimate = as.numeric(spec$est)),
      bootstrap_mean = dplyr::tibble(term = spec$names) %>%
        dplyr::left_join(boot %>% dplyr::group_by(.data$term) %>%
                           dplyr::summarize(estimate = mean(.data$value)), by = 'term'),
      bootstrap_median = dplyr::tibble(term = spec$names) %>%
        dplyr::left_join(boot %>% dplyr::group_by(.data$term) %>%
                           dplyr::summarize(estimate = median(.data$value)), by = 'term'))

    if (standard_error) {
      tab <- tab %>% dplyr::left_join(
        boot %>% dplyr::group_by(.data$term) %>% dplyr::summarize(std.error = sd(.data$value)),
        by = 'term')
    }
    if (approx_t_value) {
      tab <- tab %>% dplyr::mutate(statistic = .data$estimate / .data$std.error)
    }
    if (p_value == 'bootstrapped') {
      tab <- tab %>% dplyr::left_join(
        dplyr::left_join(boot, tab, by = 'term') %>%
          dplyr::group_by(.data$term) %>%
          dplyr::summarize(p.value = bootstrap_p_value_calculator(
            .data$value, .data$estimate[1], symmetric = symmetric_bootstrap_p)),
        by = 'term')
    } else if (p_value == 'approx') {
      tab <- tab %>% dplyr::mutate(
        p.value = 2 * pt(abs(.data$estimate / .data$std.error), df = err_df, lower.tail = FALSE))
    }
    if (confint == 'bootstrapped') {
      qs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
      ci <- boot %>%
        dplyr::group_by(.data$term) %>%
        dplyr::summarize(conf.low = stats::quantile(.data$value, qs[1], names = FALSE),
                         conf.high = stats::quantile(.data$value, qs[2], names = FALSE))
      tab <- tab %>% dplyr::left_join(ci, by = 'term')
    } else if (confint == 'approx') {
      crit <- stats::qt(1 - (1 - conf_level) / 2, df = err_df)
      tab <- tab %>% dplyr::mutate(conf.low = .data$estimate - crit * .data$std.error,
                                   conf.high = .data$estimate + crit * .data$std.error)
    }
    tab
  }

  if (component != 'all') {
    return(build_component(component))
  }

  # `y.level` follows the broom multinomial convention. As with nnet::multinom,
  # modelsummary() needs `shape = term + y.level ~ model` (or a vcov= argument)
  # to lay out the components; there is no way to make bare modelsummary(model)
  # work while keeping one row per (component, term).
  dplyr::bind_rows(lapply(y_levels, function(cn) {
    dplyr::mutate(build_component(cn), y.level = cn, .before = 1)
  })) %>%
    dplyr::mutate(y.level = factor(.data$y.level, levels = y_levels))
}
