harmonic_calc <- function(pr_count, count, pr_pareto, C, pareto_alpha) {
  pr_count * count + pr_pareto * C * (1 + pareto_alpha) / pareto_alpha
}

# audit R0.5 / 2.7: state-probability columns with the canonical names
# (pr_zero, pr_count, pr_evi) first and the deprecated duplicates
# (pr_zc, pr_pareto) after, matching what predict(type = "states") returns.
canonical_prbs <- function(prbs) {
  if ("pr_zc" %in% names(prbs)) {
    tibble::tibble(
      pr_zero = prbs$pr_zc,
      pr_count = prbs$pr_count,
      pr_evi = prbs$pr_pareto,
      pr_zc = prbs$pr_zc,
      pr_pareto = prbs$pr_pareto
    )
  } else {
    tibble::tibble(
      pr_count = prbs$pr_count,
      pr_evi = prbs$pr_pareto,
      pr_pareto = prbs$pr_pareto
    )
  }
}

explog_calc <- function(pr_count, count, pr_pareto, C, pareto_alpha) {
  pr_count * count + C * pr_pareto * exp(1 / pareto_alpha)
}

# All per-bootstrap prediction quantities for predict.evzinb() / predict.evinb()
# in a single parallel pass (round 5). Returns a list of length length(boots),
# each element list(prbs, cnts, alphs, C, q, harmonic, explog).
evinf_predict_per_boot <- function(boots, newdata, quantile, want_q, evzinb,
                                   multicore, ncores, seed) {
  nd_id <- 1:nrow(newdata)
  evinf_with_plan(multicore, ncores, {
    evinf_pmap(
      boots,
      function(b, newdata, quantile, want_q, evzinb, nd_id) {
        prob_fn  <- if (evzinb) prob_from_evzinb else prob_from_evinb
        qfn      <- if (evzinb) quantiles_from_evzinb else quantiles_from_evinb
        prbs  <- prob_fn(b, newdata = newdata)
        cnts  <- counts_from_evzinb(b, newdata = newdata)
        alphs <- fitted_alpha_from_evzinb(b, newdata = newdata)
        C     <- b$coef$C
        list(
          prbs  = dplyr::bind_cols(prbs, tibble::tibble(id = nd_id)),
          cnts  = dplyr::bind_cols(cnts, tibble::tibble(id = nd_id)),
          alphs = dplyr::bind_cols(alphs, tibble::tibble(id = nd_id)),
          C = C,
          # A single near-degenerate bootstrap fit can make the mixture-quantile
          # solver fail; drop that replicate from the interval rather than
          # aborting the whole prediction (cf. the try-error tolerance elsewhere).
          q = if (want_q) tryCatch(
            tibble::tibble(
              q = qfn(b, quantile, newdata = newdata, return_data = FALSE,
                      multicore = FALSE),
              id = nd_id),
            error = function(e) NULL) else NULL,
          harmonic = tibble::tibble(harmonic = harmonic_calc(
            prbs$pr_count, cnts$count, pr_pareto = prbs$pr_pareto,
            C = C, pareto_alpha = alphs$pareto_alpha), id = nd_id),
          explog = tibble::tibble(explog = explog_calc(
            prbs$pr_count, cnts$count, pr_pareto = prbs$pr_pareto,
            C = C, pareto_alpha = alphs$pareto_alpha), id = nd_id)
        )
      },
      newdata = newdata, quantile = quantile, want_q = want_q, evzinb = evzinb,
      nd_id = nd_id, seed = seed, label = "predict bootstrap"
    )
  })
}


# Table-driven CI / point-estimate spec for the "plain vector" prediction
# types shared by predict.evzinb() / predict.evinb() (audit0.10 F.5): each of
# these needs exactly the same three things -- the per-bootstrap tibble list
# to summarize, which column of it to summarize, and the 1- or 2-column
# output tibble to bind the point estimate onto for confint = TRUE (2 columns
# for 'zi'/'evinf', which keep the deprecated pr_zc/pr_pareto name alongside
# the canonical one). 'states'/'all' (multi-column tibbles) and 'quantile'
# (dynamic column name, and per-replicate solver failures dropped via
# tryCatch in evinf_predict_per_boot()) don't fit this shape and keep their
# own code below. evinb never produces type = 'zi' (normalize_predict_type()'s
# choices exclude it for that class), so this entry is simply unused there.
.evinf_predict_ci_spec <- list(
  harmonic = list(
    boot = "harmonic_boot", boot_col = "harmonic",
    estimate = function(e) e$harmonic,
    out_tibble = function(v) tibble::tibble(harmonic = v)
  ),
  explog = list(
    boot = "explog_boot", boot_col = "explog",
    estimate = function(e) e$explog,
    out_tibble = function(v) tibble::tibble(explog = v)
  ),
  counts = list(
    boot = "cnts_boot", boot_col = "count",
    estimate = function(e) e$cnts$count,
    out_tibble = function(v) tibble::tibble(count = v)
  ),
  pareto_alpha = list(
    boot = "alphs_boot", boot_col = "pareto_alpha",
    estimate = function(e) e$alphs$pareto_alpha,
    out_tibble = function(v) tibble::tibble(pareto_alpha = v)
  ),
  zi = list(
    boot = "prbs_boot", boot_col = "pr_zc", boot_select = "pr_zc",
    estimate = function(e) e$prbs$pr_zc,
    out_tibble = function(v) tibble::tibble(pr_zero = v, pr_zc = v)
  ),
  evinf = list(
    boot = "prbs_boot", boot_col = "pr_pareto", boot_select = "pr_pareto",
    estimate = function(e) e$prbs$pr_pareto,
    out_tibble = function(v) tibble::tibble(pr_evi = v, pr_pareto = v)
  ),
  count_state = list(
    boot = "prbs_boot", boot_col = "pr_count", boot_select = "pr_count",
    estimate = function(e) e$prbs$pr_count,
    out_tibble = function(v) tibble::tibble(pr_count = v)
  )
)

# Shared implementation of predict.evzinb() / predict.evinb() (audit0.10
# F.5): the two methods differed only in which state-probability / quantile
# functions they call and whether 'zi' is a valid type. `evzinb` picks
# between those; everything else (bootstrap-vs-original branching,
# confidence intervals, 'all'/'states' assembly) is one copy.
evinf_predict_engine <- function(object, newdata, type, pred, quantile,
                                 confint, conf_level, multicore, ncores,
                                 return_bootstraps, exclude_degenerate,
                                 evzinb) {
  valid_types <- if (evzinb) {
    c('harmonic', 'explog', 'counts', 'pareto_alpha', 'zi', 'evinf',
      'count_state', 'states', 'all', 'quantile')
  } else {
    c('harmonic', 'explog', 'counts', 'pareto_alpha', 'evinf',
      'count_state', 'states', 'all', 'quantile')
  }
  prob_fn <- if (evzinb) prob_from_evzinb else prob_from_evinb
  qfn <- if (evzinb) quantiles_from_evzinb else quantiles_from_evinb

  pred <- match.arg(pred, c('original', 'bootstrap_median', 'bootstrap_mean'))
  type <- normalize_predict_type(type, valid_types)

  if (type %in% c('states', 'all') & confint) {
    stop('Confidence interval prediction only available for vector outputs')
  }

  if (pred %in% c('bootstrap_median', 'bootstrap_mean') | confint) {
    object$bootstraps <- evinf_usable_bootstraps(object, exclude_degenerate)
    nboots <- length(object$bootstraps)
    if (is.null(newdata)) {
      newdata <- object$data$data
    }
    if (type %in% c('quantile', 'all') && type == 'quantile' && is.null(quantile)) {
      stop('quantile must be provided for quantile prediction')
    }
    per_boot <- evinf_predict_per_boot(
      object$bootstraps, newdata, quantile,
      want_q = type %in% c('quantile', 'all') && !is.null(quantile),
      evzinb = evzinb, multicore = multicore, ncores = ncores,
      seed = object$boot_seeds[[1]] %||% 1L
    )
    prbs_boot     <- purrr::map(per_boot, "prbs")
    cnts_boot     <- purrr::map(per_boot, "cnts")
    alphs_boot    <- purrr::map(per_boot, "alphs")
    C_boot        <- purrr::map_dbl(per_boot, "C")
    q_boot        <- if (type %in% c('quantile', 'all') && !is.null(quantile))
      purrr::map(per_boot, "q") else NULL
    harmonic_boot <- purrr::map(per_boot, "harmonic")
    explog_boot   <- purrr::map(per_boot, "explog")
  }

  if (pred == 'original') {
    if (type %in% c('quantile', 'all')) {
      if (type == 'quantile' & is.null(quantile)) {
        stop('quantile must be provided for quantile prediction')
      } else if (!is.null(quantile)) {
        q <- qfn(
          object,
          quantile,
          newdata = newdata,
          return_data = F,
          multicore = F
        )
      } else {
        q <- NULL
      }
    }
    ## Estimate component probabilities for all individuals
    prbs <- prob_fn(object, newdata = newdata)
    ## Estimate mu_nb for all individuals
    cnts <- counts_from_evzinb(object, newdata = newdata)
    ## Estimate pareto alpha for all individuals
    alphs <- fitted_alpha_from_evzinb(object, newdata = newdata)
    C_est <- object$coef$C

    harmonic <- harmonic_calc(
      pr_count = prbs$pr_count,
      count = cnts$count,
      pr_pareto = prbs$pr_pareto,
      C = C_est,
      pareto_alpha = alphs$pareto_alpha
    )

    explog <- explog_calc(
      pr_count = prbs$pr_count,
      count = cnts$count,
      pr_pareto = prbs$pr_pareto,
      C = C_est,
      pareto_alpha = alphs$pareto_alpha
    )
  } else if (pred == 'bootstrap_median') {
    prbs <- prbs_boot %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(median) %>%
      dplyr::select(-"id")
    cnts <- cnts_boot %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(median) %>%
      dplyr::select(-"id")
    alphs <- alphs_boot %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(median) %>%
      dplyr::select(-"id")

    if (!is.null(q_boot)) {
      q <- q_boot %>%
        dplyr::bind_rows() %>%
        dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(median) %>%
        dplyr::select(-"id") %>%
        dplyr::pull(.data$q)
    } else {
      q <- NULL
    }

    harmonic <- harmonic_boot %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize(harmonic = median(.data$harmonic)) %>%
      dplyr::pull(.data$harmonic)

    explog <- explog_boot %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize(explog = median(.data$explog)) %>%
      dplyr::pull(.data$explog)

    C_est <- median(C_boot)
  } else if (pred == 'bootstrap_mean') {
    warning(
      'Bootstrapped mean predictions are experimental and may yield infinite values'
    )
    prbs <- prbs_boot %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(mean) %>%
      dplyr::select(-"id")
    cnts <- cnts_boot %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(mean) %>%
      dplyr::select(-"id")
    alphs <- alphs_boot %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(mean) %>%
      dplyr::select(-"id")
    if (!is.null(q_boot)) {
      q <- q_boot %>%
        dplyr::bind_rows() %>%
        dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(mean) %>%
        dplyr::select(-"id") %>%
        dplyr::pull(.data$q)
    } else {
      q <- NULL
    }
    C_est <- mean(C_boot)
    harmonic <- harmonic_boot %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize(harmonic = mean(.data$harmonic)) %>%
      dplyr::pull(.data$harmonic)

    explog <- explog_boot %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$id) %>%
      dplyr::summarize(explog = mean(.data$explog)) %>%
      dplyr::pull(.data$explog)
  }

  if (confint) {
    if (type %in% c('states', 'all')) {
      stop(
        "Confidence interval prediction only available for vector predictions (not 'states' or 'all')"
      )
    }
    qs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)

    if (type %in% names(.evinf_predict_ci_spec)) {
      spec <- .evinf_predict_ci_spec[[type]]
      boot_list <- get(spec$boot)
      col <- spec$boot_col
      ci <- boot_list %>%
        dplyr::bind_rows() %>%
        dplyr::group_by(.data$id) %>%
        dplyr::summarize(
          ci_lb = quantile(.data[[col]], qs[1]),
          ci_ub = quantile(.data[[col]], qs[2])
        ) %>%
        dplyr::select(-"id")
      result <- dplyr::bind_cols(spec$out_tibble(spec$estimate(environment())), ci)
      if (return_bootstraps) {
        boots_out <- if (!is.null(spec$boot_select)) {
          purrr::map(boot_list, ~ dplyr::select(.x, spec$boot_select, "id"))
        } else {
          boot_list
        }
        return(list(ci = result, bootstraps = boots_out))
      } else {
        return(result)
      }
    }
    if (type == 'quantile') {
      warning(
        'Confidence interval prediction with Quantiles may yield unstable results'
      )
      ci <- q_boot %>%
        dplyr::bind_rows() %>%
        dplyr::group_by(.data$id) %>%
        dplyr::summarize(
          ci_lb = quantile(.data$q, qs[1]),
          ci_ub = quantile(.data$q, qs[2])
        ) %>%
        dplyr::select(-"id")
      q_name <- paste0('q', 100 * quantile)
      if (return_bootstraps) {
        return(list(
          ci = dplyr::bind_cols(tibble::tibble(!!q_name := q), ci),
          bootstraps = q_boot
        ))
      } else {
        return(dplyr::bind_cols(tibble::tibble(!!q_name := q), ci))
      }
    }
  } else {
    if (type %in% names(.evinf_predict_ci_spec)) {
      return(.evinf_predict_ci_spec[[type]]$estimate(environment()))
    }
    if (type == 'states') {
      return(canonical_prbs(prbs))
    }
    if (type == 'quantile') {
      return(q)
    }
    if (type == 'all') {
      q_name <- paste0('q', 100 * quantile)
      return(dplyr::bind_cols(
        tibble::tibble(
          harmonic = harmonic,
          explog = explog,
          !!q_name := q
        ),
        canonical_prbs(prbs),
        cnts,
        alphs
      ))
    }
  }
}

#' Predictions from evzinb object
#'
#' @details The likelihood, CDF, quantile prediction, residuals and
#'   \code{simulate()} all use the discretised Pareto distribution (the
#'   integer-valued distribution the model is actually fit on). The
#'   \code{'harmonic'} and \code{'explog'} point predictions instead use the
#'   harmonic and geometric means of the *continuous* Pareto distribution as a
#'   closed-form approximation to the corresponding discretised-Pareto moments;
#'   this keeps existing point predictions unchanged but means they are not
#'   computed from exactly the same distribution as the rest of the model.
#'
#' @param object An evzinb object for which to produce predicted values
#' @param newdata Optional new data (tibble) to produce predicted values from
#' @param type Character string, 'harmonic' for the harmonic mean and 'explog' for exponentiated expected log, 'counts' for predicted count of the negative binomial component, 'pareto_alpha' for the predicted pareto alpha value, 'states' for the predicted component states (prior), 'count_state' for predicted probability of the count state, 'evinf' for predicted probability of the pareto state,'zi' for the predicted probability of the zero state, 'all' for all predicted values, and 'quantile' for quantile prediction.
#' @param ... Other arguments passed to predict function
#' @param quantile Quantile for which to produce quantile prediction
#' @inheritParams evzinb
#' @param pred Type of prediction to be used, defaults to the original prediction from the fitted model, with alternatives being the bootstrapped median or mean. Note that bootstrap mean may yield infinite values, especially when doing quantile prediction
#' @param confint Should confidence intervals be made for the predictions? Note: only available for vector type predictions and not 'states' and 'all'.
#' @param conf_level What confidence level should be used for confidence intervals
#' @param return_bootstraps Should the bootstrapped predictions be returned as well? Useful for further custom analyses of the bootstrapped predictions.
#' @param exclude_degenerate Drop bootstrap replicates flagged degenerate (default TRUE); see the alpha_floor argument of evinf_control().
#'
#' @inheritSection evzinb Parallel processing
#' @inheritSection evzinb Reproducibility
#'
#' @return A vector of predicted values for type 'harmonic', 'explog', 'counts', 'pareto_alpha','zi','evinf', 'count_state', and 'quantile' or a tibble of predicted values for type 'states' and 'all' or if confint=T
#'
#' @importFrom rlang :=
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' predict(model)
#' predict(model, type='all', quantile = 0.9) # all available predicted values
#' }
predict.evzinb <- function(
  object,
  newdata = NULL,
  type = c(
    'harmonic',
    'explog',
    'counts',
    'pareto_alpha',
    'zi',
    'evinf',
    'count_state',
    'states',
    'all',
    'quantile'
  ),
  pred = c('original', 'bootstrap_median', 'bootstrap_mean'),
  quantile = NULL,
  confint = FALSE,
  conf_level = 0.95,
  multicore = NULL,
  ncores = NULL,
  return_bootstraps = FALSE,
  exclude_degenerate = TRUE,
  ...
) {
  evinf_predict_engine(
    object, newdata = newdata, type = type, pred = pred, quantile = quantile,
    confint = confint, conf_level = conf_level, multicore = multicore,
    ncores = ncores, return_bootstraps = return_bootstraps,
    exclude_degenerate = exclude_degenerate, evzinb = TRUE
  )
}

#' Predictions from evinb object
#'
#' @details The likelihood, CDF, quantile prediction, residuals and
#'   \code{simulate()} all use the discretised Pareto distribution (the
#'   integer-valued distribution the model is actually fit on). The
#'   \code{'harmonic'} and \code{'explog'} point predictions instead use the
#'   harmonic and geometric means of the *continuous* Pareto distribution as a
#'   closed-form approximation to the corresponding discretised-Pareto moments;
#'   this keeps existing point predictions unchanged but means they are not
#'   computed from exactly the same distribution as the rest of the model.
#'
#' @param object An evinb object for which to produce predicted values
#' @param newdata Optional new data (tibble) to produce predicted values from
#' @param type Character string, 'harmonic' for the harmonic mean and 'explog' for exponentiated expected log, 'counts' for predicted count of the negative binomial component, 'pareto_alpha' for the predicted pareto alpha value, 'states' for the predicted component states (prior), 'count_state' for predicted probability of the count state, 'evinf' for predicted probability of the pareto state, 'all' for all predicted values, and 'quantile' for quantile prediction.
#' @param ... Other arguments passed to predict function
#' @param quantile Quantile for which to produce quantile prediction
#' @inheritParams evzinb
#' @param pred Type of prediction to be used, defaults to the original prediction from the fitted model, with alternatives being the bootstrapped median or mean. Note that bootstrap mean may yield infinite values, especially when doing quantile prediction
#' @param confint Should confidence intervals be made for the predictions? Note: only available for vector type predictions and not 'states' and 'all'.
#' @param conf_level What confidence level should be used for confidence intervals
#' @param return_bootstraps Should the bootstrapped predictions be returned as well? Useful for further custom analyses of the bootstrapped predictions.
#' @param exclude_degenerate Drop bootstrap replicates flagged degenerate (default TRUE); see the alpha_floor argument of evinf_control().
#'
#' @inheritSection evzinb Parallel processing
#' @inheritSection evzinb Reproducibility
#'
#' @return A vector of predicted values for type 'harmonic', 'explog', 'counts', 'pareto_alpha','evinf', 'count_state', and 'quantile' or a tibble of predicted values for type 'states' and 'all' or if confint=T
#' @export
#'
#' @importFrom rlang :=
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' predict(model)
#' predict(model, type='all', quantile = 0.9) # all available predicted values
#' }
predict.evinb <- function(
  object,
  newdata = NULL,
  type = c(
    'harmonic',
    'explog',
    'counts',
    'pareto_alpha',
    'evinf',
    'count_state',
    'states',
    'all',
    'quantile'
  ),
  pred = c('original', 'bootstrap_median', 'bootstrap_mean'),
  quantile = NULL,
  confint = FALSE,
  conf_level = 0.95,
  multicore = NULL,
  ncores = NULL,
  return_bootstraps = FALSE,
  exclude_degenerate = TRUE,
  ...
) {
  evinf_predict_engine(
    object, newdata = newdata, type = type, pred = pred, quantile = quantile,
    confint = confint, conf_level = conf_level, multicore = multicore,
    ncores = ncores, return_bootstraps = return_bootstraps,
    exclude_degenerate = exclude_degenerate, evzinb = FALSE
  )
}


#' Random draws from a fitted evzinb model
#'
#' \strong{Superseded:} kept for backwards compatibility, with no plans for
#' removal, but new code should use \code{\link[=simulate.evzinb]{simulate()}},
#' which returns a tidy data frame (and a reproducible \code{"seed"} attribute)
#' instead of a bare vector/list.
#'
#' @param object A fitted EVZINB object
#' @param newdata Optional newdata
#' @param n_draws Number of random draws to make
#'
#' @return A vector of randomly drawn values from the fitted evzinb if n_draws == 1, or a list of length n_draws with random drawn values if n_draws > 1
#' @seealso \code{\link{simulate.evzinb}}
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3, data=genevzinb2, n_bootstraps = 5)
#' revzinb_fit(model)
#' }
revzinb_fit <- function(object, newdata = NULL, n_draws = 1) {
  ## Estimate component probabilities for all individuals
  prbs <- prob_from_evzinb(object, newdata = newdata)
  ## Estimate mu_nb for all individuals
  cnts <- counts_from_evzinb(object, newdata = newdata) %>% dplyr::pull()
  ## Estimate pareto alpha for all individuals
  alphs <- fitted_alpha_from_evzinb(object, newdata = newdata) %>% dplyr::pull()

  alpha_nb <- object$coef$Alpha.NB

  C_est <- object$coef$C

  if (is.null(newdata)) {
    n <- length(object$data$y)
  } else {
    n <- nrow(newdata)
  }
  out <- purrr::map(seq_len(n_draws), function(draw) {
      pl_draws <- rpareto_disc(n, C_est, alphs)
      count_draws <- rnbinom(n, mu = cnts, size = 1 / alpha_nb)
      state_draw <- runif(n)
      prbs %>%
        dplyr::mutate(
          rdraw = dplyr::case_when(
            state_draw <= .data$pr_zc ~ 0,
            state_draw <= .data$pr_zc + .data$pr_count ~ count_draws,
            T ~ pl_draws
          )
        ) %>%
        dplyr::pull(.data$rdraw)
    })
  if (n_draws == 1) {
    return(out[[1]])
  } else {
    return(out)
  }
}

#' Random draws from a fitted evinb model
#'
#' \strong{Superseded:} kept for backwards compatibility, with no plans for
#' removal, but new code should use \code{\link[=simulate.evinb]{simulate()}},
#' which returns a tidy data frame (and a reproducible \code{"seed"} attribute)
#' instead of a bare vector/list.
#'
#' @param object A fitted EVINB object
#' @param newdata Optional newdata
#' @param n_draws Number of random draws to make
#'
#' @return A vector of randomly drawn values from the fitted evinb if n_draws == 1, or a list of length n_draws with random drawn values if n_draws > 1
#' @seealso \code{\link{simulate.evinb}}
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evinb(y~x1+x2+x3, data=genevzinb2, n_bootstraps = 5)
#' revinb_fit(model)
#' }
#'
revinb_fit <- function(object, newdata = NULL, n_draws = 1) {
  ## Estimate component probabilities for all individuals
  prbs <- prob_from_evinb(object, newdata = newdata)
  ## Estimate mu_nb for all individuals
  cnts <- counts_from_evzinb(object, newdata = newdata) %>% dplyr::pull()
  ## Estimate pareto alpha for all individuals
  alphs <- fitted_alpha_from_evzinb(object, newdata = newdata) %>% dplyr::pull()

  alpha_nb <- object$coef$Alpha.NB

  C_est <- object$coef$C
  if (is.null(newdata)) {
    n <- length(object$data$y)
  } else {
    n <- nrow(newdata)
  }

  out <- purrr::map(seq_len(n_draws), function(draw) {
      pl_draws <- rpareto_disc(n, C_est, alphs)
      count_draws <- rnbinom(n, mu = cnts, size = 1 / alpha_nb)
      state_draw <- runif(n)
      prbs %>%
        dplyr::mutate(
          rdraw = dplyr::case_when(
            state_draw <= .data$pr_count ~ count_draws,
            T ~ pl_draws
          )
        ) %>%
        dplyr::pull(.data$rdraw)
    })

  if (n_draws == 1) {
    return(out[[1]])
  } else {
    return(out)
  }
}
