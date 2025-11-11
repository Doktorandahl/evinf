#' Extracting full mixture quantiles from an evzinb object
#'
#' @param object  An evzinb object for which to produce quantiles
#' @param quantile The quantile for which to produce predictions
#' @param newdata  Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#' @param multicore Logical: Should parallel processing be used to obtain results
#' @param ncores Number of cores if multicore is used
#'
#' @return A vector of predicted quantiles, or if return_data=T, a tibble with the predicted quantile attached last
#'
#' @noRd
quantiles_from_evzinb <- function(
  object,
  quantile,
  newdata = NULL,
  return_data = FALSE,
  multicore = FALSE,
  ncores = NULL
) {
  i <- 'temp_iter'

  ## Estimate component probabilities for all individuals
  prbs <- prob_from_evzinb(object, newdata = newdata)
  ## Estimate mu_nb for all individuals
  cnts <- counts_from_evzinb(object, newdata = newdata)
  ## Estimate pareto alpha for all individuals
  alphs <- fitted_alpha_from_evzinb(object, newdata = newdata)
  if (min(alphs$pareto_alpha) < 1e-02) {
    warning(
      'Fitted pareto alpha-values below 1e-02 detected. Setting those alphas to 1e-02 for quantile prediction'
    )
    alphs <- alphs %>%
      dplyr::mutate(
        pareto_alpha = dplyr::case_when(
          pareto_alpha < 1e-02 ~ 1e-02,
          T ~ pareto_alpha
        )
      )
  }

  all_pars <- dplyr::bind_cols(prbs, cnts, alphs) %>%
    dplyr::mutate(C = object$coef$C, alpha_nb = object$coef$Alpha.NB)

  individual_dists <- all_pars %>%
    purrr::transpose() %>%
    purrr::map(
      ~ mistr::mixdist(
        mistr::binomdist(1, 0),
        nbinomdist2(mu = .x$count, size = .x$alpha_nb),
        mistr::paretodist(scale = .x$C, shape = .x$pareto_alpha),
        weights = c(.x$pr_zc, .x$pr_count, .x$pr_pareto)
      )
    )
  if (!multicore) {
    q <- individual_dists %>%
      purrr::map(~ mistr::mistr_q(.x, quantile)) %>%
      purrr::reduce(c) %>%
      round()
  } else {
    if (is.null(ncores)) {
      doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
    } else {
      doParallel::registerDoParallel(cores = ncores)
    }
    q <- foreach::foreach(
      i = 1:length(individual_dists),
      .final = unlist
    ) %dopar%
      mistr::mistr_q(individual_dists[[i]], quantile)
  }

  if (return_data) {
    return(newdata %>% dplyr::mutate(q = q))
  } else {
    return(q)
  }
}

#' Extracting full mixture quantiles from an evinb object
#'
#' @param object  An evzinb object for which to produce quantiles
#' @param quantile The quantile for which to produce predictions
#' @param newdata  Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#' @param multicore Logical: Should parallel processing be used to obtain results
#' @param ncores Number of cores if multicore is used
#'
#' @return A vector of predicted quantiles, or if return_data=T, a tibble with the predicted quantile attached last
#'
#' @noRd
quantiles_from_evinb <- function(
  object,
  quantile,
  newdata = NULL,
  return_data = FALSE,
  multicore = TRUE,
  ncores = NULL
) {
  i <- 'temp_iter'

  ## Estimate component probabilities for all individuals
  prbs <- prob_from_evinb(object, newdata = newdata)
  ## Estimate mu_nb for all individuals
  cnts <- counts_from_evzinb(object, newdata = newdata)
  ## Estimate pareto alpha for all individuals
  alphs <- fitted_alpha_from_evzinb(object, newdata = newdata)
  # if(min(alphs)<1e-02){
  #   warning('Fitted pareto alpha-values below 1e-02 detected. Setting those alphas to 1e-02.')
  #   alphs <- alphs %>% dplyr::mutate(alpha = dplyr::case_when(alpha<1e-02 ~ 1e-02,
  #                                                             T~alpha))
  # }

  all_pars <- dplyr::bind_cols(prbs, cnts, alphs) %>%
    dplyr::mutate(C = object$coef$C, alpha_nb = object$coef$Alpha.NB)

  individual_dists <- all_pars %>%
    purrr::transpose() %>%
    purrr::map(
      ~ mistr::mixdist(
        nbinomdist2(mu = .x$count, size = .x$alpha_nb),
        mistr::paretodist(scale = .x$C, shape = .x$pareto_alpha),
        weights = c(.x$pr_count, .x$pr_pareto)
      )
    )
  if (!multicore) {
    q <- individual_dists %>%
      purrr::map(~ mistr::mistr_q(.x, quantile)) %>%
      purrr::reduce(c) %>%
      round()
  } else {
    if (is.null(ncores)) {
      doParallel::registerDoParallel(cores = parallel::detectCores() - 1)
    } else {
      doParallel::registerDoParallel(cores = ncores)
    }
    q <- foreach::foreach(
      i = 1:length(individual_dists),
      .final = unlist
    ) %dopar%
      mistr::mistr_q(individual_dists[[i]], quantile)
  }

  if (return_data) {
    return(newdata %>% dplyr::mutate(q = q))
  } else {
    return(q)
  }
}

#' Extracting state probabilities from an evzinb object
#'
#' @param object An evzinb object for which to produce probabilities
#' @param newdata Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#'
#' @return A tibble with the predicted state probabilities. If return_data=T this is appended to the data or newdata
#'
#' @noRd
prob_from_evzinb <- function(object, newdata = NULL, return_data = FALSE) {
  if (is.null(newdata)) {
    x.multinom.zc <- object$data$x.multinom.zc
    x.multinom.pl <- object$data$x.multinom.pl
  } else {
    x.multinom.zc <- as.matrix(model.frame(
      object$formulas$formula_zi[-2],
      newdata
    ))
    x.multinom.pl <- as.matrix(model.frame(
      object$formulas$formula_evi[-2],
      newdata
    ))
  }

  pr_zc <- exp(cbind(1, x.multinom.zc) %*% object$coef$Beta.multinom.ZC) /
    (1 +
      exp(cbind(1, x.multinom.zc) %*% object$coef$Beta.multinom.ZC) +
      exp(cbind(1, x.multinom.pl) %*% object$coef$Beta.multinom.PL))

  pr_pareto <- exp(cbind(1, x.multinom.pl) %*% object$coef$Beta.multinom.PL) /
    (1 +
      exp(cbind(1, x.multinom.zc) %*% object$coef$Beta.multinom.ZC) +
      exp(cbind(1, x.multinom.pl) %*% object$coef$Beta.multinom.PL))

  pr_count <- 1 - pr_zc - pr_pareto

  if (min(pr_count) < 0) {
    if (min(pr_count < -1e20)) {
      stop('Error in prediction, negative probabilities produced')
    } else {
      pr_count[which(pr_count < 0)] <- 0
    }
  }
  out <- tibble::tibble(
    pr_zc = as.numeric(pr_zc),
    pr_count = as.numeric(pr_count),
    pr_pareto = as.numeric(pr_pareto)
  )

  if (return_data) {
    if (is.null(newdata)) {
      out <- dplyr::bind_cols(object$data$data, out)
    } else {
      out <- dplyr::bind_cols(newdata, out)
    }
  }

  return(out)
}

#' Extracting state probabilities from an evzinb object
#'
#' @param object An evzinb object for which to produce probabilities
#' @param newdata Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#'
#' @return A tibble with the predicted state probabilities. If return_data=T this is appended to the data or newdata
#'
#' @noRd
prob_from_evinb <- function(object, newdata = NULL, return_data = FALSE) {
  if (is.null(newdata)) {
    x.multinom.pl <- object$data$x.multinom.pl
  } else {
    x.multinom.pl <- as.matrix(model.frame(
      object$formulas$formula_evi[-2],
      newdata
    ))
  }

  # pr_zc <- exp(cbind(1,x.multinom.zc)%*%object$coef$Beta.multinom.ZC)/
  #   (1+exp(cbind(1,x.multinom.zc)%*%object$coef$Beta.multinom.ZC)+
  #      exp(cbind(1,x.multinom.pl)%*%object$coef$Beta.multinom.PL))
  #

  pr_pareto <- exp(cbind(1, x.multinom.pl) %*% object$coef$Beta.multinom.PL) /
    (1 + exp(cbind(1, x.multinom.pl) %*% object$coef$Beta.multinom.PL))

  pr_count <- 1 - pr_pareto

  out <- tibble::tibble(
    pr_count = as.numeric(pr_count),
    pr_pareto = as.numeric(pr_pareto)
  )

  if (return_data) {
    if (is.null(newdata)) {
      out <- dplyr::bind_cols(object$data$data, out)
    } else {
      out <- dplyr::bind_cols(newdata, out)
    }
  }

  return(out)
}

#' Extracting fitted count values of the NB component of an evzinb object
#'
#' @param object An evzinb object for which to produce counts
#' @param newdata Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#'
#' @return A tibble with the predicted nb counts. If return_data=T this is appended to the data or newdata
#'
#' @noRd
counts_from_evzinb <- function(object, newdata = NULL, return_data = FALSE) {
  if (is.null(newdata)) {
    x.nb <- object$data$x.nb
  } else {
    x.nb <- as.matrix(model.frame(object$formulas$formula_nb[-2], newdata))
  }

  count <- exp(cbind(1, x.nb) %*% object$coef$Beta.NB)

  out <- tibble::tibble(count = as.numeric(count))

  if (return_data) {
    if (is.null(newdata)) {
      out <- dplyr::bind_cols(object$data$data, out)
    } else {
      out <- dplyr::bind_cols(newdata, out)
    }
  }
  return(out)
}

#' Extracting fitted alpha values of the pareto component of an evzinb object
#'
#' @param object An evzinb or evinb object for which to produce pareto_alphas
#' @param newdata Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#'
#' @return A tibble with the predicted pareto alpha values. If return_data=T this is appended to the data or newdata
#'
#' @noRd
fitted_alpha_from_evzinb <- function(
  object,
  newdata = NULL,
  return_data = FALSE
) {
  if (is.null(newdata)) {
    x.pl <- object$data$x.pl
  } else {
    x.pl <- as.matrix(model.frame(object$formulas$formula_pareto[-2], newdata))
  }

  alpha <- exp(cbind(1, x.pl) %*% object$coef$Beta.PL)

  out <- tibble::tibble(pareto_alpha = as.numeric(alpha))

  if (return_data) {
    if (is.null(newdata)) {
      out <- dplyr::bind_cols(object$data$data, out)
    } else {
      out <- dplyr::bind_cols(newdata, out)
    }
  }
  return(out)
}
