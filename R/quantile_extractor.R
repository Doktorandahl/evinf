# Clamp a fitted Pareto shape (alpha_pl) away from 0 (round9 0.1, review §2;
# originally audit0.10 §1.13/D.5 for the quantile path alone). Several
# downstream quantities involve exp(1/alpha_pl) or 1/alpha_pl --
# predict(type = "harmonic")/"explog", the continuous-Pareto mixture
# quantile's bisection, and the $fitted Pareto-tail summaries -- and silently
# return Inf or an astronomically large finite number once alpha_pl collapses
# toward 0. Shared by harmonic_calc(), explog_calc(), em_fitted_values() and
# quantiles_from_evzinb()/quantiles_from_evinb(). Warns once per call (unless
# `warn = FALSE`), naming how many observations were clamped, so the warning
# is diagnostic even when only a handful of rows are affected. em_fitted_values()
# passes warn = FALSE: it runs on every EM fit, including every bootstrap
# replicate, so warning there would fire routinely rather than only when a
# user actually queries an affected prediction; the collapse is still visible
# without digging via glance()$min_alpha_pl and the print() note.
evinf_clamp_alpha_pl <- function(alpha, floor = 0.01, context = "", warn = TRUE) {
  bad <- alpha < floor
  n_bad <- sum(bad)
  if (n_bad > 0L) {
    if (warn) {
      warning(
        "evinf", if (nzchar(context)) paste0(" (", context, ")") else "",
        ": ", n_bad, " fitted Pareto alpha value",
        if (n_bad != 1L) "s" else "",
        " below the floor (", floor, "); clamping to the floor. Set ",
        "`alpha_pl_floor` in evinf_control() to change this.",
        call. = FALSE
      )
    }
    alpha[bad] <- floor
  }
  alpha
}

#' Extracting full mixture quantiles from an evzinb object
#'
#' @param object  An evzinb object for which to produce quantiles
#' @param quantile The quantile for which to produce predictions
#' @param newdata  Optional new data (tibble) to produce predicted quantiles for
#' @param return_data Logical: Should the data be returned in the object
#' @param multicore,ncores Retained for back-compatibility; ignored. The mixture
#'   quantile is now a vectorised integer bisection on \code{mixture_p()} that
#'   runs in one pass over all rows.
#' @param round Logical: return the integer mixture quantile (the sensible
#'   default for a count-valued prediction). Pass \code{FALSE} for the
#'   continuous root of the linearly-interpolated mixture CDF, e.g. for
#'   numerical differentiation in \code{marginal_effects()}.
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
  ncores = NULL,
  round = TRUE
) {

  prbs <- prob_from_evzinb(object, newdata = newdata)
  cnts <- counts_from_evzinb(object, newdata = newdata)
  alphs <- fitted_alpha_from_evzinb(object, newdata = newdata)
  alphs$pareto_alpha <- evinf_clamp_alpha_pl(
    alphs$pareto_alpha, floor = object$control$alpha_pl_floor %||% 0.01,
    context = "quantile prediction"
  )

  q <- mixture_quantile(
    quantile,
    pl_alphas = alphs$pareto_alpha,
    C = object$coef$C,
    nb_mu = cnts$count,
    nb_alpha = object$coef$Alpha.NB,
    probabilities = cbind(prbs$pr_zc, prbs$pr_count, prbs$pr_pareto),
    continuous = !round,
    family_count = (object$family %||% evinf_family())$count
  )

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
#' @param multicore,ncores Retained for back-compatibility; ignored.
#' @param round Logical: return the integer mixture quantile. Pass \code{FALSE}
#'   for the continuous root of the linearly-interpolated mixture CDF.
#'
#' @return A vector of predicted quantiles, or if return_data=T, a tibble with the predicted quantile attached last
#'
#' @noRd
quantiles_from_evinb <- function(
  object,
  quantile,
  newdata = NULL,
  return_data = FALSE,
  multicore = FALSE,
  ncores = NULL,
  round = TRUE
) {

  prbs <- prob_from_evinb(object, newdata = newdata)
  cnts <- counts_from_evzinb(object, newdata = newdata)
  alphs <- fitted_alpha_from_evzinb(object, newdata = newdata)
  alphs$pareto_alpha <- evinf_clamp_alpha_pl(
    alphs$pareto_alpha, floor = object$control$alpha_pl_floor %||% 0.01,
    context = "quantile prediction"
  )

  q <- mixture_quantile(
    quantile,
    pl_alphas = alphs$pareto_alpha,
    C = object$coef$C,
    nb_mu = cnts$count,
    nb_alpha = object$coef$Alpha.NB,
    probabilities = cbind(0, prbs$pr_count, prbs$pr_pareto),
    continuous = !round,
    family_count = (object$family %||% evinf_family())$count
  )

  if (return_data) {
    return(newdata %>% dplyr::mutate(q = q))
  } else {
    return(q)
  }
}

# Stable 3-category softmax (implicit zero logit for the count/baseline
# category), shared by prob_from_evzinb() / prob_from_evinb() and em_fit()'s
# final state-probability recompute (audit0.10 §1.11, D.3). eta_z / eta_pl are
# numeric vectors of linear predictors (already X %*% coef); subtracting the
# row max (including the implicit 0) before exponentiating means the largest
# exp() argument is always 0, so a large coefficient (coef_limit allows up to
# 50) cannot overflow it. Passing eta_z = -Inf gives the 2-category (evinb)
# softmax as a special case (its "zero" column is then exactly 0 everywhere).
evinf_stable_props3 <- function(eta_z, eta_pl) {
  m <- pmax(0, eta_z, eta_pl)
  base <- exp(-m)
  d_z <- exp(eta_z - m)
  d_pl <- exp(eta_pl - m)
  denom <- base + d_z + d_pl
  cbind(zero = d_z / denom, count = base / denom, evi = d_pl / denom)
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
    offset_zc <- object$offset_zc %||% rep(0, nrow(x.multinom.zc))
    offset_pl_mult <- object$offset_pl_mult %||% rep(0, nrow(x.multinom.pl))
  } else {
    x.multinom.zc <- evinf_design_newdata(
      object$terms$zi, object$xlevels$zi, newdata
    )
    x.multinom.pl <- evinf_design_newdata(
      object$terms$evi, object$xlevels$evi, newdata
    )
    # round9 D.1: an offset() term requires and uses the offset variable from
    # newdata, exactly like the count component's offset already does.
    offset_zc <- evinf_offset_newdata(object$terms$zi, newdata) %||% rep(0, nrow(x.multinom.zc))
    offset_pl_mult <- evinf_offset_newdata(object$terms$evi, newdata) %||% rep(0, nrow(x.multinom.pl))
  }

  eta_zc <- as.numeric(cbind(1, x.multinom.zc) %*% object$coef$Beta.multinom.ZC) + offset_zc
  eta_pl <- as.numeric(cbind(1, x.multinom.pl) %*% object$coef$Beta.multinom.PL) + offset_pl_mult

  # pr_count is computed directly (never by subtraction), so it is a proper
  # probability by construction and cannot come out negative.
  sp <- evinf_stable_props3(eta_zc, eta_pl)

  out <- tibble::tibble(
    pr_zc = sp[, "zero"],
    pr_count = sp[, "count"],
    pr_pareto = sp[, "evi"]
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
    offset_pl_mult <- object$offset_pl_mult %||% rep(0, nrow(x.multinom.pl))
  } else {
    x.multinom.pl <- evinf_design_newdata(
      object$terms$evi, object$xlevels$evi, newdata
    )
    offset_pl_mult <- evinf_offset_newdata(object$terms$evi, newdata) %||% rep(0, nrow(x.multinom.pl))
  }

  # pr_zc <- exp(cbind(1,x.multinom.zc)%*%object$coef$Beta.multinom.ZC)/
  #   (1+exp(cbind(1,x.multinom.zc)%*%object$coef$Beta.multinom.ZC)+
  #      exp(cbind(1,x.multinom.pl)%*%object$coef$Beta.multinom.PL))
  #

  # audit0.10 §1.11: stable 2-category softmax, see evinf_stable_props3().
  eta_pl <- as.numeric(cbind(1, x.multinom.pl) %*% object$coef$Beta.multinom.PL) + offset_pl_mult
  sp <- evinf_stable_props3(rep(-Inf, length(eta_pl)), eta_pl)

  out <- tibble::tibble(
    pr_count = sp[, "count"],
    pr_pareto = sp[, "evi"]
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
    offset_nb <- object$offset_nb
  } else {
    x.nb <- evinf_design_newdata(object$terms$nb, object$xlevels$nb, newdata)
    offset_nb <- evinf_offset_newdata(object$terms$nb, newdata)  # audit 4.4
  }
  if (is.null(offset_nb)) {
    offset_nb <- rep(0, nrow(x.nb))
  }

  count <- exp(cbind(1, x.nb) %*% object$coef$Beta.NB + offset_nb)

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
    x.pl <- evinf_design_newdata(
      object$terms$pareto, object$xlevels$pareto, newdata
    )
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
