# round10 H.2-H.6 (audit sec5.8, dev/claude_code_round10_prompt.md): forecasting
# helpers built on evinf_pmf()/evinf_cdf() (R/evinf_pmf.R). New dedicated
# functions rather than extensions of evinf_predict_engine() (R/predict_evzinb.R)
# -- that engine is shared by every existing type and is the highest-risk file
# to touch; predict.evzinb()/predict.evinb() dispatch to these first and fall
# through to the untouched engine otherwise. This also satisfies H.3's own
# backward-compat requirement directly: a length-1 `quantile` never reaches
# this file at all, so it is unconditionally identical to today.
#
# Every function below takes an explicit `evzinb` boolean from the caller's
# own S3 dispatch (predict.evzinb()/predict.evinb() hardcode TRUE/FALSE, the
# same pattern evinf_predict_per_boot()/predict_from_boot() already use) and
# threads it through as `is_zinb` to evinf_dist_params()/evinf_pmf()/
# evinf_cdf() -- never inherits() -- since a bootstrap replicate does not
# reliably carry the "evzinb"/"evinb" class (see R/evinf_pmf.R).

.evinf_predict_types_evzinb <- c(
  'harmonic', 'explog', 'counts', 'pareto_alpha', 'zi', 'evinf',
  'count_state', 'states', 'all', 'quantile',
  'distribution', 'exceedance', 'draws'
)
.evinf_predict_types_evinb <- c(
  'harmonic', 'explog', 'counts', 'pareto_alpha', 'evinf',
  'count_state', 'states', 'all', 'quantile',
  'distribution', 'exceedance', 'draws'
)

# Shared entry point for predict.evzinb()/predict.evinb(): the three new
# H.2/H.4/H.5 types, and type = "quantile" with length(quantile) > 1 (H.3),
# are handled by the dedicated functions in this file; every other type
# (including a length-1 or NULL `quantile`) falls through to the untouched
# evinf_predict_engine() (R/predict_evzinb.R), so existing behavior for
# every pre-round10 type is exactly what it was before this file existed.
evinf_predict_dispatch <- function(object, newdata, type, pred, quantile,
                                   confint, conf_level, multicore, ncores,
                                   return_bootstraps, exclude_degenerate,
                                   support, max_support, format, threshold,
                                   n_draws, parameter_uncertainty, seed, keep,
                                   evzinb) {
  choices <- if (evzinb) .evinf_predict_types_evzinb else .evinf_predict_types_evinb
  type <- normalize_predict_type(type, choices)

  if (identical(type, 'quantile') && !is.null(quantile) && length(quantile) > 1L) {
    return(evinf_predict_quantile_vec(
      object, newdata = newdata, quantile = quantile, confint = confint,
      conf_level = conf_level, multicore = multicore, ncores = ncores,
      return_bootstraps = return_bootstraps, exclude_degenerate = exclude_degenerate,
      keep = keep, evzinb = evzinb
    ))
  }
  if (identical(type, 'distribution')) {
    if (isTRUE(confint)) {
      stop(
        "Confidence interval prediction is not available for type = ",
        "'distribution' (the full predictive distribution is already the ",
        "quantity a CI would summarize; use type = 'quantile' or ",
        "'exceedance' with confint = TRUE instead).", call. = FALSE
      )
    }
    return(evinf_predict_distribution(
      object, newdata = newdata, support = support, max_support = max_support,
      format = format, keep = keep, evzinb = evzinb
    ))
  }
  if (identical(type, 'exceedance')) {
    return(evinf_predict_exceedance(
      object, newdata = newdata, threshold = threshold, confint = confint,
      conf_level = conf_level, multicore = multicore, ncores = ncores,
      return_bootstraps = return_bootstraps, exclude_degenerate = exclude_degenerate,
      keep = keep, evzinb = evzinb
    ))
  }
  if (identical(type, 'draws')) {
    return(evinf_predict_draws(
      object, newdata = newdata, n_draws = n_draws,
      parameter_uncertainty = parameter_uncertainty, seed = seed, keep = keep,
      evzinb = evzinb
    ))
  }

  evinf_predict_engine(
    object, newdata = newdata, type = type, pred = pred, quantile = quantile,
    confint = confint, conf_level = conf_level, multicore = multicore,
    ncores = ncores, return_bootstraps = return_bootstraps,
    exclude_degenerate = exclude_degenerate, evzinb = evzinb
  )
}

# round10 H.6: carry identifier columns from newdata into a result, for a
# join key on panel data rather than row-order matching. `idx` maps each
# output row to its newdata row (defaults to 1:1 for one-row-per-observation
# results; a long result with several rows per observation, e.g. H.2's
# distribution or H.5's draws, passes the repeated index it already built).
evinf_bind_keep <- function(out, newdata, keep, idx = NULL) {
  if (is.null(keep)) {
    return(out)
  }
  missing_cols <- setdiff(keep, names(newdata))
  if (length(missing_cols) > 0L) {
    stop(
      "`keep` column", if (length(missing_cols) > 1L) "s " else " ",
      "not found in newdata: ", paste(sQuote(missing_cols), collapse = ", "),
      call. = FALSE
    )
  }
  if (is.null(idx)) {
    idx <- seq_len(nrow(out))
  }
  dplyr::bind_cols(tibble::as_tibble(newdata)[idx, keep, drop = FALSE], out)
}

#' The full predictive distribution: a long tibble or an n x K matrix
#' (round10 H.2)
#'
#' @param support Optional shared support (a vector of y values) for every
#'   row. Defaults to \code{0:K}, \code{K} the ceiling of the 0.999 mixture
#'   quantile across all rows (so the heaviest-tailed row sets the support
#'   for every row), capped at \code{max_support} (a message names the cap
#'   when it bites).
#' @param max_support Upper bound on the default support's \code{K} (default
#'   1e5); ignored when \code{support} is supplied explicitly.
#' @param format \code{"long"} (default; a tibble with \code{.row}, \code{y},
#'   \code{prob}) or \code{"matrix"} (the raw \eqn{n \times K} matrix
#'   \code{evinf_pmf()} returns, columns named by \code{support}).
#' @param keep Optional character vector of \code{newdata} column names to
#'   carry into the \code{"long"} result (round10 H.6); ignored (with a
#'   message) for \code{format = "matrix"}, which has nowhere to put them.
#' @noRd
evinf_predict_distribution <- function(object, newdata, support = NULL,
                                       max_support = 1e5,
                                       format = c("long", "matrix"),
                                       keep = NULL, evzinb) {
  format <- match.arg(format)
  if (is.null(newdata)) {
    newdata <- object$data$data
  }
  n <- nrow(newdata)

  if (is.null(support)) {
    p <- evinf_dist_params(object, newdata = newdata, is_zinb = evzinb)
    q999 <- mixture_quantile(
      0.999, p$pl_alpha, p$C, p$nb_mu, p$nb_alpha, p$probabilities,
      family_count = p$family$count, family_zero = p$family$zero
    )
    K <- max(q999[is.finite(q999)])
    if (K > max_support) {
      message(
        "evinf (predict distribution): the default support's upper bound (",
        K, ", the 0.999 mixture quantile of the heaviest-tailed row) ",
        "exceeds max_support (", max_support, "); truncating there. Pass ",
        "`support = ` explicitly for a wider (or narrower) support."
      )
      K <- max_support
    }
    support <- 0:K
  }

  pmf <- evinf_pmf(object, newdata = newdata, support = support, is_zinb = evzinb)

  if (format == "matrix") {
    if (!is.null(keep)) {
      message("evinf (predict distribution): `keep` is ignored for format = \"matrix\".")
    }
    dimnames(pmf) <- list(NULL, as.character(support))
    return(pmf)
  }

  idx <- rep(seq_len(n), each = length(support))
  out <- tibble::tibble(
    .row = idx,
    y = rep(support, times = n),
    prob = as.vector(t(pmf))
  )
  evinf_bind_keep(out, newdata, keep, idx = idx)
}

# Shared bootstrap-CI machinery for the "vector of scalar predictions per
# row" shape (H.3's several quantiles, H.4's several thresholds): `fn(dp)`
# computes an n x length(cols) matrix (or a length-n vector for a single
# column) of point estimates from one evinf_dist_params() result; this runs
# it once on `object` for the point estimate and, if `confint`, once per
# usable bootstrap replicate for the interval.
evinf_boot_ci_matrix <- function(object, newdata, fn, col_names, confint,
                                 conf_level, multicore, ncores,
                                 return_bootstraps, exclude_degenerate, evzinb) {
  p <- evinf_dist_params(object, newdata = newdata, is_zinb = evzinb)
  est <- fn(p)
  if (!is.matrix(est)) {
    est <- matrix(est, ncol = 1)
  }
  colnames(est) <- col_names
  out <- tibble::as_tibble(est)

  if (!confint) {
    return(list(out = out, bootstraps = NULL))
  }

  boots <- evinf_usable_bootstraps(object, exclude_degenerate)
  if (length(boots) == 0L) {
    stop("No usable bootstrap replicates for a confidence interval.", call. = FALSE)
  }
  per_boot <- evinf_with_plan(multicore, ncores, {
    evinf_pmap(
      boots,
      function(b, newdata, fn, col_names, evzinb) {
        pb <- evinf_dist_params(b, newdata = newdata, is_zinb = evzinb)
        eb <- fn(pb)
        if (!is.matrix(eb)) {
          eb <- matrix(eb, ncol = 1)
        }
        colnames(eb) <- col_names
        eb
      },
      newdata = newdata, fn = fn, col_names = col_names, evzinb = evzinb,
      seed = object$boot_seeds[[1]] %||% 1L, label = "predict bootstrap"
    )
  })

  qs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
  for (j in seq_along(col_names)) {
    col_boot <- vapply(per_boot, function(m) m[, j], numeric(nrow(newdata)))
    out[[paste0(col_names[j], "_lo")]] <- apply(col_boot, 1, stats::quantile, probs = qs[1])
    out[[paste0(col_names[j], "_hi")]] <- apply(col_boot, 1, stats::quantile, probs = qs[2])
  }
  list(out = out, bootstraps = if (return_bootstraps) per_boot else NULL)
}

#' Several mixture quantiles at once, with optional bootstrap CI (round10 H.3)
#' @noRd
evinf_predict_quantile_vec <- function(object, newdata, quantile, confint,
                                       conf_level, multicore, ncores,
                                       return_bootstraps, exclude_degenerate,
                                       keep, evzinb) {
  if (is.null(newdata)) {
    newdata <- object$data$data
  }
  col_names <- paste0("q", 100 * quantile)
  fn <- function(dp) {
    mixture_quantile(quantile, dp$pl_alpha, dp$C, dp$nb_mu, dp$nb_alpha,
                     dp$probabilities, family_count = dp$family$count,
                     family_zero = dp$family$zero)
  }
  res <- evinf_boot_ci_matrix(
    object, newdata, fn, col_names, confint, conf_level, multicore, ncores,
    return_bootstraps, exclude_degenerate, evzinb
  )
  out <- evinf_bind_keep(res$out, newdata, keep)
  if (return_bootstraps) list(ci = out, bootstraps = res$bootstraps) else out
}

# P(Y >= k) = 1 - F(k - 1); mixture_p()/evinf_cdf() were never designed for
# x = -1, so threshold = 0 ("P(Y >= 0)", which is 1 for every count-valued
# model here) is handled directly instead of evaluated at k - 1 = -1.
evinf_exceedance_one <- function(k, dp) {
  if (k <= 0) {
    return(rep(1, length(dp$nb_mu)))
  }
  Fkm1 <- mixture_p(k - 1, dp$pl_alpha, dp$C, dp$nb_mu, dp$nb_alpha,
                    dp$probabilities, family_count = dp$family$count,
                    family_zero = dp$family$zero)
  1 - Fkm1
}

#' Exceedance probabilities P(Y >= threshold), with optional bootstrap CI
#' (round10 H.4)
#' @noRd
evinf_predict_exceedance <- function(object, newdata, threshold, confint,
                                     conf_level, multicore, ncores,
                                     return_bootstraps, exclude_degenerate,
                                     keep, evzinb) {
  if (is.null(newdata)) {
    newdata <- object$data$data
  }
  if (is.null(threshold)) {
    stop("`threshold` must be provided for exceedance prediction.", call. = FALSE)
  }
  col_names <- paste0("p_ge_", threshold)
  fn <- function(dp) {
    vapply(threshold, evinf_exceedance_one, numeric(length(dp$nb_mu)), dp = dp)
  }
  res <- evinf_boot_ci_matrix(
    object, newdata, fn, col_names, confint, conf_level, multicore, ncores,
    return_bootstraps, exclude_degenerate, evzinb
  )
  out <- evinf_bind_keep(res$out, newdata, keep)
  if (return_bootstraps) list(ci = out, bootstraps = res$bootstraps) else out
}

#' Predictive draws, long tibble (round10 H.5)
#'
#' Shares the sampler with \code{\link{simulate.evzinb}}/
#' \code{\link{simulate.evinb}} -- \code{revzinb_fit()}/\code{revinb_fit()}
#' (R/predict_evzinb.R) -- rather than a second one.
#'
#' @param parameter_uncertainty If \code{TRUE}, each draw uses a randomly
#'   chosen usable bootstrap replicate's parameters instead of the
#'   full-sample estimate.
#' @noRd
evinf_predict_draws <- function(object, newdata, n_draws, parameter_uncertainty,
                                seed, keep, evzinb) {
  if (is.null(newdata)) {
    newdata <- object$data$data
  }
  n <- nrow(newdata)
  draw_fn <- if (evzinb) revzinb_fit else revinb_fit

  rng_state <- evinf_pre_draw_rng_state(seed)

  if (!isTRUE(parameter_uncertainty)) {
    draws <- evinf_with_seed(seed, draw_fn(object, newdata = newdata, n_draws = n_draws))
    if (n_draws == 1L) {
      draws <- list(draws)
    }
  } else {
    boots <- evinf_usable_bootstraps(object, exclude_degenerate = TRUE)
    if (length(boots) == 0L) {
      stop("No usable bootstrap replicates for parameter_uncertainty = TRUE.",
          call. = FALSE)
    }
    draws <- evinf_with_seed(seed, {
      boot_idx <- sample.int(length(boots), n_draws, replace = TRUE)
      lapply(boot_idx, function(i) draw_fn(boots[[i]], newdata = newdata, n_draws = 1))
    })
  }

  idx <- rep(seq_len(n), times = n_draws)
  out <- tibble::tibble(
    .row = idx,
    .draw = rep(seq_len(n_draws), each = n),
    y = unlist(draws, use.names = FALSE)
  )
  out <- evinf_bind_keep(out, newdata, keep, idx = idx)
  attr(out, "seed") <- rng_state
  out
}
