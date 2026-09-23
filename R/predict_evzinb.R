# round9 E.2: zero-truncated draws for a hurdle count state -- rejection
# sampling redraws only the still-zero subset each round (expected number of
# rounds is 1/(1-f0), f0 = P(Y=0) under the untruncated distribution), so
# this stays cheap even when f0 is not small.
evinf_rztrunc_pois <- function(mu) {
  n <- length(mu)
  out <- rpois(n, lambda = mu)
  bad <- out == 0
  while (any(bad)) {
    out[bad] <- rpois(sum(bad), lambda = mu[bad])
    bad[bad] <- out[bad] == 0
  }
  out
}
evinf_rztrunc_nbinom <- function(mu, alpha) {
  n <- length(mu)
  out <- rnbinom(n, mu = mu, size = 1 / alpha)
  bad <- out == 0
  while (any(bad)) {
    out[bad] <- rnbinom(sum(bad), mu = mu[bad], size = 1 / alpha)
    bad[bad] <- out[bad] == 0
  }
  out
}

# round9 E.2: E[Y | count state] under a hurdle process is mu/(1-f0) (the
# zero-truncated mean), not the raw mu -- applied before harmonic_calc()/
# explog_calc() (both otherwise family-agnostic, operating on mu/count only)
# so their aggregate point predictions reflect the count state actually
# being zero-truncated. predict(type = "counts") is unaffected -- it reports
# the count state's own model parameter mu, not this conditional mean.
evinf_hurdle_mean <- function(mu, alpha_nb, family) {
  if (!identical(family$zero, "hurdle")) {
    return(mu)
  }
  f0 <- if (family$count == "poisson") {
    exp(-mu)
  } else {
    (1 + alpha_nb * mu)^(-1 / alpha_nb)
  }
  mu / (1 - f0)
}

harmonic_calc <- function(pr_count, count, pr_pareto, C, pareto_alpha,
                          floor = 0.01, warn = TRUE) {
  pareto_alpha <- evinf_clamp_alpha_pl(pareto_alpha, floor = floor,
                                       context = "harmonic prediction",
                                       warn = warn)
  # round9 0.1 follow-up: (1 + alpha) / alpha is algebraically identical to
  # 1/alpha + 1, but the former is Inf/Inf = NaN once a fitted Pareto shape
  # overflows exp() (an extreme Beta.PL coefficient on a sparse factor level
  # can push the linear predictor past ~709.78, e.g. genevzinb2_factor()'s
  # y ~ x1 + g fit). 1/alpha + 1 is exact for finite alpha and gives the
  # mathematically correct limit (1) as alpha -> Inf, so no floor/ceiling on
  # alpha is needed for this specific failure mode.
  pr_count * count + pr_pareto * C * (1 / pareto_alpha + 1)
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

# round11 A4 (dev/claude_code_round10_prompt.md §0.6, as amended): the
# explog path is exempt from the global alpha_pl_floor (unlike harmonic and
# quantile, which keep using it) -- `clamp_alpha_pl` is its own, separate
# knob. `FALSE` (default) uses the fitted alpha_pl as-is and warns when any
# value is below 0.1 (where the prediction is effectively undefined, well
# before a floor would bite); `TRUE` clamps at 0.1; a positive number clamps
# there instead. `warn = FALSE` silences the per-call warning/message (used
# for the per-bootstrap-replicate calls in evinf_predict_per_boot(), so a
# confint = TRUE call warns once instead of once per replicate).
explog_calc <- function(pr_count, count, pr_pareto, C, pareto_alpha,
                        clamp_alpha_pl = FALSE, warn = TRUE) {
  if (isFALSE(clamp_alpha_pl)) {
    bad <- pareto_alpha < 0.1
    n_bad <- sum(bad, na.rm = TRUE)
    if (warn && n_bad > 0L) {
      warning(
        "evinf (explog prediction): ", n_bad, " fitted Pareto alpha ",
        "value", if (n_bad != 1L) "s" else "", " below 0.1 (smallest: ",
        signif(min(pareto_alpha, na.rm = TRUE), 3), "); the geometric-mean ",
        "prediction C * exp(1 / alpha_pl) is effectively undefined there ",
        "(e.g. exp(10) ~= 2.2e4) and may be Inf. Consider ",
        "predict(..., type = \"explog\", clamp_alpha_pl = TRUE).",
        call. = FALSE
      )
    }
  } else {
    clamp_floor <- if (isTRUE(clamp_alpha_pl)) 0.1 else clamp_alpha_pl
    n_clamped <- sum(pareto_alpha < clamp_floor, na.rm = TRUE)
    if (warn && n_clamped > 0L) {
      message(
        "evinf (explog prediction): ", n_clamped, " fitted Pareto alpha ",
        "value", if (n_clamped != 1L) "s" else "", " clamped to ",
        clamp_floor, " (clamp_alpha_pl)."
      )
    }
    pareto_alpha <- evinf_clamp_alpha_pl(pareto_alpha, floor = clamp_floor,
                                         context = "explog prediction",
                                         warn = FALSE)
  }
  pr_count * count + C * pr_pareto * exp(1 / pareto_alpha)
}

# All per-bootstrap prediction quantities for predict.evzinb() / predict.evinb()
# in a single parallel pass (round 5). Returns a list of length length(boots),
# each element list(prbs, cnts, alphs, C, q, harmonic, explog).
evinf_predict_per_boot <- function(boots, newdata, quantile, want_q, evzinb,
                                   multicore, ncores, seed,
                                   clamp_alpha_pl = FALSE) {
  nd_id <- 1:nrow(newdata)
  evinf_with_plan(multicore, ncores, {
    evinf_pmap(
      boots,
      function(b, newdata, quantile, want_q, evzinb, nd_id, clamp_alpha_pl) {
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
          # warn = FALSE for the same reason as explog below: once per
          # bootstrap replicate would flood the user with the same warning.
          harmonic = tibble::tibble(harmonic = harmonic_calc(
            prbs$pr_count,
            evinf_hurdle_mean(cnts$count, b$coef$Alpha.NB, b$family %||% evinf_family()),
            pr_pareto = prbs$pr_pareto,
            C = C, pareto_alpha = alphs$pareto_alpha,
            floor = b$control$alpha_pl_floor %||% 0.01, warn = FALSE), id = nd_id),
          # round11 A4: warn = FALSE -- this runs once per bootstrap replicate,
          # so the user-facing warning/message comes from the single
          # point-estimate call in evinf_predict_engine() instead. Every
          # replicate is clamped identically to the same clamp_alpha_pl the
          # point estimate uses.
          explog = tibble::tibble(explog = explog_calc(
            prbs$pr_count,
            evinf_hurdle_mean(cnts$count, b$coef$Alpha.NB, b$family %||% evinf_family()),
            pr_pareto = prbs$pr_pareto,
            C = C, pareto_alpha = alphs$pareto_alpha,
            clamp_alpha_pl = clamp_alpha_pl, warn = FALSE), id = nd_id)
        )
      },
      newdata = newdata, quantile = quantile, want_q = want_q, evzinb = evzinb,
      nd_id = nd_id, clamp_alpha_pl = clamp_alpha_pl, seed = seed,
      label = "predict bootstrap"
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
                                 evzinb, clamp_alpha_pl = FALSE) {
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
  # round11 A4: the clamp actually applied, recorded as an attribute on every
  # explog-bearing result (FALSE when unclamped, else the numeric floor).
  resolved_clamp_alpha_pl <- if (isFALSE(clamp_alpha_pl)) {
    FALSE
  } else if (isTRUE(clamp_alpha_pl)) {
    0.1
  } else {
    clamp_alpha_pl
  }

  if (type %in% c('states', 'all') & confint) {
    stop('Confidence interval prediction only available for vector outputs')
  }

  if (pred %in% c('bootstrap_median', 'bootstrap_mean') | confint) {
    object$bootstraps <- evinf_usable_bootstraps(object, exclude_degenerate)
    nboots <- length(object$bootstraps)
    # round11 (found while testing A4 on hks): with zero usable replicates,
    # every per-boot list below is empty and dplyr::bind_rows() on an empty
    # list drops the `id` column entirely, so the group_by(id) further down
    # failed with a confusing "column `id` not found" instead of a clear
    # message (the same failure evinf_boot_ci_matrix() in
    # R/predict_distribution.R already guards against for the H.2-H.6 types).
    if (nboots == 0L) {
      stop(
        "No usable bootstrap replicates for pred = ", sQuote(pred),
        if (confint) " / confint = TRUE" else "", ". Every replicate was ",
        "excluded (error, or degenerate; see exclude_degenerate and ",
        "failed_bootstraps()).", call. = FALSE
      )
    }
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
      seed = object$boot_seeds[[1]] %||% 1L, clamp_alpha_pl = clamp_alpha_pl
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

    # round11 A4: only warn for the quantity the caller actually asked for --
    # harmonic_calc()'s own warning used to fire even for type = "explog" (and
    # vice versa), since both were always computed here regardless of `type`
    # (needed for type = "all"), which would have made clamp_alpha_pl = TRUE
    # still show a (harmonic-side) warning on an unrelated call.
    harmonic <- harmonic_calc(
      pr_count = prbs$pr_count,
      count = evinf_hurdle_mean(cnts$count, object$coef$Alpha.NB,
                                object$family %||% evinf_family()),
      pr_pareto = prbs$pr_pareto,
      C = C_est,
      pareto_alpha = alphs$pareto_alpha,
      floor = object$control$alpha_pl_floor %||% 0.01,
      warn = type %in% c('harmonic', 'all')
    )

    explog <- explog_calc(
      pr_count = prbs$pr_count,
      count = evinf_hurdle_mean(cnts$count, object$coef$Alpha.NB,
                                object$family %||% evinf_family()),
      pr_pareto = prbs$pr_pareto,
      C = C_est,
      pareto_alpha = alphs$pareto_alpha,
      clamp_alpha_pl = clamp_alpha_pl,
      warn = type %in% c('explog', 'all')
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
      if (type == 'explog') {
        attr(result, "clamp_alpha_pl") <- resolved_clamp_alpha_pl
      }
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
      out <- .evinf_predict_ci_spec[[type]]$estimate(environment())
      if (type == 'explog') {
        attr(out, "clamp_alpha_pl") <- resolved_clamp_alpha_pl
      }
      return(out)
    }
    if (type == 'states') {
      return(canonical_prbs(prbs))
    }
    if (type == 'quantile') {
      return(q)
    }
    if (type == 'all') {
      q_name <- paste0('q', 100 * quantile)
      out_all <- dplyr::bind_cols(
        tibble::tibble(
          harmonic = harmonic,
          explog = explog,
          !!q_name := q
        ),
        canonical_prbs(prbs),
        cnts,
        alphs
      )
      attr(out_all, "clamp_alpha_pl") <- resolved_clamp_alpha_pl
      return(out_all)
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
#'   \code{type = 'explog'} is \eqn{C \cdot \exp(1/\alpha_{PL})}, which grows
#'   explosively as the fitted Pareto shape \eqn{\alpha_{PL}} approaches 0.
#'   Unlike \code{'harmonic'} and \code{'quantile'}, this path is exempt from
#'   the global \code{evinf_control(alpha_pl_floor = )}; \code{clamp_alpha_pl}
#'   is its own, separate knob (round11 A4). With \code{clamp_alpha_pl = FALSE}
#'   (the default), the fitted \eqn{\alpha_{PL}} is used as-is and a warning
#'   fires whenever any value is below 0.1, since the prediction is
#'   effectively undefined there (\eqn{e^{10} \approx 2.2 \times 10^4}, and it
#'   may be \code{Inf}); check \code{glance()$min_alpha_pl}, and either pass
#'   \code{clamp_alpha_pl = TRUE} (clamps at 0.1) or a positive number (clamps
#'   there instead), or use \code{type = 'harmonic'}. The clamp actually
#'   applied (\code{FALSE}, or the numeric floor) is recorded as a
#'   \code{"clamp_alpha_pl"} attribute on the result.
#'
#' @param object An evzinb object for which to produce predicted values
#' @param newdata Optional new data (tibble) to produce predicted values from
#' @param type Character string, 'harmonic' for the harmonic mean and 'explog' for exponentiated expected log, 'counts' for predicted count of the negative binomial component, 'pareto_alpha' for the predicted pareto alpha value, 'states' for the predicted component states (prior), 'count_state' for predicted probability of the count state, 'evinf' for predicted probability of the pareto state,'zi' for the predicted probability of the zero state, 'all' for all predicted values, 'quantile' for quantile prediction (scalar or vector, see `quantile`), 'distribution' for the full predictive distribution, 'exceedance' for exceedance probabilities (see `threshold`), and 'draws' for predictive draws (see `n_draws`).
#' @param ... Not used; any name here is an unknown argument and errors (round11 A5), naming the valid arguments for the requested `type`. A known argument supplied for a `type` it does nothing for (e.g. `threshold` with `type = "harmonic"`) warns instead, since it has a named formal below and never reaches here.
#' @param quantile Quantile(s) for which to produce quantile prediction. A single value keeps the existing scalar behavior exactly (a vector, e.g. `c(.5, .9, .99)` (round10 H.3), returns a tibble with one `qXX` column per probability (`qXX_lo`/`qXX_hi` too, with `confint = TRUE`)).
#' @inheritParams evzinb
#' @param pred Type of prediction to be used, defaults to the original prediction from the fitted model, with alternatives being the bootstrapped median or mean. Note that bootstrap mean may yield infinite values, especially when doing quantile prediction
#' @param confint Should confidence intervals be made for the predictions? Note: only available for vector type predictions and not 'states', 'all' or 'distribution'.
#' @param conf_level What confidence level should be used for confidence intervals
#' @param return_bootstraps Should the bootstrapped predictions be returned as well? Useful for further custom analyses of the bootstrapped predictions.
#' @param exclude_degenerate Drop bootstrap replicates flagged degenerate (default TRUE); see the alpha_floor argument of evinf_control().
#' @param support (round10 H.2) `type = 'distribution'` only: shared support (vector of `y` values) for every row. `NULL` (default) uses `0:K`, `K` the ceiling of the 0.999 mixture quantile of the heaviest-tailed row, capped at `max_support`.
#' @param max_support (round10 H.2) Upper bound on the default `support`'s `K` (default 1e5); ignored when `support` is given explicitly.
#' @param format (round10 H.2) `type = 'distribution'` only: `'long'` (default; a tibble with `.row`, `y`, `prob`) or `'matrix'` (the raw n x K matrix; `keep` is ignored).
#' @param threshold (round10 H.4) `type = 'exceedance'` only: a vector of thresholds; result columns are `p_ge_<threshold>`, `P(Y >= threshold)`.
#' @param n_draws (round10 H.5) `type = 'draws'` only: number of predictive draws per row.
#' @param parameter_uncertainty (round10 H.5) `type = 'draws'` only: if `TRUE`, each draw uses a randomly chosen usable bootstrap replicate's parameters instead of the full-sample estimate.
#' @param seed (round10 H.5) `type = 'draws'` only: optional seed; the caller's `.Random.seed` is left untouched either way (same convention as `simulate()`).
#' @param keep (round10 H.6) Optional character vector of `newdata` column names to carry into the result (a join key for panel data); supported for `type` in `'quantile'` (vector), `'distribution'` (`format = 'long'`), `'exceedance'` and `'draws'`.
#' @param clamp_alpha_pl (round11 A4) `type = 'explog'` (and `'all'`) only: `FALSE` (default) uses the fitted `alpha_pl` as-is and warns below 0.1; `TRUE` clamps at 0.1; a positive number clamps there instead. Exempt from the global `alpha_pl_floor` (which stays in force for `'harmonic'`/`'quantile'`). The clamp used is recorded as a `"clamp_alpha_pl"` attribute on the result.
#'
#' @inheritSection evzinb Parallel processing
#' @inheritSection evzinb Reproducibility
#'
#' @return A vector of predicted values for type 'harmonic', 'explog', 'counts', 'pareto_alpha','zi','evinf', 'count_state', and 'quantile' (scalar), or a tibble for type 'states', 'all', 'quantile' (vector), 'distribution' (`format = 'long'`), 'exceedance', 'draws', or if confint=T; a matrix for 'distribution' with `format = 'matrix'`.
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
#' predict(model, type='quantile', quantile = c(.5, .9, .99)) # several quantiles (H.3)
#' predict(model, type='distribution') # the full predictive distribution (H.2)
#' predict(model, type='exceedance', threshold = c(10, 100)) # exceedance probs (H.4)
#' predict(model, type='draws', n_draws = 100) # predictive draws (H.5)
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
    'quantile',
    'distribution',
    'exceedance',
    'draws'
  ),
  pred = c('original', 'bootstrap_median', 'bootstrap_mean'),
  quantile = NULL,
  confint = FALSE,
  conf_level = 0.95,
  multicore = NULL,
  ncores = NULL,
  return_bootstraps = FALSE,
  exclude_degenerate = TRUE,
  support = NULL,
  max_support = 1e5,
  format = c('long', 'matrix'),
  threshold = NULL,
  n_draws = 1000,
  parameter_uncertainty = FALSE,
  seed = NULL,
  keep = NULL,
  clamp_alpha_pl = FALSE,
  ...
) {
  type_checked <- normalize_predict_type(type, .evinf_predict_types_evzinb)
  evinf_validate_predict_args(
    type_checked, list(...),
    supplied = c(
      quantile = !missing(quantile), support = !missing(support),
      max_support = !missing(max_support), format = !missing(format),
      threshold = !missing(threshold), n_draws = !missing(n_draws),
      parameter_uncertainty = !missing(parameter_uncertainty),
      seed = !missing(seed), keep = !missing(keep),
      clamp_alpha_pl = !missing(clamp_alpha_pl)
    )
  )
  evinf_predict_dispatch(
    object, newdata = newdata, type = type, pred = pred, quantile = quantile,
    confint = confint, conf_level = conf_level, multicore = multicore,
    ncores = ncores, return_bootstraps = return_bootstraps,
    exclude_degenerate = exclude_degenerate, support = support,
    max_support = max_support, format = format, threshold = threshold,
    n_draws = n_draws, parameter_uncertainty = parameter_uncertainty,
    seed = seed, keep = keep, clamp_alpha_pl = clamp_alpha_pl, evzinb = TRUE
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
#'   \code{type = 'explog'} is \eqn{C \cdot \exp(1/\alpha_{PL})}, which grows
#'   explosively as the fitted Pareto shape \eqn{\alpha_{PL}} approaches 0.
#'   Unlike \code{'harmonic'} and \code{'quantile'}, this path is exempt from
#'   the global \code{evinf_control(alpha_pl_floor = )}; \code{clamp_alpha_pl}
#'   is its own, separate knob (round11 A4). With \code{clamp_alpha_pl = FALSE}
#'   (the default), the fitted \eqn{\alpha_{PL}} is used as-is and a warning
#'   fires whenever any value is below 0.1, since the prediction is
#'   effectively undefined there (\eqn{e^{10} \approx 2.2 \times 10^4}, and it
#'   may be \code{Inf}); check \code{glance()$min_alpha_pl}, and either pass
#'   \code{clamp_alpha_pl = TRUE} (clamps at 0.1) or a positive number (clamps
#'   there instead), or use \code{type = 'harmonic'}. The clamp actually
#'   applied (\code{FALSE}, or the numeric floor) is recorded as a
#'   \code{"clamp_alpha_pl"} attribute on the result.
#'
#' @param object An evinb object for which to produce predicted values
#' @param newdata Optional new data (tibble) to produce predicted values from
#' @param type Character string, 'harmonic' for the harmonic mean and 'explog' for exponentiated expected log, 'counts' for predicted count of the negative binomial component, 'pareto_alpha' for the predicted pareto alpha value, 'states' for the predicted component states (prior), 'count_state' for predicted probability of the count state, 'evinf' for predicted probability of the pareto state, 'all' for all predicted values, 'quantile' for quantile prediction (scalar or vector, see `quantile`), 'distribution' for the full predictive distribution, 'exceedance' for exceedance probabilities (see `threshold`), and 'draws' for predictive draws (see `n_draws`).
#' @param ... Not used; any name here is an unknown argument and errors (round11 A5), naming the valid arguments for the requested `type`. A known argument supplied for a `type` it does nothing for (e.g. `threshold` with `type = "harmonic"`) warns instead, since it has a named formal below and never reaches here.
#' @param quantile Quantile(s) for which to produce quantile prediction. A single value keeps the existing scalar behavior exactly (a vector, e.g. `c(.5, .9, .99)` (round10 H.3), returns a tibble with one `qXX` column per probability (`qXX_lo`/`qXX_hi` too, with `confint = TRUE`)).
#' @inheritParams evzinb
#' @param pred Type of prediction to be used, defaults to the original prediction from the fitted model, with alternatives being the bootstrapped median or mean. Note that bootstrap mean may yield infinite values, especially when doing quantile prediction
#' @param confint Should confidence intervals be made for the predictions? Note: only available for vector type predictions and not 'states', 'all' or 'distribution'.
#' @param conf_level What confidence level should be used for confidence intervals
#' @param return_bootstraps Should the bootstrapped predictions be returned as well? Useful for further custom analyses of the bootstrapped predictions.
#' @param exclude_degenerate Drop bootstrap replicates flagged degenerate (default TRUE); see the alpha_floor argument of evinf_control().
#' @param support (round10 H.2) `type = 'distribution'` only: shared support (vector of `y` values) for every row. `NULL` (default) uses `0:K`, `K` the ceiling of the 0.999 mixture quantile of the heaviest-tailed row, capped at `max_support`.
#' @param max_support (round10 H.2) Upper bound on the default `support`'s `K` (default 1e5); ignored when `support` is given explicitly.
#' @param format (round10 H.2) `type = 'distribution'` only: `'long'` (default; a tibble with `.row`, `y`, `prob`) or `'matrix'` (the raw n x K matrix; `keep` is ignored).
#' @param threshold (round10 H.4) `type = 'exceedance'` only: a vector of thresholds; result columns are `p_ge_<threshold>`, `P(Y >= threshold)`.
#' @param n_draws (round10 H.5) `type = 'draws'` only: number of predictive draws per row.
#' @param parameter_uncertainty (round10 H.5) `type = 'draws'` only: if `TRUE`, each draw uses a randomly chosen usable bootstrap replicate's parameters instead of the full-sample estimate.
#' @param seed (round10 H.5) `type = 'draws'` only: optional seed; the caller's `.Random.seed` is left untouched either way (same convention as `simulate()`).
#' @param keep (round10 H.6) Optional character vector of `newdata` column names to carry into the result (a join key for panel data); supported for `type` in `'quantile'` (vector), `'distribution'` (`format = 'long'`), `'exceedance'` and `'draws'`.
#' @param clamp_alpha_pl (round11 A4) `type = 'explog'` (and `'all'`) only: `FALSE` (default) uses the fitted `alpha_pl` as-is and warns below 0.1; `TRUE` clamps at 0.1; a positive number clamps there instead. Exempt from the global `alpha_pl_floor` (which stays in force for `'harmonic'`/`'quantile'`). The clamp used is recorded as a `"clamp_alpha_pl"` attribute on the result.
#'
#' @inheritSection evzinb Parallel processing
#' @inheritSection evzinb Reproducibility
#'
#' @return A vector of predicted values for type 'harmonic', 'explog', 'counts', 'pareto_alpha','evinf', 'count_state', and 'quantile' (scalar), or a tibble for type 'states', 'all', 'quantile' (vector), 'distribution' (`format = 'long'`), 'exceedance', 'draws', or if confint=T; a matrix for 'distribution' with `format = 'matrix'`.
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
#' predict(model, type='quantile', quantile = c(.5, .9, .99)) # several quantiles (H.3)
#' predict(model, type='distribution') # the full predictive distribution (H.2)
#' predict(model, type='exceedance', threshold = c(10, 100)) # exceedance probs (H.4)
#' predict(model, type='draws', n_draws = 100) # predictive draws (H.5)
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
    'quantile',
    'distribution',
    'exceedance',
    'draws'
  ),
  pred = c('original', 'bootstrap_median', 'bootstrap_mean'),
  quantile = NULL,
  confint = FALSE,
  conf_level = 0.95,
  multicore = NULL,
  ncores = NULL,
  return_bootstraps = FALSE,
  exclude_degenerate = TRUE,
  support = NULL,
  max_support = 1e5,
  format = c('long', 'matrix'),
  threshold = NULL,
  n_draws = 1000,
  parameter_uncertainty = FALSE,
  seed = NULL,
  keep = NULL,
  clamp_alpha_pl = FALSE,
  ...
) {
  type_checked <- normalize_predict_type(type, .evinf_predict_types_evinb)
  evinf_validate_predict_args(
    type_checked, list(...),
    supplied = c(
      quantile = !missing(quantile), support = !missing(support),
      max_support = !missing(max_support), format = !missing(format),
      threshold = !missing(threshold), n_draws = !missing(n_draws),
      parameter_uncertainty = !missing(parameter_uncertainty),
      seed = !missing(seed), keep = !missing(keep),
      clamp_alpha_pl = !missing(clamp_alpha_pl)
    )
  )
  evinf_predict_dispatch(
    object, newdata = newdata, type = type, pred = pred, quantile = quantile,
    confint = confint, conf_level = conf_level, multicore = multicore,
    ncores = ncores, return_bootstraps = return_bootstraps,
    exclude_degenerate = exclude_degenerate, support = support,
    max_support = max_support, format = format, threshold = threshold,
    n_draws = n_draws, parameter_uncertainty = parameter_uncertainty,
    seed = seed, keep = keep, clamp_alpha_pl = clamp_alpha_pl, evzinb = FALSE
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
  fam <- object$family %||% evinf_family()
  family_count <- fam$count
  is_hurdle <- identical(fam$zero, "hurdle")

  C_est <- object$coef$C

  if (is.null(newdata)) {
    n <- length(object$data$y)
  } else {
    n <- nrow(newdata)
  }
  out <- purrr::map(seq_len(n_draws), function(draw) {
      pl_draws <- rpareto_disc(n, C_est, alphs)
      # round9 E.1: alpha_nb is NULL for a Poisson count state.
      # round9 E.2: a hurdle count state never draws 0 -- zero-truncated
      # rejection sampling instead of the raw distribution.
      count_draws <- if (family_count == "poisson") {
        if (is_hurdle) evinf_rztrunc_pois(cnts) else rpois(n, lambda = cnts)
      } else {
        if (is_hurdle) evinf_rztrunc_nbinom(cnts, alpha_nb) else rnbinom(n, mu = cnts, size = 1 / alpha_nb)
      }
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
  family_count <- (object$family %||% evinf_family())$count

  C_est <- object$coef$C
  if (is.null(newdata)) {
    n <- length(object$data$y)
  } else {
    n <- nrow(newdata)
  }

  out <- purrr::map(seq_len(n_draws), function(draw) {
      pl_draws <- rpareto_disc(n, C_est, alphs)
      # round9 E.1: alpha_nb is NULL for a Poisson count state.
      count_draws <- if (family_count == "poisson") {
        rpois(n, lambda = cnts)
      } else {
        rnbinom(n, mu = cnts, size = 1 / alpha_nb)
      }
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
