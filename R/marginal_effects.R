# audit 4.8 - average marginal effects with bootstrap confidence intervals.

# Prediction from a single (full or bootstrap) fit. `mod` may be a bootstrap
# object, which carries $coef / $terms / $xlevels / $formulas.
predict_from_boot <- function(mod, newdata, type, quantile = NULL, evzinb = TRUE) {
  pf <- if (evzinb) predict.evzinb else predict.evinb
  if (type == "states") {
    st <- pf(mod, newdata = newdata, type = "states")
    cols <- intersect(c("pr_zero", "pr_count", "pr_evi"), names(st))
    return(as.matrix(st[, cols, drop = FALSE]))
  }
  if (type == "quantile") {
    # continuous (interpolated) mixture quantile: the integer step function has
    # a zero derivative almost everywhere (audit N1).
    qf <- if (evzinb) quantiles_from_evzinb else quantiles_from_evinb
    return(as.numeric(qf(mod, quantile, newdata = newdata,
                         return_data = FALSE, multicore = FALSE, round = FALSE)))
  }
  as.numeric(pf(mod, newdata = newdata, type = type))
}

.me_colmeans <- function(x) if (is.matrix(x)) colMeans(x) else mean(x)

# sample.int() with an optional seed, leaving the caller's RNG state untouched.
evinf_seeded_sample <- function(n, size, seed = NULL) {
  if (is.null(seed)) {
    return(sample.int(n, size))
  }
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) saved <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", saved, envir = .GlobalEnv)
    } else {
      suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
    }
  }, add = TRUE)
  set.seed(seed)
  sample.int(n, size)
}

# AME of one variable for one fit; returns a named vector (length 1 for scalar
# types, length 3 for "states").
evinf_ame_one <- function(mod, newdata, variable, type, quantile, method, eps,
                          delta, evzinb) {
  col <- newdata[[variable]]
  if (is.numeric(col)) {
    if (method == "difference") {
      # average change from a `delta`-unit increase in the covariate
      nd1 <- newdata; nd1[[variable]] <- col + delta
      d <- predict_from_boot(mod, nd1, type, quantile, evzinb) -
        predict_from_boot(mod, newdata, type, quantile, evzinb)
    } else {
      # central finite difference (derivative)
      ndp <- newdata; ndp[[variable]] <- col + eps
      ndm <- newdata; ndm[[variable]] <- col - eps
      d <- (predict_from_boot(mod, ndp, type, quantile, evzinb) -
              predict_from_boot(mod, ndm, type, quantile, evzinb)) / (2 * eps)
    }
    out <- .me_colmeans(d)
    return(stats::setNames(out, if (length(out) == 1L) "dydx" else names(out)))
  }
  # factor / character: difference of each level vs the reference (first) level
  lv <- if (is.factor(col)) levels(col) else sort(unique(as.character(col)))
  ref <- lv[1]
  res <- lapply(lv[-1], function(l) {
    ndr <- newdata; ndl <- newdata
    ndr[[variable]] <- if (is.factor(col)) factor(ref, lv) else ref
    ndl[[variable]] <- if (is.factor(col)) factor(l, lv) else l
    d <- predict_from_boot(mod, ndl, type, quantile, evzinb) -
      predict_from_boot(mod, ndr, type, quantile, evzinb)
    .me_colmeans(d)
  })
  names(res) <- paste0(lv[-1], " - ", ref)
  res
}

#' Average marginal effects for an evzinb / evinb model
#'
#' Average, over the data, of the numerical effect of each covariate on the
#' chosen predicted quantity, with bootstrap confidence intervals obtained by
#' recomputing the effect on every bootstrap fit.
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model with bootstraps.
#' @param variables Covariates to compute effects for (default: all raw
#'   covariates).
#' @param type \code{"harmonic"}, \code{"states"} or \code{"quantile"}.
#' @param quantile Quantile for \code{type = "quantile"}.
#' @param at Named list of values at which to hold covariates before
#'   differencing.
#' @param method For numeric covariates, \code{"derivative"} (central finite
#'   difference of step \code{eps}) or \code{"difference"} (average change from
#'   a \code{delta}-unit increase in the covariate). Defaults to
#'   \code{"derivative"}, except for \code{type = "quantile"} where it defaults
#'   to \code{"difference"} (a one-unit change is what a quantile effect means
#'   in practice). Factor covariates always use level-vs-reference differences.
#' @param eps Step size for the central finite difference
#'   (\code{method = "derivative"}).
#' @param delta Covariate increase for \code{method = "difference"} (default
#'   \code{1} unit; pass e.g. \code{sd(x)} for a one-SD change).
#' @param conf_level Confidence level for the bootstrap intervals.
#' @param newdata Data to average over (default: the estimation data).
#' @inheritParams evzinb
#' @param n_max For \code{type = "quantile"} only: cap the number of rows the
#'   effect is averaged over (the per-observation quantile machinery is slow).
#'   The cap is applied automatically when \code{nrow(newdata) * n_bootstraps >
#'   2e5}, or whenever \code{n_max} is passed explicitly; \code{n_max = Inf}
#'   disables it. The subsample is reproducible from the model's bootstrap seed.
#' @param exclude_degenerate Drop bootstrap replicates flagged degenerate
#'   (default \code{TRUE}); see the \code{alpha_floor} argument of
#'   \code{\link{evinf_control}}.
#'
#' @inheritSection evzinb Parallel processing
#'
#' @details
#' The confidence interval is a percentile interval: \code{conf.low} /
#' \code{conf.high} are the empirical quantiles of the effect recomputed on each
#' bootstrap fit. \code{std.error} is the standard deviation of those recomputed
#' effects, except for \code{type = "quantile"} where the bootstrap distribution
#' is heavy-tailed and \code{std.error} is instead a robust scale read off the
#' interval, \code{(conf.high - conf.low) / (2 * qnorm(1 - (1 - conf_level)/2))}.
#' It is reported for reference and is not used to build the interval.
#'
#' \strong{Harmonic-mean effects and bootstrap intervals.} The harmonic-mean
#' prediction is \eqn{C(1 + \alpha)/\alpha}, which blows up on bootstrap fits
#' whose Pareto shape \eqn{\alpha} is very small. The point estimate for
#' \code{type = "harmonic"} is well behaved, but the percentile interval can be
#' extremely wide because a handful of bootstrap draws are enormous. For
#' inference about how a covariate shifts the outcome, prefer
#' \code{type = "quantile"} (effect on a predicted quantile) or
#' \code{type = "states"} (effect on the state probabilities), whose bootstrap
#' distributions are bounded.
#'
#' \strong{Quantile effects.} The mixture quantile is a step function of a
#' count, so its derivative is either 0 or huge; \code{type = "quantile"}
#' therefore defaults to \code{method = "difference"} (the average change in the
#' predicted quantile from a \code{delta}-unit increase in the covariate). The
#' derivative is still available with \code{method = "derivative"} but its
#' bootstrap distribution can be very heavy-tailed. The quantile effect is
#' evaluated per observation by bisection on the mixture CDF, which is still
#' the slowest prediction type, so it is averaged over at most \code{n_max}
#' rows (see that argument).
#'
#' @return A tibble with columns \code{variable}, \code{contrast}, \code{type},
#'   \code{estimate}, \code{std.error}, \code{conf.low}, \code{conf.high}.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 10)
#' marginal_effects(model, variables = "x1")
#' }
marginal_effects <- function(object, variables = NULL,
                             type = c("harmonic", "states", "quantile"),
                             quantile = NULL, at = list(),
                             method = c("derivative", "difference"),
                             eps = 1e-4, delta = 1, conf_level = 0.9,
                             newdata = NULL, n_max = 500,
                             exclude_degenerate = TRUE,
                             multicore = NULL, ncores = NULL) {
  type <- match.arg(type)
  # type = "quantile" defaults to a one-unit difference (the derivative of a
  # quantile is numerically unstable on bootstrap fits with tiny Pareto alpha).
  method <- if (missing(method) && type == "quantile") "difference"
            else match.arg(method)
  if (!inherits(object, c("evzinb", "evinb"))) {
    stop("`object` must be a fitted evzinb / evinb model.", call. = FALSE)
  }
  if (is.null(object$bootstraps)) {
    stop("marginal_effects() needs a model fitted with bootstrap = TRUE.",
         call. = FALSE)
  }
  if (type == "quantile" && is.null(quantile)) {
    stop("`quantile` is required for type = \"quantile\".", call. = FALSE)
  }
  evzinb <- inherits(object, "evzinb")

  nd <- if (is.null(newdata)) object$data$data else newdata
  for (v in names(at)) nd[[v]] <- at[[v]]

  raw_vars <- setdiff(
    unique(unlist(lapply(object$formulas, all.vars))),
    all.vars(object$formulas$formula_nb)[1]
  )
  if (is.null(variables)) variables <- raw_vars
  variables <- intersect(variables, raw_vars)

  boots <- evinf_usable_bootstraps(object, exclude_degenerate)
  qs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)

  # audit N8: the per-observation mixture-quantile machinery is slow; subsample
  # the rows the quantile AME averages over -- automatically when the work is
  # large, or whenever `n_max` was set explicitly. `n_max = Inf` disables it.
  # The subsample is reproducible from the model's bootstrap seed.
  n_max_set <- !missing(n_max)
  if (type == "quantile" && is.finite(n_max) && nrow(nd) > n_max &&
      (n_max_set || nrow(nd) * length(boots) > 2e5)) {
    seed <- if (length(object$boot_seeds)) object$boot_seeds[[1]] else NULL
    keep <- sort(evinf_seeded_sample(nrow(nd), n_max, seed))
    nd <- nd[keep, , drop = FALSE]
    message("marginal_effects(): averaging the quantile effect over a random ",
            n_max, "-row subsample of `newdata` (", length(boots),
            " bootstraps). Pass n_max = Inf to use all rows.")
  }

  # full-model effect per variable (sequential; one call each)
  est_by_var <- stats::setNames(
    lapply(variables, function(v)
      evinf_ame_one(object, nd, v, type, quantile, method, eps, delta, evzinb)),
    variables
  )

  # per-bootstrap effects, over the variable x bootstrap grid, in parallel
  grid <- expand.grid(vi = seq_along(variables), bi = seq_along(boots))
  boot_seed <- if (length(object$boot_seeds)) object$boot_seeds[[1]] %||% 1L else 1L
  flat <- if (nrow(grid)) {
    evinf_with_plan(multicore, ncores, {
      evinf_pmap(
        seq_len(nrow(grid)),
        function(k, grid, variables, boots, nd, type, quantile, method, eps,
                 delta, evzinb) {
          v <- variables[grid$vi[k]]
          tryCatch(
            evinf_ame_one(boots[[grid$bi[k]]], nd, v, type, quantile, method,
                          eps, delta, evzinb),
            error = function(e) NULL)
        },
        grid = grid, variables = variables, boots = boots, nd = nd,
        type = type, quantile = quantile, method = method, eps = eps,
        delta = delta, evzinb = evzinb,
        seed = boot_seed, label = "marginal_effects bootstrap"
      )
    })
  } else {
    list()
  }

  rows <- list()
  for (vi in seq_along(variables)) {
    v <- variables[vi]
    est <- est_by_var[[v]]
    boot_vals <- flat[grid$vi == vi]
    boot_vals <- boot_vals[!vapply(boot_vals, is.null, logical(1))]

    contrasts <- if (is.list(est)) names(est) else "dydx"
    est_list <- if (is.list(est)) est else stats::setNames(list(est), "dydx")

    for (ct in contrasts) {
      e <- est_list[[ct]]
      types <- if (length(e) == 1L) type else names(e)
      bmat <- if (length(boot_vals)) {
        do.call(rbind, lapply(boot_vals, function(bv) {
          x <- if (is.list(bv)) bv[[ct]] else bv
          if (is.null(x)) rep(NA_real_, length(e)) else x
        }))
      } else {
        matrix(NA_real_, nrow = 0, ncol = length(e))
      }
      lo <- if (nrow(bmat) > 1) apply(bmat, 2L, stats::quantile, qs[1],
                                     na.rm = TRUE, names = FALSE)
            else rep(NA_real_, length(e))
      hi <- if (nrow(bmat) > 1) apply(bmat, 2L, stats::quantile, qs[2],
                                     na.rm = TRUE, names = FALSE)
            else rep(NA_real_, length(e))
      se <- if (nrow(bmat) <= 1) {
        rep(NA_real_, length(e))
      } else if (type == "quantile") {
        # the bootstrap distribution of a quantile effect is heavy-tailed
        # (a few near-degenerate fits dominate the variance); report a robust
        # scale from the percentile interval instead of the raw sd.
        (hi - lo) / (2 * stats::qnorm(qs[2]))
      } else {
        apply(bmat, 2L, stats::sd, na.rm = TRUE)
      }

      rows[[length(rows) + 1L]] <- tibble::tibble(
        variable = v, contrast = ct, type = types,
        estimate = as.numeric(e), std.error = as.numeric(se),
        conf.low = as.numeric(lo), conf.high = as.numeric(hi)
      )
    }
  }
  dplyr::bind_rows(rows)
}
