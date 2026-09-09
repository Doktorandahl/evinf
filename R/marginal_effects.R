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
    return(as.numeric(pf(mod, newdata = newdata, type = "quantile",
                         quantile = quantile)))
  }
  as.numeric(pf(mod, newdata = newdata, type = type))
}

.me_colmeans <- function(x) if (is.matrix(x)) colMeans(x) else mean(x)

# AME of one variable for one fit; returns a named vector (length 1 for scalar
# types, length 3 for "states").
evinf_ame_one <- function(mod, newdata, variable, type, quantile, method, eps,
                          evzinb) {
  col <- newdata[[variable]]
  if (is.numeric(col)) {
    ndp <- newdata; ndp[[variable]] <- col + eps
    ndm <- newdata; ndm[[variable]] <- col - eps
    d <- (predict_from_boot(mod, ndp, type, quantile, evzinb) -
            predict_from_boot(mod, ndm, type, quantile, evzinb)) / (2 * eps)
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
#' @param method \code{"derivative"} (central difference, numeric variables) or
#'   \code{"difference"}; factors always use level-vs-reference differences.
#' @param eps Step size for the central difference.
#' @param conf_level Confidence level for the bootstrap intervals.
#' @param newdata Data to average over (default: the estimation data).
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
                             eps = 1e-4, conf_level = 0.9, newdata = NULL) {
  type <- match.arg(type)
  method <- match.arg(method)
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

  boots <- object$bootstraps[!vapply(object$bootstraps, inherits, logical(1),
                                     "try-error")]
  qs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)

  rows <- list()
  for (v in variables) {
    est <- evinf_ame_one(object, nd, v, type, quantile, method, eps, evzinb)
    boot_vals <- lapply(boots, function(b)
      tryCatch(evinf_ame_one(b, nd, v, type, quantile, method, eps, evzinb),
               error = function(e) NULL))
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
      se <- if (nrow(bmat) > 1) apply(bmat, 2L, stats::sd, na.rm = TRUE)
            else rep(NA_real_, length(e))
      lo <- if (nrow(bmat) > 1) apply(bmat, 2L, stats::quantile, qs[1],
                                     na.rm = TRUE, names = FALSE)
            else rep(NA_real_, length(e))
      hi <- if (nrow(bmat) > 1) apply(bmat, 2L, stats::quantile, qs[2],
                                     na.rm = TRUE, names = FALSE)
            else rep(NA_real_, length(e))

      rows[[length(rows) + 1L]] <- tibble::tibble(
        variable = v, contrast = ct, type = types,
        estimate = as.numeric(e), std.error = as.numeric(se),
        conf.low = as.numeric(lo), conf.high = as.numeric(hi)
      )
    }
  }
  dplyr::bind_rows(rows)
}
