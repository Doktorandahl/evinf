# audit 4.6 - prediction grids for effect displays.

# Held-fixed value for a covariate: mean/median for numerics, modal level for
# factors/characters/logicals.
evinf_hold_fixed <- function(x, fixed) {
  if (is.numeric(x)) {
    return(if (fixed == "median") stats::median(x, na.rm = TRUE)
           else mean(x, na.rm = TRUE))
  }
  tab <- sort(table(x), decreasing = TRUE)
  lev <- names(tab)[1]
  if (is.factor(x)) factor(lev, levels = levels(x)) else methods::as(lev, class(x)[1])
}

#' Build a prediction grid and evaluate a model over it
#'
#' Varies one covariate over a range (or set of values) while holding the others
#' at representative values, and returns the model predictions in long form -
#' the raw material for an effect plot.
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model.
#' @param variable Name of the covariate to vary.
#' @param values Optional explicit values for \code{variable}; otherwise an
#'   evenly spaced sequence (numeric) or all levels (factor).
#' @param n Number of grid points for a numeric \code{variable}.
#' @param at Named list of values at which to hold specific other covariates.
#' @param fixed How to hold the remaining numeric covariates: \code{"mean"} or
#'   \code{"median"}.
#' @param type Prediction type: \code{"states"}, \code{"harmonic"},
#'   \code{"quantile"}, \code{"counts"} or \code{"pareto_alpha"}.
#' @param quantile Quantile for \code{type = "quantile"}.
#' @param confint Add bootstrap confidence intervals (not for \code{"states"}).
#' @param conf_level Confidence level.
#'
#' @return A tibble with columns \code{variable}, \code{value}, \code{type},
#'   \code{estimate}, \code{conf.low}, \code{conf.high} (the last two \code{NA}
#'   when \code{confint = FALSE}).
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' predict_grid(model, "x1", type = "harmonic")
#' }
predict_grid <- function(object, variable, values = NULL, n = 50, at = list(),
                         fixed = c("mean", "median"),
                         type = c("states", "harmonic", "quantile", "counts",
                                  "pareto_alpha"),
                         quantile = NULL, confint = FALSE, conf_level = 0.9) {
  fixed <- match.arg(fixed)
  type <- match.arg(type)
  if (!inherits(object, c("evzinb", "evinb"))) {
    stop("`object` must be a fitted evzinb / evinb model.", call. = FALSE)
  }
  dat <- object$data$data
  raw_vars <- setdiff(
    unique(unlist(lapply(object$formulas, all.vars))),
    all.vars(object$formulas$formula_nb)[1]
  )
  if (!variable %in% raw_vars) {
    stop("`variable` must be one of: ", paste(raw_vars, collapse = ", "), ".",
         call. = FALSE)
  }

  if (is.null(values)) {
    col <- dat[[variable]]
    values <- if (is.numeric(col)) {
      seq(min(col, na.rm = TRUE), max(col, na.rm = TRUE), length.out = n)
    } else if (is.factor(col)) {
      factor(levels(col), levels = levels(col))
    } else {
      sort(unique(col))
    }
  }

  grid <- tibble::tibble(!!variable := values)
  for (v in setdiff(raw_vars, variable)) {
    grid[[v]] <- if (v %in% names(at)) at[[v]] else
      evinf_hold_fixed(dat[[v]], fixed)
  }

  if (type == "states") {
    if (confint) {
      message("predict_grid(): confidence intervals are not available for ",
              "type = \"states\".")
    }
    st <- stats::predict(object, newdata = grid, type = "states")
    keep <- c("pr_zero", "pr_count", "pr_evi")
    keep <- intersect(keep, names(st))
    long <- tidyr::pivot_longer(
      dplyr::bind_cols(tibble::tibble(value = values), st[, keep, drop = FALSE]),
      cols = dplyr::all_of(keep), names_to = "type", values_to = "estimate"
    )
    return(tibble::tibble(
      variable = variable, value = long$value, type = long$type,
      estimate = long$estimate, conf.low = NA_real_, conf.high = NA_real_
    ))
  }

  args <- list(object, newdata = grid, type = type)
  if (type == "quantile") {
    if (is.null(quantile)) {
      stop("`quantile` is required for type = \"quantile\".", call. = FALSE)
    }
    args$quantile <- quantile
  }
  if (confint) {
    args$confint <- TRUE
    args$conf_level <- conf_level
  }
  pred <- do.call(stats::predict, args)

  if (confint) {
    est_col <- setdiff(names(pred), c("ci_lb", "ci_ub"))[1]
    out <- tibble::tibble(
      variable = variable, value = values, type = type,
      estimate = pred[[est_col]], conf.low = pred$ci_lb, conf.high = pred$ci_ub
    )
  } else {
    out <- tibble::tibble(
      variable = variable, value = values, type = type,
      estimate = as.numeric(pred), conf.low = NA_real_, conf.high = NA_real_
    )
  }
  out
}
