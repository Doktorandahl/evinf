# audit 4.8 - marginaleffects compatibility.

#' \code{marginaleffects} support for evzinb / evinb models
#'
#' On load the package registers the \code{"evzinb"} / \code{"evinb"} classes
#' with \pkg{marginaleffects} and provides \code{get_coef()}, \code{set_coef()},
#' \code{get_vcov()} and \code{get_predict()} methods, so
#' \code{marginaleffects::avg_slopes()}, \code{predictions()} etc. work.
#' Standard errors come from the delta method applied to the bootstrap
#' covariance matrix (\code{\link{vcov.evzinb}}).
#'
#' @details The coefficient vector includes \code{c_ev}, the extreme-value
#'   threshold, which is estimated on the grid of unique observed response
#'   values rather than by a smooth optimiser. \pkg{marginaleffects} perturbs
#'   every coefficient continuously when it builds the delta-method Jacobian, so
#'   standard errors for quantities that depend strongly on \code{c_ev} should
#'   be treated as approximate. The recommended route for uncertainty on
#'   covariate effects is the bootstrap-based \code{\link{marginal_effects}()}.
#'
#' @name marginaleffects-methods
#' @keywords internal
NULL

evinf_set_flat_coef <- function(model, coefs) {
  comps <- evinf_components(model)
  for (cn in comps) {
    slot <- .evinf_coef_slots[[cn]]
    nm <- names(model$coef[[slot]])
    model$coef[[slot]] <- stats::setNames(unname(coefs[paste0(cn, "_", nm)]), nm)
  }
  # round9 E.1: "alpha_nb" isn't in `coefs` at all for a Poisson count state
  # (evinf_flatten_coef() omits it entirely); leave model$coef$Alpha.NB
  # absent rather than setting it to NA.
  if ("alpha_nb" %in% names(coefs)) {
    model$coef$Alpha.NB <- unname(coefs[["alpha_nb"]])
  }
  model$coef$C <- unname(coefs[["c_ev"]])
  model
}

#' @exportS3Method marginaleffects::get_coef evzinb
get_coef.evzinb <- function(model, ...) coef(model, "all")
#' @exportS3Method marginaleffects::get_coef evinb
get_coef.evinb <- function(model, ...) coef(model, "all")

#' @exportS3Method marginaleffects::set_coef evzinb
set_coef.evzinb <- function(model, coefs, ...) evinf_set_flat_coef(model, coefs)
#' @exportS3Method marginaleffects::set_coef evinb
set_coef.evinb <- function(model, coefs, ...) evinf_set_flat_coef(model, coefs)

#' @exportS3Method marginaleffects::get_vcov evzinb
get_vcov.evzinb <- function(model, ...) vcov(model)
#' @exportS3Method marginaleffects::get_vcov evinb
get_vcov.evinb <- function(model, ...) vcov(model)

# round10 I.1 (audit §5.9): "states" (a long group/estimate frame, one group
# per component, the marginaleffects::get_predict.multinom() convention --
# group = rep(colnames, each = nrow), estimate = c(pred)), "quantile" (a
# single quantile passed through `...`, defaulting to the median so
# avg_predictions()/avg_slopes() work without it) and "exceedance" (a single
# threshold passed through `...`; several thresholds fan out into groups the
# same way "states" does) join the pre-round10 scalar types.
evinf_get_predict <- function(model, newdata, type, ...) {
  if (identical(type, "states")) {
    st <- stats::predict(model, newdata = newdata, type = "states")
    cols <- intersect(c("pr_zero", "pr_count", "pr_evi"), names(st))
    return(data.frame(
      rowid = rep(seq_len(nrow(st)), times = length(cols)),
      group = rep(sub("^pr_", "", cols), each = nrow(st)),
      estimate = unlist(st[cols], use.names = FALSE)
    ))
  }
  if (identical(type, "quantile")) {
    quantile <- list(...)$quantile %||% 0.5
    if (length(quantile) != 1L) {
      stop("get_predict(..., type = \"quantile\"): pass a single `quantile` ",
          "(marginaleffects perturbs coefficients one at a time and expects ",
          "one estimate per row); use predict() directly for several at once.",
          call. = FALSE)
    }
    # audit N1's own precedent (predict_from_boot(), used by the
    # bootstrap-based marginal_effects()): the integer mixture quantile has a
    # zero derivative almost everywhere, which starves marginaleffects'
    # numeric-differentiation delta method (NA std. errors) -- the
    # continuous (linearly-interpolated) surrogate quantiles_from_evzinb()/
    # quantiles_from_evinb() already provide via round = FALSE fixes that.
    qfn <- if (inherits(model, "evzinb")) quantiles_from_evzinb else quantiles_from_evinb
    p <- qfn(model, quantile, newdata = newdata, return_data = FALSE,
            multicore = FALSE, round = FALSE)
    return(data.frame(rowid = seq_len(NROW(p)), estimate = as.numeric(p)))
  }
  if (identical(type, "exceedance")) {
    threshold <- list(...)$threshold
    if (is.null(threshold)) {
      stop("get_predict(..., type = \"exceedance\"): `threshold` must be ",
          "provided.", call. = FALSE)
    }
    p <- stats::predict(model, newdata = newdata, type = "exceedance", threshold = threshold)
    if (length(threshold) == 1L) {
      return(data.frame(rowid = seq_len(nrow(p)), estimate = p[[1]]))
    }
    return(data.frame(
      rowid = rep(seq_len(nrow(p)), times = length(threshold)),
      group = rep(names(p), each = nrow(p)),
      estimate = unlist(p, use.names = FALSE)
    ))
  }
  ok <- c("harmonic", "explog", "counts", "pareto_alpha")
  if (length(type) != 1L || is.na(type) || !type %in% ok) {
    type <- "harmonic"
  }
  # round11 A4: clamp_alpha_pl only means anything for type = "explog"; predict()
  # itself would warn if it were forwarded for any other type (round11 A5).
  p <- if (identical(type, "explog") && !is.null(list(...)$clamp_alpha_pl)) {
    stats::predict(model, newdata = newdata, type = type,
                   clamp_alpha_pl = list(...)$clamp_alpha_pl)
  } else {
    stats::predict(model, newdata = newdata, type = type)
  }
  data.frame(rowid = seq_len(NROW(p)), estimate = as.numeric(p))
}

#' @exportS3Method marginaleffects::get_predict evzinb
get_predict.evzinb <- function(model, newdata, type = "harmonic", ...) {
  evinf_get_predict(model, newdata, type, ...)
}
#' @exportS3Method marginaleffects::get_predict evinb
get_predict.evinb <- function(model, newdata, type = "harmonic", ...) {
  evinf_get_predict(model, newdata, type, ...)
}

# insight (used by marginaleffects) tries to recover the model frame from the
# fitting environment first and warns when it cannot. The full model frame,
# with the response and every raw covariate, is stored on the object, so hand
# it straight back (issue 3.3).
#' @exportS3Method insight::get_data evzinb
get_data.evzinb <- function(x, ...) x$data$data
#' @exportS3Method insight::get_data evinb
get_data.evinb <- function(x, ...) x$data$data
