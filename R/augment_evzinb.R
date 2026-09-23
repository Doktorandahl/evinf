# round10 I.2 (audit §5.9): a broom-style augment() method, reusing
# classify_states() (round10 H.1's refactor of it), predict(type =
# "harmonic") and residuals(type = "quantile") (round10 I.2's newdata =
# extension of it) rather than a fourth copy of the same prior/posterior/
# prediction machinery.

# `newdata_arg` is what predict()/classify_states()/residuals() are called
# with (NULL means "the estimation data", exactly like every other method in
# this package); `nd` is the actual data frame the new columns are bound
# onto, so the two stay in sync without comparing data frames by reference.
evinf_augment_engine <- function(x, data, newdata, seed, evzinb) {
  newdata_arg <- newdata %||% data
  nd <- newdata_arg %||% x$data$data

  cs <- classify_states(x, newdata = newdata_arg)
  fitted_vals <- as.numeric(stats::predict(x, newdata = newdata_arg, type = "harmonic"))

  resp_name <- all.vars(x$formulas$formula_nb)[1]
  has_response <- resp_name %in% names(nd)

  prior_cols <- if (evzinb) {
    tibble::tibble(.prob_zero = cs$prior_zero, .prob_count = cs$prior_count,
                   .prob_evi = cs$prior_evi)
  } else {
    tibble::tibble(.prob_count = cs$prior_count, .prob_evi = cs$prior_evi)
  }

  out <- dplyr::bind_cols(
    tibble::as_tibble(nd),
    tibble::tibble(.fitted = fitted_vals),
    prior_cols,
    tibble::tibble(.state = cs$map_prior)
  )

  if (has_response) {
    resid_vals <- stats::residuals(x, newdata = newdata_arg, type = "quantile", seed = seed)
    post_cols <- if (evzinb) {
      tibble::tibble(.post_zero = cs$posterior_zero, .post_count = cs$posterior_count,
                     .post_evi = cs$posterior_evi)
    } else {
      tibble::tibble(.post_count = cs$posterior_count, .post_evi = cs$posterior_evi)
    }
    out <- dplyr::bind_cols(
      out, tibble::tibble(.resid = resid_vals), post_cols,
      tibble::tibble(.post_state = cs$map_posterior)
    )
  }

  out
}

#' Augment data with fitted values, residuals and state classification
#'
#' @param x A fitted \code{evzinb} model.
#' @param data,newdata Data to augment. \code{newdata} takes precedence over
#'   \code{data} when both are given (the \pkg{broom} convention); both
#'   default to \code{NULL}, in which case the model's own estimation data is
#'   used.
#' @param seed Optional seed for the randomized quantile residual
#'   (\code{.resid}); see \code{residuals(type = "quantile")}.
#' @param ... Unused.
#'
#' @return \code{data}/\code{newdata} (or the estimation data) with columns
#'   appended: \code{.fitted} (the \code{predict(type = "harmonic")} point
#'   prediction), \code{.prob_zero}/\code{.prob_count}/\code{.prob_evi} (prior
#'   state probabilities -- \code{augment.evinb()} has no zero state, so
#'   \code{.prob_zero} is dropped there) and \code{.state} (the MAP prior
#'   state, from \code{\link{classify_states}()}). When the response column is
#'   present, three more columns are added: \code{.resid} (a seeded
#'   randomized quantile residual, \code{residuals(type = "quantile")}),
#'   \code{.post_zero}/\code{.post_count}/\code{.post_evi} (posterior state
#'   probabilities) and \code{.post_state} (the MAP posterior state).
#' @export
#' @importFrom generics augment
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3, data=genevzinb2, n_bootstraps = 5)
#' augment(model)
#' }
augment.evzinb <- function(x, data = NULL, newdata = NULL, seed = NULL, ...) {
  evinf_augment_engine(x, data = data, newdata = newdata, seed = seed, evzinb = TRUE)
}

#' Augment data with fitted values, residuals and state classification
#'
#' @param x A fitted \code{evinb} model.
#' @inheritParams augment.evzinb
#'
#' @return See \code{\link{augment.evzinb}}; \code{evinb} has no zero state,
#'   so \code{.prob_zero}/\code{.post_zero} are not produced.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evinb(y~x1+x2+x3, data=genevzinb2, n_bootstraps = 5)
#' augment(model)
#' }
augment.evinb <- function(x, data = NULL, newdata = NULL, seed = NULL, ...) {
  evinf_augment_engine(x, data = data, newdata = newdata, seed = seed, evzinb = FALSE)
}
