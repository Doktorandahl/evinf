# audit 4.10: posterior-state tools.

# Discretised Pareto pmf on the integers >= C: P(Y = y) = (C/y)^a - (C/(y+1))^a.
evinf_pareto_pmf <- function(y, C, a) {
  out <- (C / y)^a - (C / (y + 1))^a
  out[y < C] <- 0
  out[!is.finite(out)] <- 0
  out
}

# Component responsibilities (posterior state probabilities) from eq. 8 of the
# ISQ appendix: resp_k proportional to prior_k * density_k(y).
# prior is an n x 3 matrix of prior state probabilities (zero, count, evi).
evinf_responsibilities <- function(y, mu_nb, alpha_nb, pl_alpha, C, prior) {
  d_zero <- as.numeric(y == 0)
  d_count <- stats::dnbinom(y, mu = mu_nb, size = 1 / alpha_nb)
  d_evi <- evinf_pareto_pmf(y, C, pl_alpha)

  num <- cbind(prior[, 1] * d_zero,
               prior[, 2] * d_count,
               prior[, 3] * d_evi)
  denom <- rowSums(num)
  denom[denom == 0] <- NA_real_
  resp <- num / denom
  colnames(resp) <- c("zero", "count", "evi")
  resp
}

.evinf_state_levels <- c("zero", "count", "evi")

evinf_map <- function(prob_mat) {
  idx <- max.col(prob_mat, ties.method = "first")
  idx[!stats::complete.cases(prob_mat)] <- NA_integer_
  factor(.evinf_state_levels[idx], levels = .evinf_state_levels)
}

evinf_threshold <- function(prob_mat, threshold) {
  mx <- apply(prob_mat, 1L, max)
  idx <- max.col(prob_mat, ties.method = "first")
  idx[mx <= threshold | !stats::complete.cases(prob_mat)] <- NA_integer_
  factor(.evinf_state_levels[idx], levels = .evinf_state_levels)
}

#' Prior and posterior state classification
#'
#' For each observation, returns the prior state probabilities (from the
#' inflation components) and, for the estimation data, the posterior state
#' probabilities (responsibilities), together with the maximum-a-posteriori state
#' under each.
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model.
#' @param rule \code{"map"} (assign to the most probable state) or
#'   \code{"threshold"} (assign only when the top probability exceeds
#'   \code{threshold}, otherwise \code{NA}).
#' @param threshold Threshold for \code{rule = "threshold"}.
#' @param newdata Optional data. Prior probabilities are always available;
#'   posterior probabilities need the response column and are \code{NA} otherwise.
#'
#' @return A tibble with \code{prior_zero} / \code{prior_count} / \code{prior_evi},
#'   \code{posterior_*}, \code{map_prior}, \code{map_posterior} (factors with
#'   levels \code{zero}, \code{count}, \code{evi}) and \code{y}.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' classify_states(model)
#' }
classify_states <- function(object, rule = c("map", "threshold"),
                            threshold = 0.5, newdata = NULL) {
  rule <- match.arg(rule)
  is_zinb <- inherits(object, "evzinb")

  # --- prior probabilities ---
  if (is.null(newdata)) {
    prior <- as.matrix(object$props)
    if (!is_zinb) {
      prior <- cbind(zero = 0, prior)
    }
    colnames(prior) <- c("zero", "count", "evi")
    y <- object$data$y
  } else {
    prb <- if (is_zinb) prob_from_evzinb(object, newdata = newdata) else
      prob_from_evinb(object, newdata = newdata)
    prior <- cbind(
      zero = if (is_zinb) prb$pr_zc else 0,
      count = prb$pr_count,
      evi = prb$pr_pareto
    )
    resp_name <- all.vars(object$formulas$formula_nb)[1]
    y <- if (resp_name %in% names(newdata)) newdata[[resp_name]] else NULL
  }

  # --- posterior probabilities ---
  if (is.null(newdata)) {
    posterior <- as.matrix(object$resp)
    if (!is_zinb) {
      posterior <- cbind(zero = 0, posterior)
    }
    colnames(posterior) <- c("zero", "count", "evi")
  } else if (!is.null(y)) {
    cnt <- counts_from_evzinb(object, newdata = newdata)$count
    alph <- fitted_alpha_from_evzinb(object, newdata = newdata)$pareto_alpha
    posterior <- evinf_responsibilities(y, cnt, object$coef$Alpha.NB, alph,
                                        object$coef$C, prior)
  } else {
    posterior <- matrix(NA_real_, nrow = nrow(prior), ncol = 3,
                        dimnames = list(NULL, c("zero", "count", "evi")))
  }

  assign_state <- function(m) {
    if (rule == "map") evinf_map(m) else evinf_threshold(m, threshold)
  }

  tibble::tibble(
    prior_zero = prior[, "zero"],
    prior_count = prior[, "count"],
    prior_evi = prior[, "evi"],
    posterior_zero = posterior[, "zero"],
    posterior_count = posterior[, "count"],
    posterior_evi = posterior[, "evi"],
    map_prior = assign_state(prior),
    map_posterior = assign_state(posterior),
    y = if (is.null(y)) NA_real_ else as.numeric(y)
  )
}

#' Cross-tabulate prior vs posterior state classification
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model.
#' @return An object of class \code{evinf_state_table}: the \code{map_prior} x
#'   \code{map_posterior} contingency table with row percentages.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' state_table(model)
#' }
state_table <- function(object) {
  cs <- classify_states(object)
  tab <- table(prior = cs$map_prior, posterior = cs$map_posterior,
               useNA = "ifany")
  pct <- prop.table(tab, margin = 1L) * 100
  structure(list(counts = tab, row_pct = pct), class = "evinf_state_table")
}

#' @export
print.evinf_state_table <- function(x, ...) {
  cat("Prior x posterior state classification (counts)\n")
  print(x$counts)
  cat("\nRow percentages\n")
  print(round(x$row_pct, 1))
  invisible(x)
}
