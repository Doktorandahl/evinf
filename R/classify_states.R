# audit 4.10: posterior-state tools.

# Component responsibilities (posterior state probabilities) from eq. 8 of the
# ISQ appendix: resp_k proportional to prior_k * density_k(y).
# prior is an n x 3 matrix of prior state probabilities (zero, count, evi).
#
# The evi-state density uses the discretised Pareto pmf dpareto_disc() (audit0.10
# §1.12, D.4) -- the same distribution used by the likelihood, CDF, quantiles,
# residuals and simulation -- rather than a second, independently maintained
# copy of the same formula.
#
# round9 E.1: family$count selects the count-state density (dpois() for
# "poisson", the alpha_nb -> 0 limit of the NB one); alpha_nb is unused (and
# may be any placeholder) for a Poisson count state.
# round9 E.2: family$zero == "hurdle" zero-truncates the count-state density
# for y>0 (divide by 1-f0, f0 = P(Y=0) under the untruncated count
# distribution) and forces y=0 rows to the zero state -- the count state
# cannot produce a zero at all under a hurdle, mirroring the E-step override
# in update_bfgs_fun() (src/evinf.cpp) exactly.
evinf_responsibilities <- function(y, mu_nb, alpha_nb, pl_alpha, C, prior,
                                   family = evinf_family()) {
  # round10 H.1: the three per-state densities are now the one shared
  # building block (R/evinf_pmf.R), also behind evinf_pmf()/evinf_dmix().
  d <- evinf_state_densities(y, mu_nb, alpha_nb, pl_alpha, C, prior, family = family)
  num <- cbind(d$zero, d$count, d$evi)
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
    # round10 H.1: shared with evinf_pmf()/evinf_cdf() (R/evinf_pmf.R) --
    # also clamps alpha_pl (context = "distribution"), which this newdata
    # path did not do before this refactor (predict() already clamped for
    # every other type).
    dp <- evinf_dist_params(object, newdata = newdata)
    prior <- dp$probabilities
    colnames(prior) <- c("zero", "count", "evi")
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
    posterior <- evinf_responsibilities(y, dp$nb_mu, dp$nb_alpha, dp$pl_alpha,
                                        dp$C, prior, family = dp$family)
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
