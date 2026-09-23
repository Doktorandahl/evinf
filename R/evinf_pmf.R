# The predictive distribution: one shared definition (round10 H.1, audit
# §5.8, following dev/plan_extensions_5.1_5.2_5.4.md §B.1).
#
# evinf_state_densities() is the actual shared building block: the n x 3
# (or n x K x 3, for a support matrix) prior-weighted per-state density at
# given y values. Everything else in this file, and every refactored
# consumer (mixture_p()/mixture_quantile() in R/mixture_quantile_extractor.R,
# residuals(type = "quantile"), classify_states()), is built on it or on
# mixture_p() (already the package's one CDF implementation, closed-form via
# pnbinom()/ppois()/ppareto_disc() rather than summing a pmf).
#
# `y` may be a length-n vector (one value per row) or an n x K matrix (K
# support values per row, one row per observation) -- R's usual recycling of
# a length-n vector (mu_nb, alpha_nb, pl_alpha, and each prior column) against
# an n x K matrix aligns exactly column-by-column, so every function here is
# "vectorised over rows" and "vectorised over a support grid" for free,
# without a separate code path for either.

# Per-state prior-weighted densities (zero, count, evi) at y (round9's
# evinf_responsibilities() factored into building block + normalisation).
# round9 E.1: family$count selects the count-state density (dpois() for
# "poisson"); alpha_nb is unused for a Poisson count state.
# round9 E.2: family$zero == "hurdle" zero-truncates the count-state density
# for y > 0 and forces y = 0 to the zero state entirely.
evinf_state_densities <- function(y, mu_nb, alpha_nb, pl_alpha, C, prior,
                                  family = evinf_family()) {
  d_zero <- (y == 0) * 1
  d_count <- if (family$count == "poisson") {
    stats::dpois(y, lambda = mu_nb)
  } else {
    stats::dnbinom(y, mu = mu_nb, size = 1 / alpha_nb)
  }
  d_evi <- dpareto_disc(y, C, pl_alpha)

  if (identical(family$zero, "hurdle")) {
    f0 <- if (family$count == "poisson") {
      exp(-mu_nb)
    } else {
      (1 + alpha_nb * mu_nb)^(-1 / alpha_nb)
    }
    d_count <- ifelse(y == 0, 0, d_count / (1 - f0))
  }

  list(zero = prior[, 1] * d_zero, count = prior[, 2] * d_count, evi = prior[, 3] * d_evi)
}

# The mixture pmf itself: P(Y = y) row-wise (or an n x K matrix, for a
# support matrix y) -- the sum of the three (prior-weighted) state densities.
evinf_dmix <- function(y, mu_nb, alpha_nb, pl_alpha, C, prior,
                       family = evinf_family()) {
  d <- evinf_state_densities(y, mu_nb, alpha_nb, pl_alpha, C, prior, family = family)
  d$zero + d$count + d$evi
}

# Shared parameter extraction for evinf_pmf()/evinf_cdf() (and, through them,
# every refactored consumer): the same prob_from_evzinb()/prob_from_evinb()/
# counts_from_evzinb()/fitted_alpha_from_evzinb() calls
# quantiles_from_evzinb()/quantiles_from_evinb() already made, in one place.
# alpha_pl is clamped (context = "distribution") exactly like harmonic/
# explog/quantile prediction already clamp it -- residuals(type =
# "quantile")/classify_states() did not clamp before this refactor; they now
# do too, for the same reason predict() does (a collapsed alpha_pl gives a
# degenerate discretised-Pareto tail otherwise).
evinf_dist_params <- function(object, newdata = NULL) {
  is_zinb <- inherits(object, "evzinb")
  prbs <- if (is_zinb) {
    prob_from_evzinb(object, newdata = newdata)
  } else {
    prob_from_evinb(object, newdata = newdata)
  }
  cnts <- counts_from_evzinb(object, newdata = newdata)
  alphs <- fitted_alpha_from_evzinb(object, newdata = newdata)
  alphs$pareto_alpha <- evinf_clamp_alpha_pl(
    alphs$pareto_alpha, floor = object$control$alpha_pl_floor %||% 0.01,
    context = "distribution"
  )
  list(
    nb_mu = cnts$count,
    nb_alpha = object$coef$Alpha.NB,
    pl_alpha = alphs$pareto_alpha,
    C = object$coef$C,
    probabilities = cbind(if (is_zinb) prbs$pr_zc else 0, prbs$pr_count, prbs$pr_pareto),
    family = object$family %||% evinf_family()
  )
}

# Broadcast a shared support vector into an n x K matrix, one row per
# observation -- what evinf_pmf()/evinf_cdf() hand to evinf_dmix()/
# mixture_p() for `support = `.
evinf_support_matrix <- function(support, n) {
  matrix(rep(support, each = n), nrow = n, ncol = length(support))
}

#' The predictive distribution of a fitted evzinb / evinb model
#'
#' \code{evinf_pmf()} / \code{evinf_cdf()} are the package's one definition
#' of the predictive distribution -- family-aware across all four (count x
#' zero) combinations -- that \code{predict(type = "distribution"/"quantile"/
#' "exceedance"/"draws")}, \code{residuals(type = "quantile")} and
#' \code{classify_states()} all build on.
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model.
#' @param newdata Optional new data; \code{NULL} uses the estimation data.
#' @param y Optional numeric vector, one value per row of \code{newdata} (or
#'   the estimation data): the pmf/cdf at each row's own \code{y}. Exactly
#'   one of \code{y} or \code{support} must be given.
#' @param support Optional numeric vector, the same support for every row:
#'   returns an \eqn{n \times K} matrix. Exactly one of \code{y} or
#'   \code{support} must be given.
#'
#' @return A numeric vector (with \code{y}) or an \eqn{n \times K} numeric
#'   matrix (with \code{support}).
#' @keywords internal
evinf_pmf <- function(object, newdata = NULL, y = NULL, support = NULL) {
  if (is.null(y) == is.null(support)) {
    stop("evinf_pmf(): exactly one of `y` or `support` must be given.", call. = FALSE)
  }
  p <- evinf_dist_params(object, newdata = newdata)
  yy <- y %||% evinf_support_matrix(support, length(p$nb_mu))
  evinf_dmix(yy, p$nb_mu, p$nb_alpha, p$pl_alpha, p$C, p$probabilities, family = p$family)
}

#' @rdname evinf_pmf
#' @keywords internal
evinf_cdf <- function(object, newdata = NULL, y = NULL, support = NULL) {
  if (is.null(y) == is.null(support)) {
    stop("evinf_cdf(): exactly one of `y` or `support` must be given.", call. = FALSE)
  }
  p <- evinf_dist_params(object, newdata = newdata)
  yy <- y %||% evinf_support_matrix(support, length(p$nb_mu))
  mixture_p(yy, p$pl_alpha, p$C, p$nb_mu, p$nb_alpha, p$probabilities,
           family_count = p$family$count, family_zero = p$family$zero)
}
