# round10 J.3 (audit §5.11, decision D8): a benchmark data generator rather
# than a bundled data file, so inst/bench/bench_evinf.R (and anyone doing
# their own scale/stress testing) can ask for whatever n they need -- a
# 100-row and a 20000-row benchmark case, say -- without the installed
# package carrying an extra .rda per size.

#' Synthetic data from an EVZINB model with known true parameters
#'
#' Generates \code{n} observations from an EVZINB mixture with fixed,
#' hand-chosen true parameters (see Details) -- \code{formula_nb = y ~ x1 +
#' x2 + x3}, \code{formula_zi}/\code{formula_evi} sharing the same three
#' covariates, \code{formula_pareto} depending on \code{x1} alone. Meant for
#' benchmarking (\code{inst/bench/bench_evinf.R}) and scale testing, not as a
#' realistic applied example -- see \code{\link{genevzinb2}} for that.
#'
#' @param n Number of observations.
#' @param seed Optional seed; the caller's \code{.Random.seed} is left
#'   untouched either way (the same convention as \code{simulate()}/
#'   \code{predict(type = "draws")}).
#'
#' @details True parameters (on each component's own scale, the same
#'   parameterisation \code{coef()} returns):
#'   \itemize{
#'     \item Zero-state log-odds against the count-state baseline:
#'       \eqn{-0.5 + 0.8 x_1 - 0.5 x_2}.
#'     \item Evi-state log-odds against the count-state baseline:
#'       \eqn{-2.5 + 0.5 x_1 + 0.3 x_3} (a deliberately small share, matching
#'       an extreme-value-inflation process being the rare state).
#'     \item Count state (negative binomial, log link):
#'       \eqn{\mu = \exp(1.5 + 0.4 x_1 - 0.3 x_2 + 0.2 x_3)}, dispersion
#'       \eqn{\alpha_{NB} = 0.8}.
#'     \item Evi state (the discretised Pareto tail the package fits
#'       throughout): threshold \eqn{C = 50}, shape \eqn{\alpha_{PL} =
#'       \exp(0.6 + 0.15 x_1)} (kept well above 0 for any simulated
#'       \code{x1}, so the tail never collapses).
#'   }
#'
#' @return A tibble with columns \code{y}, \code{x1}, \code{x2}, \code{x3}.
#' @export
#'
#' @examples
#' d <- evinf_bench_data(200, seed = 1)
#' \donttest{
#' m <- evzinb(y ~ x1 + x2 + x3, data = d, bootstrap = FALSE)
#' }
evinf_bench_data <- function(n = 1000, seed = NULL) {
  evinf_with_seed(seed, {
    x1 <- stats::rnorm(n)
    x2 <- stats::rnorm(n)
    x3 <- stats::rnorm(n)

    eta_zero <- -0.5 + 0.8 * x1 - 0.5 * x2
    eta_evi <- -2.5 + 0.5 * x1 + 0.3 * x3
    sp <- evinf_stable_props3(eta_zero, eta_evi)

    mu <- exp(1.5 + 0.4 * x1 - 0.3 * x2 + 0.2 * x3)
    alpha_nb <- 0.8
    C <- 50
    alpha_pl <- exp(0.6 + 0.15 * x1)

    state_draw <- stats::runif(n)
    count_draws <- stats::rnbinom(n, mu = mu, size = 1 / alpha_nb)
    pl_draws <- rpareto_disc(n, C, alpha_pl)

    y <- dplyr::case_when(
      state_draw <= sp[, "zero"] ~ 0,
      state_draw <= sp[, "zero"] + sp[, "count"] ~ count_draws,
      TRUE ~ pl_draws
    )

    tibble::tibble(y = as.numeric(y), x1 = x1, x2 = x2, x3 = x3)
  })
}
