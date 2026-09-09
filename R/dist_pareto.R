# Discretised Pareto distribution for the extreme-value component of the
# EVZINB / EVINB mixture (round 5 Part B).
#
# For an integer outcome y, real scale C > 0 and shape alpha > 0 this is the
# floor() of the continuous Pareto with minimum C and shape alpha, so that
# P(Y >= y) = (C / y)^alpha exactly and
#
#   pmf  f(y) = (C / y)^alpha - (C / (y + 1))^alpha   for y >= C, else 0
#   cdf  F(y) = 1 - (C / (y + 1))^alpha               for y >= C, else 0
#   qf   Q(p) = smallest integer y with F(y) >= p
#   rng  floor(C * U^(-1/alpha)),  U ~ Uniform(0, 1)
#
# This is the single definition the whole package uses (the likelihood, the
# mixture CDF, residuals, prediction quantiles and simulation). Every argument
# is vectorised and recycled (mistr's Pareto helpers are not vectorised over
# `shape`, which is why the package no longer depends on mistr).

#' @keywords internal
ppareto_disc <- function(q, C, alpha) {
  out <- 1 - (C / (q + 1))^alpha
  out[q < C] <- 0
  out[!is.finite(out)] <- 0
  out[out < 0] <- 0
  out
}

#' @keywords internal
dpareto_disc <- function(x, C, alpha) {
  ppareto_disc(x, C, alpha) - ppareto_disc(x - 1, C, alpha)
}

#' @keywords internal
qpareto_disc <- function(p, C, alpha) {
  # F(y) >= p  <=>  (C / (y + 1))^alpha <= 1 - p  <=>  y >= C / (1 - p)^(1/alpha) - 1
  y <- ceiling(C / (1 - p)^(1 / alpha) - 1)
  m <- ceiling(C)                       # support minimum
  y <- pmax(y, m)
  y[p <= 0] <- m
  y[p >= 1] <- Inf
  y
}

#' @keywords internal
rpareto_disc <- function(n, C, alpha) {
  floor(C * stats::runif(n)^(-1 / alpha))
}
