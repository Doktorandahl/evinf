# CDF of the EVZINB / EVINB mixture at integer x, vectorised over the rows
# (pl_alphas, nb_mu and the rows of `probabilities`). The extreme-value piece is
# the discretised Pareto (see R/dist_pareto.R). `probabilities` has columns
# (zero, count, extreme-value); the zero column is 0 for evinb.
mixture_p <- function(x, pl_alphas, C, nb_mu, nb_alpha, probabilities) {
  probabilities[, 1] +
    probabilities[, 2] * pnbinom(x, mu = nb_mu, size = 1 / nb_alpha) +
    probabilities[, 3] * ppareto_disc(x, C, pl_alphas)
}

# Smallest integer q with the mixture CDF F(q) >= p, found per row by vectorised
# bisection on mixture_p(). With `continuous = TRUE` the CDF step is linearly
# interpolated, giving a monotone real-valued surrogate for numerical
# differentiation in marginal_effects() (audit N1). `p` is a single probability.
mixture_quantile <- function(p, pl_alphas, C, nb_mu, nb_alpha, probabilities,
                             continuous = FALSE) {
  n <- length(nb_mu)
  Fp <- function(x) {
    mixture_p(x, pl_alphas, C, nb_mu, nb_alpha, probabilities)
  }
  max_q <- 1e15
  lo <- rep(0, n)
  hi <- rep(1, n)
  need <- Fp(hi) < p
  while (any(need) && max(hi) < max_q) {
    hi[need] <- pmin(hi[need] * 2, max_q)
    need <- Fp(hi) < p
  }
  # Binary search for the smallest integer in [lo, hi] whose CDF reaches p.
  while (any((hi - lo) > 1)) {
    mid <- floor((lo + hi) / 2)
    ge <- Fp(mid) >= p
    hi[ge] <- mid[ge]
    lo[!ge] <- mid[!ge]
  }
  q_hi <- ifelse(Fp(lo) >= p, lo, hi)
  q_hi[Fp(q_hi) < p] <- Inf              # p unreachable within max_q
  if (!continuous) {
    return(q_hi)
  }
  f_hi <- Fp(q_hi)
  q_lo <- q_hi - 1
  f_lo <- Fp(q_lo)
  f_lo[q_lo < 0] <- 0
  denom <- f_hi - f_lo
  frac <- ifelse(denom > 0, (p - f_lo) / denom, 0)
  frac <- pmin(pmax(frac, 0), 1)
  out <- pmax(q_lo + frac, 0)
  out[!is.finite(q_hi)] <- Inf
  out
}
