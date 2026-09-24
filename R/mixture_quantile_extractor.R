# CDF of the EVZINB / EVINB mixture at integer x, vectorised over the rows
# (pl_alphas, nb_mu and the rows of `probabilities`). The extreme-value piece is
# the discretised Pareto (see R/dist_pareto.R). `probabilities` has columns
# (zero, count, extreme-value); the zero column is 0 for evinb.
# round9 E.1: family_count = "poisson" uses ppois() (nb_alpha is unused --
# and NULL -- in that case, exactly like object$coef$Alpha.NB).
# round9 E.2: family_zero = "hurdle" zero-truncates the count state's CDF
# contribution: F_trunc(x) = (F(x) - f0) / (1 - f0), which is exactly 0 at
# x = 0 (the count state has no mass there under a hurdle) -- pmax(.,0)
# only guards floating-point noise at that boundary, not a real negative.
mixture_p <- function(x, pl_alphas, C, nb_mu, nb_alpha, probabilities,
                      family_count = c("nbinom", "poisson"),
                      family_zero = c("mixture", "hurdle")) {
  family_count <- match.arg(family_count)
  family_zero <- match.arg(family_zero)
  count_cdf <- if (family_count == "poisson") {
    ppois(x, lambda = nb_mu)
  } else {
    pnbinom(x, mu = nb_mu, size = 1 / nb_alpha)
  }
  if (family_zero == "hurdle") {
    f0 <- if (family_count == "poisson") {
      exp(-nb_mu)
    } else {
      (1 + nb_alpha * nb_mu)^(-1 / nb_alpha)
    }
    count_cdf <- pmax((count_cdf - f0) / (1 - f0), 0)
  }
  probabilities[, 1] +
    probabilities[, 2] * count_cdf +
    probabilities[, 3] * ppareto_disc(x, C, pl_alphas)
}

# Smallest integer q with the mixture CDF F(q) >= p, found per row by vectorised
# bisection on mixture_p(). With `continuous = TRUE` the CDF step is linearly
# interpolated, giving a monotone real-valued surrogate for numerical
# differentiation in marginal_effects() (audit N1).
#
# round10 H.3: `p` may be a vector of several probabilities, sharing one
# bisection over an n x length(p) grid rather than one call per probability.
# The broadcast is the same trick evinf_support_matrix() already uses
# (R/evinf_pmf.R) in the other direction: pl_alphas/nb_mu/nb_alpha/
# probabilities stay length-n vectors and recycle column-by-column against
# an n x length(p) matrix (lo/hi/mid/q_hi), exactly like they already do
# against evinf_pmf()'s n x K support matrix -- mixture_p() needed no change
# to accept a matrix `x`. For length(p) == 1, lo/hi/mid are n x 1 matrices
# instead of length-n vectors, but every floating-point operation on them is
# identical elementwise to the pre-H.3 vector code; drop() at the end
# restores the original plain-vector return so this case is byte-identical
# to before.
mixture_quantile <- function(p, pl_alphas, C, nb_mu, nb_alpha, probabilities,
                             continuous = FALSE,
                             family_count = c("nbinom", "poisson"),
                             family_zero = c("mixture", "hurdle")) {
  family_count <- match.arg(family_count)
  family_zero <- match.arg(family_zero)
  n <- length(nb_mu)
  m <- length(p)
  Fp <- function(x) {
    mixture_p(x, pl_alphas, C, nb_mu, nb_alpha, probabilities,
              family_count = family_count, family_zero = family_zero)
  }
  p_mat <- matrix(rep(p, each = n), nrow = n, ncol = m)
  max_q <- 1e15
  lo <- matrix(0, nrow = n, ncol = m)
  hi <- matrix(1, nrow = n, ncol = m)
  need <- Fp(hi) < p_mat
  while (any(need) && max(hi) < max_q) {
    hi[need] <- pmin(hi[need] * 2, max_q)
    need <- Fp(hi) < p_mat
  }
  # Binary search for the smallest integer in [lo, hi] whose CDF reaches p.
  while (any((hi - lo) > 1)) {
    mid <- floor((lo + hi) / 2)
    ge <- Fp(mid) >= p_mat
    hi[ge] <- mid[ge]
    lo[!ge] <- mid[!ge]
  }
  q_hi <- ifelse(Fp(lo) >= p_mat, lo, hi)
  q_hi[Fp(q_hi) < p_mat] <- Inf          # p unreachable within max_q
  # Only a single probability collapses back to a plain vector (byte-identical
  # to the pre-H.3 return); several probabilities stay an n x length(p) matrix
  # even when n == 1 (a single-row newdata), so column j is always "quantile j".
  as_out <- if (m == 1) function(x) drop(x) else function(x) x
  if (!continuous) {
    return(as_out(q_hi))
  }
  f_hi <- Fp(q_hi)
  q_lo <- q_hi - 1
  f_lo <- Fp(q_lo)
  f_lo[q_lo < 0] <- 0
  denom <- f_hi - f_lo
  frac <- ifelse(denom > 0, (p_mat - f_lo) / denom, 0)
  frac <- pmin(pmax(frac, 0), 1)
  out <- pmax(q_lo + frac, 0)
  out[!is.finite(q_hi)] <- Inf
  as_out(out)
}
