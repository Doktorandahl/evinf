# audit N3: mixture_p() / ppareto_vec() must use the *discretised* Pareto CDF
# that matches the pmf the likelihood uses: (C/y)^a - (C/(y+1))^a for integer
# y >= C, so F(y) = 1 - (C/(y+1))^a.

test_that("ppareto_vec() is the CDF of the discretised Pareto pmf", {
  C <- 12
  a <- 0.8
  K <- C + 50
  pmf <- vapply(C:K, function(y) (C / y)^a - (C / (y + 1))^a, numeric(1))
  # sum of pmf from C..K == F(K) - F(C-1)
  expect_equal(
    sum(pmf),
    evinf:::ppareto_vec(K, C, a) - evinf:::ppareto_vec(C - 1, C, a),
    tolerance = 1e-12
  )
  expect_identical(evinf:::ppareto_vec(C - 1, C, a), 0)
  expect_true(evinf:::ppareto_vec(C, C, a) > 0)
})

test_that("ppareto_vec() is vectorised over the shape", {
  q <- c(20, 20, 20)
  a <- c(0.3, 0.8, 1.5)
  out <- evinf:::ppareto_vec(q, 12, a)
  expect_length(out, 3L)
  expect_equal(out, vapply(a, function(ai) evinf:::ppareto_vec(20, 12, ai),
                           numeric(1)))
})

test_that("mixture_p() is a proper (non-decreasing, -> 1) CDF", {
  probs <- matrix(rep(c(0.2, 0.5, 0.3), each = 6), ncol = 3)
  x <- 0:40
  p <- vapply(x, function(xi) {
    evinf:::mixture_p(rep(xi, 6), pl_alphas = rep(0.9, 6), C = 12,
                      nb_mu = rep(4, 6), nb_alpha = rep(0.5, 6),
                      probabilities = probs)[1]
  }, numeric(1))
  expect_true(all(diff(p) >= -1e-12))
  expect_gt(evinf:::mixture_p(rep(1e9, 6), rep(0.9, 6), 12, rep(4, 6),
                              rep(0.5, 6), probs)[1], 0.999)
})

test_that("residuals(type = 'quantile') stay finite / reproducible after N3", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  rq <- residuals(m, type = "quantile", seed = 1)
  expect_length(rq, nobs(m))
  expect_true(all(is.finite(rq)))
  expect_equal(rq, residuals(m, type = "quantile", seed = 1))
})
