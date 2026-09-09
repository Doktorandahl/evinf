# audit N3 / round 5 Part B: the mixture uses the *discretised* Pareto whose pmf
# matches the likelihood: (C/y)^a - (C/(y+1))^a for integer y >= C, so
# F(y) = 1 - (C/(y+1))^a. R/dist_pareto.R is the single definition.

test_that("ppareto_disc() is the CDF of the discretised Pareto pmf", {
  C <- 12
  a <- 0.8
  K <- C + 50
  pmf <- vapply(C:K, function(y) (C / y)^a - (C / (y + 1))^a, numeric(1))
  expect_equal(
    sum(pmf),
    evinf:::ppareto_disc(K, C, a) - evinf:::ppareto_disc(C - 1, C, a),
    tolerance = 1e-12
  )
  expect_identical(evinf:::ppareto_disc(C - 1, C, a), 0)
  expect_true(evinf:::ppareto_disc(C, C, a) > 0)
})

test_that("dist_pareto d/p/q/r are mutually consistent", {
  C <- 8
  a <- 1.3
  y <- C:(C + 40)
  # pmf sums to the CDF increment
  expect_equal(evinf:::dpareto_disc(y, C, a),
               evinf:::ppareto_disc(y, C, a) - evinf:::ppareto_disc(y - 1, C, a))
  # quantile function inverts the CDF
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  q <- evinf:::qpareto_disc(p, C, a)
  expect_true(all(evinf:::ppareto_disc(q, C, a) >= p - 1e-9))
  expect_true(all(evinf:::ppareto_disc(q - 1, C, a) < p + 1e-9))
  # rng: P(Y >= y) = (C/y)^a  -> tail matches to Monte-Carlo error
  set.seed(1)
  draws <- evinf:::rpareto_disc(2e5, C, a)
  expect_true(all(draws >= C))
  expect_equal(mean(draws >= 20), (C / 20)^a, tolerance = 0.02)
})

test_that("ppareto_disc() is vectorised over the shape", {
  q <- c(20, 20, 20)
  a <- c(0.3, 0.8, 1.5)
  out <- evinf:::ppareto_disc(q, 12, a)
  expect_length(out, 3L)
  expect_equal(out, vapply(a, function(ai) evinf:::ppareto_disc(20, 12, ai),
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

test_that("mixture_quantile() inverts mixture_p() (integer and continuous)", {
  probs <- matrix(rep(c(0.15, 0.55, 0.30), each = 5), ncol = 3)
  args <- list(pl_alphas = rep(0.9, 5), C = 12, nb_mu = c(2, 4, 6, 8, 10),
               nb_alpha = rep(0.4, 5), probabilities = probs)
  for (p in c(0.2, 0.5, 0.8, 0.95)) {
    qi <- do.call(evinf:::mixture_quantile, c(list(p = p), args))
    Fq  <- do.call(evinf:::mixture_p, c(list(x = qi), args))
    Fq1 <- do.call(evinf:::mixture_p, c(list(x = qi - 1), args))
    expect_true(all(Fq >= p - 1e-9))
    expect_true(all(Fq1 < p + 1e-9))
    qc <- do.call(evinf:::mixture_quantile, c(list(p = p, continuous = TRUE), args))
    expect_true(all(qc > qi - 1 - 1e-9 & qc <= qi + 1e-9))
  }
})

test_that("residuals(type = 'quantile') stay finite / reproducible after N3", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  rq <- residuals(m, type = "quantile", seed = 1)
  expect_length(rq, nobs(m))
  expect_true(all(is.finite(rq)))
  expect_equal(rq, residuals(m, type = "quantile", seed = 1))
})
