# round9 E.2 (audit §5.5): family = evinf_family(zero = "hurdle") for evzinb().
# The zero state owns every zero (no mixing with the count state there); the
# count state is zero-truncated for y > 0. See src/evinf.cpp's
# hurdle_trunc_derivs_fun() for the derivative math this exercises.

# --- Finite-difference verification of the new derivatives -----------------
#
# G(mu, alpha) = -log(1 - f0), f0 = P(Y=0) under the untruncated count
# distribution, is what a hurdle zero process adds to the count state's
# log-density for y > 0. These are pure-R re-implementations of the formulas
# committed in hurdle_trunc_derivs_fun() (src/evinf.cpp), checked against
# central finite differences over a grid of mu and alpha before that C++ was
# written -- kept here as a standing regression test on the algebra itself,
# independent of the EM machinery.

logf0_nb <- function(mu, alpha) -(1 / alpha) * log(1 + alpha * mu)
G_nb <- function(mu, alpha) -log(-expm1(logf0_nb(mu, alpha)))
G_pois <- function(mu) -log(-expm1(-mu))

dGdmu_nb <- function(mu, alpha) {
  oam <- 1 + alpha * mu
  f0 <- exp(logf0_nb(mu, alpha))
  (f0 / (1 - f0)) * (-1 / oam)
}
dGdalpha_nb <- function(mu, alpha) {
  oam <- 1 + alpha * mu
  f0 <- exp(logf0_nb(mu, alpha))
  q_alpha <- log(oam) / alpha^2 - mu / (alpha * oam)
  (f0 / (1 - f0)) * q_alpha
}
d2Gdmu2_nb <- function(mu, alpha) {
  oam <- 1 + alpha * mu
  f0 <- exp(logf0_nb(mu, alpha))
  p <- f0 / (1 - f0)
  hpp <- p * (1 + p)
  q_mu <- -1 / oam
  Q_mumu <- alpha / oam^2
  hpp * q_mu^2 + p * Q_mumu
}
d2Gdmudalpha_nb <- function(mu, alpha) {
  oam <- 1 + alpha * mu
  f0 <- exp(logf0_nb(mu, alpha))
  p <- f0 / (1 - f0)
  hpp <- p * (1 + p)
  q_mu <- -1 / oam
  q_alpha <- log(oam) / alpha^2 - mu / (alpha * oam)
  Q_mualpha <- mu / oam^2
  hpp * q_mu * q_alpha + p * Q_mualpha
}
d2Gdalpha2_nb <- function(mu, alpha) {
  oam <- 1 + alpha * mu
  f0 <- exp(logf0_nb(mu, alpha))
  p <- f0 / (1 - f0)
  hpp <- p * (1 + p)
  q_alpha <- log(oam) / alpha^2 - mu / (alpha * oam)
  Q_alphaalpha <- mu * (2 + 3 * alpha * mu) / (alpha^2 * oam^2) - 2 * log(oam) / alpha^3
  hpp * q_alpha^2 + p * Q_alphaalpha
}
dGdmu_pois <- function(mu) {
  f0 <- exp(-mu)
  -(f0 / (1 - f0))
}
d2Gdmu2_pois <- function(mu) {
  f0 <- exp(-mu)
  p <- f0 / (1 - f0)
  p * (1 + p)
}

fd1 <- function(f, x, h = 1e-6) (f(x + h) - f(x - h)) / (2 * h)
fd2 <- function(f, x, h = 1e-4) (f(x + h) - 2 * f(x) + f(x - h)) / h^2
fd2_cross <- function(f, x, y, hx = 1e-4, hy = 1e-4) {
  (f(x + hx, y + hy) - f(x + hx, y - hy) - f(x - hx, y + hy) + f(x - hx, y - hy)) / (4 * hx * hy)
}

# round10 0.9 (review §7): a plain central second difference balances O(h^2)
# truncation error against O(1/h^2) roundoff amplification, so no single h
# gets much past ~1e-2 relative accuracy. Richardson extrapolation --
# combining the estimates at h and h/2 as (4*D(h/2) - D(h))/3 -- cancels the
# leading h^2 error term, giving O(h^4) accuracy without shrinking h enough
# to hit the roundoff floor, and reaches ~1e-6 (verified against the raw
# central difference above, which this replaces rather than supplements).
fd2_richardson <- function(f, x, h = 1e-2) {
  (4 * fd2(f, x, h / 2) - fd2(f, x, h)) / 3
}
fd2_cross_richardson <- function(f, x, y, hx = 1e-2, hy = 1e-2) {
  (4 * fd2_cross(f, x, y, hx / 2, hy / 2) - fd2_cross(f, x, y, hx, hy)) / 3
}

test_that("hurdle zero-truncation derivatives (NB) match central finite differences", {
  mus <- c(0.5, 2, 10, 50)
  alphas <- c(0.05, 0.3, 1, 3)
  for (mu in mus) for (alpha in alphas) {
    expect_equal(dGdmu_nb(mu, alpha), fd1(function(m) G_nb(m, alpha), mu), tolerance = 1e-4)
    expect_equal(dGdalpha_nb(mu, alpha), fd1(function(a) G_nb(mu, a), alpha), tolerance = 1e-4)
    expect_equal(d2Gdmu2_nb(mu, alpha), fd2_richardson(function(m) G_nb(m, alpha), mu), tolerance = 1e-6)
    expect_equal(d2Gdmudalpha_nb(mu, alpha), fd2_cross_richardson(G_nb, mu, alpha), tolerance = 1e-6)
    expect_equal(d2Gdalpha2_nb(mu, alpha), fd2_richardson(function(a) G_nb(mu, a), alpha), tolerance = 1e-6)
  }
})

test_that("hurdle zero-truncation derivatives (Poisson) match central finite differences", {
  for (mu in c(0.5, 2, 10, 50)) {
    expect_equal(dGdmu_pois(mu), fd1(G_pois, mu), tolerance = 1e-4)
    expect_equal(d2Gdmu2_pois(mu), fd2_richardson(G_pois, mu), tolerance = 1e-6)
  }
})

# --- End-to-end verification against pscl::hurdle() ------------------------
#
# Pin C_EV far above every observed count (the extreme-value state becomes
# unreachable) so the fitted model reduces exactly to a two-part hurdle model
# -- pscl::hurdle()'s model -- letting the full EM pipeline (E-step forced
# assignment, M-step gradient/Hessian corrections, and log_lik_fun()'s hurdle
# branch together) be checked against an independent implementation. pscl's
# zero-hurdle component models logit(P(Y>0)) (the complementary event to this
# package's logit(P(zero state))), so its coefficients are the exact negation
# of this package's zero-inflation coefficients -- not a discrepancy.
hurdle_em_fixed_c <- function(y, x, family_count, alpha_init = 0.5) {
  dd <- evinf:::evinf_design(y ~ x, data.frame(y = y, x = x))
  np <- ncol(dd$X) + 1L
  xo <- list(X.multinom.ZC = dd$X, X.multinom.PL = dd$X, X.NB = dd$X, X.PL = dd$X,
             offset.nb = rep(0, nrow(dd$X)))
  ini <- list(
    Beta.multinom.ZC = rep(0, np),
    Beta.multinom.PL = c(-1000, rep(0, np - 1L)),
    Beta.NB = rep(0, np),
    Alpha.NB = alpha_init,
    Beta.PL = rep(0, np),
    C = 1e6
  )
  ctrl <- evinf_control(c.lim = c(50, 1e6), init.C = 1e6, max.diff.par = 1e-9,
                        max.no.em.steps = 4000)
  suppressMessages(evinf:::em_fit_fixed_c(
    y, xo, ini, ctrl, fixed_zc = FALSE,
    family = evinf_family(count = family_count, zero = "hurdle")
  ))
}

test_that("hurdle NB matches pscl::hurdle(dist = 'negbin') to numerical precision", {
  skip_if_not_installed("pscl")
  skip_on_cran()
  set.seed(42)
  n <- 500
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  zc_eta <- 0.3 - 0.4 * x1 + 0.2 * x2
  mu <- exp(1.0 + 0.4 * x1 - 0.3 * x2)
  alpha_true <- 0.7
  rztnb <- function(n, mu, alpha) {
    out <- integer(n)
    for (i in seq_len(n)) repeat {
      v <- rnbinom(1, mu = mu[i], size = 1 / alpha)
      if (v > 0) { out[i] <- v; break }
    }
    out
  }
  is_zero <- runif(n) < plogis(zc_eta)
  y <- integer(n)
  y[!is_zero] <- rztnb(sum(!is_zero), mu[!is_zero], alpha_true)
  d <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- hurdle_em_fixed_c(y, cbind(x1, x2), family_count = "nbinom")
  expect_true(fit$converge)

  h <- pscl::hurdle(y ~ x1 + x2 | x1 + x2, data = d, dist = "negbin",
                    zero.dist = "binomial")

  expect_equal(as.numeric(fit$par.mat$Beta.NB), unname(coef(h, model = "count")),
               tolerance = 1e-3)
  expect_equal(unname(fit$par.mat$Alpha.NB), unname(1 / h$theta), tolerance = 1e-3)
  expect_equal(as.numeric(fit$par.mat$Beta.multinom.ZC),
               -unname(coef(h, model = "zero")), tolerance = 1e-3)

  ll <- evinf:::log_lik_fun(
    fit$par.mat$Beta.multinom.ZC, fit$par.mat$Beta.multinom.PL, fit$par.mat$Beta.NB,
    fit$par.mat$Alpha.NB, fit$par.mat$Beta.PL, fit$par.mat$C,
    cbind(1, x1, x2), cbind(1, x1, x2), cbind(1, x1, x2), cbind(1, x1, x2),
    y, rep(0, n), rep(0, n), rep(0, n), rep(1, n), family_count = 0L, family_zero = 1L
  )
  expect_equal(ll, as.numeric(logLik(h)), tolerance = 1e-4)
})

test_that("hurdle Poisson matches pscl::hurdle(dist = 'poisson') to numerical precision", {
  skip_if_not_installed("pscl")
  skip_on_cran()
  set.seed(7)
  n <- 500
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  zc_eta <- 0.2 - 0.5 * x1 + 0.3 * x2
  mu <- exp(0.8 + 0.3 * x1 - 0.2 * x2)
  rztpois <- function(n, mu) {
    out <- integer(n)
    for (i in seq_len(n)) repeat {
      v <- rpois(1, mu[i])
      if (v > 0) { out[i] <- v; break }
    }
    out
  }
  is_zero <- runif(n) < plogis(zc_eta)
  y <- integer(n)
  y[!is_zero] <- rztpois(sum(!is_zero), mu[!is_zero])
  d <- data.frame(y = y, x1 = x1, x2 = x2)

  fit <- hurdle_em_fixed_c(y, cbind(x1, x2), family_count = "poisson")
  expect_true(fit$converge)

  h <- pscl::hurdle(y ~ x1 + x2 | x1 + x2, data = d, dist = "poisson",
                    zero.dist = "binomial")

  expect_equal(as.numeric(fit$par.mat$Beta.NB), unname(coef(h, model = "count")),
               tolerance = 1e-3)
  expect_equal(as.numeric(fit$par.mat$Beta.multinom.ZC),
               -unname(coef(h, model = "zero")), tolerance = 1e-3)

  ll <- evinf:::log_lik_fun(
    fit$par.mat$Beta.multinom.ZC, fit$par.mat$Beta.multinom.PL, fit$par.mat$Beta.NB,
    fit$par.mat$Alpha.NB, fit$par.mat$Beta.PL, fit$par.mat$C,
    cbind(1, x1, x2), cbind(1, x1, x2), cbind(1, x1, x2), cbind(1, x1, x2),
    y, rep(0, n), rep(0, n), rep(0, n), rep(1, n), family_count = 1L, family_zero = 1L
  )
  expect_equal(ll, as.numeric(logLik(h)), tolerance = 1e-4)
})

# --- evinb() guard -----------------------------------------------------------

test_that("evinb() rejects a hurdle zero process unconditionally", {
  expect_error(
    fit_evinb_fast(family = evinf_family(zero = "hurdle"), bootstrap = FALSE),
    "no zero state to hurdle over"
  )
})

# --- Structural / S3 coverage on a full evzinb(family = hurdle) fit --------

fit_evzinb_hurdle <- function(...) {
  set.seed(11)
  n <- 300
  x1 <- rnorm(n)
  zc_eta <- 0.2 - 0.4 * x1
  mu <- exp(0.9 + 0.3 * x1)
  alpha_true <- 0.6
  rztnb <- function(n, mu, alpha) {
    out <- integer(n)
    for (i in seq_len(n)) repeat {
      v <- rnbinom(1, mu = mu[i], size = 1 / alpha)
      if (v > 0) { out[i] <- v; break }
    }
    out
  }
  is_zero <- runif(n) < plogis(zc_eta)
  y <- integer(n)
  y[!is_zero] <- rztnb(sum(!is_zero), mu[!is_zero], alpha_true)
  y[1:6] <- y[1:6] + 40L  # a handful of extreme-value-state rows
  d <- data.frame(y = y, x1 = x1)
  suppressMessages(suppressWarnings(evzinb(
    y ~ x1, data = d, family = evinf_family(zero = "hurdle"), bootstrap = FALSE,
    verbose = FALSE, control = evinf_control(c.lim = c(20, 1000), init.C = 30),
    ...
  )))
}

test_that("a hurdle fit assigns every y=0 row to the zero state exactly", {
  m <- fit_evzinb_hurdle()
  cs <- classify_states(m)
  zero_rows <- m$data$y == 0
  expect_true(all(cs$posterior_zero[zero_rows] == 1))
  expect_true(all(cs$posterior_count[zero_rows] == 0))
  expect_true(all(cs$map_posterior[zero_rows] == "zero"))
})

test_that("Alpha.NB is present (unlike a Poisson count state) for a hurdle NB fit", {
  m <- fit_evzinb_hurdle()
  expect_false(is.null(m$coef$Alpha.NB))
  expect_true("alpha_nb" %in% names(coef(m, "all")))
})

test_that("glance()/print()/summary() report a hurdle fit's family", {
  m <- fit_evzinb_hurdle()
  expect_equal(glance(m)$family, "nbinom/hurdle")
  out <- capture.output(print(m))
  expect_true(any(grepl("nbinom count, hurdle zero", out)))
})

test_that("residuals(), predict(type='quantile') and simulate() run on a hurdle fit", {
  m <- fit_evzinb_hurdle()
  n <- nrow(m$data$x.nb)
  expect_length(residuals(m, type = "quantile"), n)
  q <- predict(m, type = "quantile", quantile = 0.75)
  expect_length(q, n)
  sim <- simulate(m, nsim = 3, seed = 1)
  expect_equal(dim(sim), c(n, 3))
  # round9 E.2: a hurdle count state never draws 0.
  prbs <- prob_from_evzinb(m)
  count_state_rows <- which(prbs$pr_count > 0.99)
  if (length(count_state_rows) > 0) {
    expect_true(all(sim[count_state_rows, 1] != 0 | m$data$y[count_state_rows] == 0))
  }
})

test_that("default (mixture) fits are unaffected by hurdle support", {
  m1 <- fit_evzinb_fast(bootstrap = FALSE)
  m2 <- fit_evzinb_fast(family = evinf_family(), bootstrap = FALSE)
  expect_identical(m1$coef$Beta.NB, m2$coef$Beta.NB)
  expect_identical(m1$family$zero, "mixture")
})

test_that("on data simulated without a separate zero-inflation process, mixture and hurdle fits give close log-likelihoods and the mixture fit's zero-state prior is near 0", {
  skip_on_cran()
  set.seed(99)
  n <- 400
  x1 <- rnorm(n)
  mu <- exp(1.2 + 0.3 * x1)
  y <- rnbinom(n, mu = mu, size = 1 / 0.5)  # ordinary NB, zeros arise only from the count process
  d <- data.frame(y = y, x1 = x1)

  m_mix <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1, data = d, bootstrap = FALSE, verbose = FALSE
  )))
  m_hurdle <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1, data = d, family = evinf_family(zero = "hurdle"), bootstrap = FALSE,
    verbose = FALSE
  )))

  expect_equal(m_mix$log.lik, m_hurdle$log.lik, tolerance = 0.05 * abs(m_hurdle$log.lik))
  expect_true(mean(m_mix$props[, 1]) < 0.1)
})

test_that("an intercept-only hurdle zero probability equals the observed zero share at tight convergence (round10 0.9)", {
  # review §8: at exact convergence, an intercept-only zero logit's score
  # equation is sum(r_i0 - pi_i0) = 0; a hurdle's zero-state responsibility
  # r_i0 is exactly 1 for y = 0 rows (forced in the E-step) and 0 otherwise,
  # so mean(pi_i0) must equal the observed zero share exactly. The default
  # max.diff.par = 0.01 stops early enough that this only holds to ~1e-4
  # (review's own reproduction: 0.392159 fitted vs 0.3925 observed); a tight
  # max.diff.par makes it a free, exact test of the hurdle E/M-step algebra.
  # fit_evzinb_hurdle() already hard-codes `control =`, so this can't reuse
  # it via ... (duplicate named argument) -- inline the same data generator.
  set.seed(11)
  n <- 300
  x1 <- rnorm(n)
  zc_eta <- 0.2 - 0.4 * x1
  mu <- exp(0.9 + 0.3 * x1)
  alpha_true <- 0.6
  rztnb <- function(n, mu, alpha) {
    out <- integer(n)
    for (i in seq_len(n)) repeat {
      v <- rnbinom(1, mu = mu[i], size = 1 / alpha)
      if (v > 0) { out[i] <- v; break }
    }
    out
  }
  is_zero <- runif(n) < plogis(zc_eta)
  y <- integer(n)
  y[!is_zero] <- rztnb(sum(!is_zero), mu[!is_zero], alpha_true)
  y[1:6] <- y[1:6] + 40L
  d <- data.frame(y = y, x1 = x1)

  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1, formula_zi = ~1, data = d, family = evinf_family(zero = "hurdle"),
    bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(20, 1000), init.C = 30,
                            max.diff.par = 1e-8, max.no.em.steps = 5000)
  )))
  expect_true(m$converge)
  observed_zero_share <- mean(m$data$y == 0)
  expect_equal(mean(m$props[, "zero"]), observed_zero_share, tolerance = 1e-6)
})
