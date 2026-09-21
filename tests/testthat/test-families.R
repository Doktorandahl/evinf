# round9 E.0/E.1 (audit §5.5): evinf_family() and the Poisson count state.
# Default family = evinf_family() reproduces today's NB/mixture model exactly
# (see test-em-identity.R); these tests cover the constructor, the
# family="poisson" fit path, and every downstream consumer that had to learn
# to treat Alpha.NB as *absent* (not NA) for a Poisson count state.

test_that("evinf_family() validates its arguments and defaults to nbinom/mixture", {
  fam <- evinf_family()
  expect_s3_class(fam, "evinf_family")
  expect_identical(fam$count, "nbinom")
  expect_identical(fam$zero, "mixture")

  fam2 <- evinf_family(count = "poisson", zero = "hurdle")
  expect_identical(fam2$count, "poisson")
  expect_identical(fam2$zero, "hurdle")

  expect_error(evinf_family(count = "gaussian"), "should be one of")
  expect_error(evinf_family(zero = "truncated"), "should be one of")
})

test_that("print.evinf_family() reports both slots", {
  out <- capture.output(print(evinf_family(count = "poisson")))
  expect_true(any(grepl("poisson", out)))
  expect_true(any(grepl("mixture", out)))
})

test_that("evinf_resolve_family() accepts an object or a bare string", {
  expect_identical(evinf:::evinf_resolve_family(evinf_family(count = "poisson")),
                   evinf_family(count = "poisson"))
  expect_identical(evinf:::evinf_resolve_family("poisson"), evinf_family(count = "poisson"))
  expect_identical(evinf:::evinf_resolve_family("nbinom"), evinf_family())
  expect_error(evinf:::evinf_resolve_family(123), "single string")
})

test_that("evzinb()/evinb() accept family = as a bare string, equivalent to evinf_family(count = )", {
  m_str <- fit_evzinb_fast(family = "poisson", bootstrap = FALSE)
  m_obj <- fit_evzinb_fast(family = evinf_family(count = "poisson"), bootstrap = FALSE)
  expect_equal(unname(m_str$coef$Beta.NB), unname(m_obj$coef$Beta.NB))
  expect_null(m_str$coef$Alpha.NB)
})

test_that("evinb() rejects a hurdle zero process: it has no zero state to hurdle over", {
  expect_error(
    fit_evinb_fast(family = evinf_family(zero = "hurdle"), bootstrap = FALSE),
    "no zero state to hurdle over"
  )
})

# --- M-step correctness -----------------------------------------------------
#
# Rather than compare a Poisson fit against an NB fit pinned near alpha = 0
# (unstable close to that boundary for reasons unrelated to family="poisson",
# see NEWS/the round9 report), fix C_EV far above every observed count so the
# extreme-value state is unreachable and hard-pin the zero-inflation intercept
# so far below zero (-1000, not just -100: dpois(0, mu) itself underflows to
# an honest ~1e-51 for genevzinb2's largest fitted mu, so a merely "very
# negative" -100 pin (exp(-100) ~ 3.7e-44) can still numerically *beat* the
# count state's exact-zero emission density and steal responsibility for a
# handful of y=0 rows -- this is genuine zero-inflation-mixture behaviour, not
# a bug, but it contaminates this comparison; -1000 underflows pi_zero itself
# to an exact double 0, so the zero state can never win regardless of how
# small the count state's own density gets) that the zero state's prior mass
# underflows to an exact double 0. What's left is (to numerical precision) a
# plain Poisson GLM on the count formula, which lets the EM's Poisson M-step
# be checked directly against glm(family = poisson())'s closed-form MLE -- a
# much more robust ground truth than a second EM fit.
poisson_em_vs_glm <- function(formula_nb, data, tol = 1e-6) {
  d <- evinf:::evinf_design(formula_nb, data)
  y <- data[[all.vars(formula_nb)[1]]]
  np <- ncol(d$X) + 1L
  xo <- list(
    X.multinom.ZC = d$X, X.multinom.PL = d$X, X.NB = d$X, X.PL = d$X,
    offset.nb = rep(0, nrow(d$X))
  )
  ini <- list(
    Beta.multinom.ZC = c(-1000, rep(0, np - 1L)),
    Beta.multinom.PL = c(-1000, rep(0, np - 1L)),
    Beta.NB = rep(0, np),
    Alpha.NB = 1,  # inert placeholder: decoupled for family = poisson
    Beta.PL = rep(0, np),
    C = 1e6
  )
  ctrl <- evinf_control(c.lim = c(50, 1e6), init.C = 1e6,
                        max.diff.par = 1e-10, max.no.em.steps = 3000)
  fit <- suppressMessages(evinf:::em_fit_fixed_c(
    y, xo, ini, ctrl, fixed_zc = TRUE, family = evinf_family(count = "poisson")
  ))
  glm_fit <- glm(formula_nb, data = data, family = poisson())
  list(em = as.numeric(fit$par.mat$Beta.NB), glm = unname(coef(glm_fit)),
       converge = fit$converge)
}

test_that("the Poisson M-step matches glm(family = poisson())'s MLE to a tight tolerance", {
  data(genevzinb2, package = "evinf", envir = environment())
  cmp <- poisson_em_vs_glm(y ~ x1 + x2 + x3, genevzinb2)
  expect_true(cmp$converge)
  expect_equal(unname(cmp$em), cmp$glm, tolerance = 1e-6)
})

test_that("the Poisson M-step matches glm() on a second, independently simulated dataset", {
  set.seed(101)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.4)
  mu <- exp(0.4 + 0.5 * x1 - 0.3 * x2)
  d <- data.frame(y = rpois(n, mu), x1 = x1, x2 = x2)
  cmp <- poisson_em_vs_glm(y ~ x1 + x2, d)
  expect_true(cmp$converge)
  expect_equal(unname(cmp$em), cmp$glm, tolerance = 1e-6)
})

# --- Downstream consumers ----------------------------------------------------

fit_evzinb_poisson <- function(...) {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$y[1:8] <- d$y[1:8] + 40L  # keep a handful of rows in the extreme-value state
  suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d, family = "poisson", bootstrap = FALSE,
    verbose = FALSE, control = evinf_control(c.lim = c(50, 1000), init.C = 200),
    ...
  )))
}

fit_evinb_poisson <- function(...) {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$y[1:8] <- d$y[1:8] + 40L
  suppressMessages(suppressWarnings(evinb(
    y ~ x1 + x2 + x3, data = d, family = "poisson", bootstrap = FALSE,
    verbose = FALSE, control = evinf_control(c.lim = c(50, 1000), init.C = 200),
    ...
  )))
}

test_that("a Poisson fit has no Alpha.NB anywhere it would appear for NB (absent, not NA)", {
  m <- fit_evzinb_poisson()
  expect_null(m$coef$Alpha.NB)
  expect_null(m$par.mat$Alpha.NB)
  expect_false("alpha_nb" %in% names(coef(m, "all")))
  expect_true(is.na(summary(m)$model_statistics$Alpha_nb[["Alpha_NB"]]))

  mi <- fit_evinb_poisson()
  expect_null(mi$coef$Alpha.NB)
  expect_false("alpha_nb" %in% names(coef(mi, "all")))
})

test_that("npar drops by exactly one for a Poisson fit vs. the matching NB fit", {
  data(genevzinb2, package = "evinf", envir = environment())
  ctrl <- evinf_control(c.lim = c(50, 1000), init.C = 200)
  m_nb <- suppressMessages(evzinb(y ~ x1 + x2, data = genevzinb2, bootstrap = FALSE,
                                  verbose = FALSE, control = ctrl))
  m_pois <- suppressMessages(evzinb(y ~ x1 + x2, data = genevzinb2, family = "poisson",
                                    bootstrap = FALSE, verbose = FALSE, control = ctrl))
  expect_equal(length(m_pois$par.all), length(m_nb$par.all) - 1L)
})

test_that("glance()/print()/summary() report the family and tolerate an absent alpha", {
  m <- fit_evzinb_poisson()
  g <- glance(m)
  expect_equal(g$family, "poisson/mixture")
  expect_true(is.na(g$alpha))

  out <- capture.output(print(m))
  expect_true(any(grepl("poisson count, mixture zero", out)))

  s <- summary(m)
  expect_true(is.na(s$model_statistics$Alpha_nb[["Alpha_NB"]]))
  expect_true(any(grepl("Count component \\(Poisson\\)", capture.output(print(s)))))
})

test_that("tidy() omits alpha_nb from a Poisson fit's coefficient table", {
  m <- fit_evzinb_poisson()
  td <- tidy(m)
  expect_false("alpha_nb" %in% td$term)
})

test_that("residuals(), predict(type='quantile') and simulate() all run on a Poisson fit", {
  m <- fit_evzinb_poisson()
  n <- nrow(m$data$x.nb)

  r_resp <- residuals(m, type = "response")
  expect_length(r_resp, n)
  expect_true(all(is.finite(r_resp)))

  r_quant <- residuals(m, type = "quantile")
  expect_length(r_quant, n)

  q <- predict(m, type = "quantile", quantile = 0.75)
  expect_length(q, n)
  expect_true(all(is.finite(q) | is.infinite(q)))

  sim <- simulate(m, nsim = 3, seed = 1)
  expect_equal(dim(sim), c(n, 3))
  expect_true(all(vapply(sim, is.numeric, logical(1))))

  mi <- fit_evinb_poisson()
  n_i <- nrow(mi$data$x.nb)
  expect_length(residuals(mi, type = "quantile"), n_i)
  expect_length(predict(mi, type = "quantile", quantile = 0.75), n_i)
  sim_i <- simulate(mi, nsim = 1, seed = 1)
  expect_equal(nrow(sim_i), n_i)
})

test_that("marginaleffects get_coef()/set_coef() round-trip a Poisson fit without alpha_nb", {
  skip_if_not_installed("marginaleffects")
  m <- fit_evzinb_poisson()
  cf <- marginaleffects::get_coef(m)
  expect_false("alpha_nb" %in% names(cf))
  m2 <- marginaleffects::set_coef(m, cf)
  expect_equal(m2$coef$Beta.NB, m$coef$Beta.NB)
  expect_null(m2$coef$Alpha.NB)
})

test_that("bootstrap replicates of a Poisson fit inherit the Poisson family (round9 E.1 bug fix)", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$y[1:8] <- d$y[1:8] + 40L
  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2, data = d, family = "poisson", n_bootstraps = 3, boot_seed = 42,
    multicore = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200)
  )))
  expect_true(length(m$bootstraps) > 0)
  for (b in m$bootstraps) {
    if (inherits(b, "try-error")) next
    expect_null(b$coef$Alpha.NB)
    expect_identical(b$family$count, "poisson")
  }
})

test_that("default (nbinom/mixture) fits are numerically unaffected by family support", {
  m1 <- fit_evzinb_fast(bootstrap = FALSE)
  m2 <- fit_evzinb_fast(family = evinf_family(), bootstrap = FALSE)
  expect_identical(m1$coef$Beta.NB, m2$coef$Beta.NB)
  expect_identical(m1$coef$Alpha.NB, m2$coef$Alpha.NB)
  expect_identical(m1$family$count, "nbinom")
  expect_identical(m1$family$zero, "mixture")
})
