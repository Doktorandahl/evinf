# round9 D.2 (audit §5.6): weights = for evzinb() / evinb(). Frequency-weight
# semantics: every observation's contribution to the log-likelihood and the
# M-step accumulations is multiplied by its weight; nobs() (and therefore
# AIC/BIC/the approximate t p-values) use sum(weights).

tight_ctrl_w <- function(max.diff.par = 1e-6, max.no.em.steps = 3000, ...) {
  evinf_control(c.lim = c(50, 1000), init.C = 200,
                max.diff.par = max.diff.par, max.no.em.steps = max.no.em.steps, ...)
}

test_that("a fit with integer weights equals a fit on the row-duplicated data", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(5)
  w_int <- sample(1:3, nrow(genevzinb2), replace = TRUE)
  d <- genevzinb2
  d$w <- w_int

  m_w <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d, weights = w,
    bootstrap = FALSE, verbose = FALSE, control = tight_ctrl_w()
  )))

  idx <- rep(seq_len(nrow(d)), w_int)
  d_dup <- d[idx, ]
  m_dup <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d_dup,
    bootstrap = FALSE, verbose = FALSE, control = tight_ctrl_w()
  )))

  expect_equal(unname(m_w$coef$Beta.NB), unname(m_dup$coef$Beta.NB), tolerance = 1e-8)
  expect_equal(unname(m_w$coef$Beta.multinom.ZC), unname(m_dup$coef$Beta.multinom.ZC),
               tolerance = 1e-8)
  expect_equal(unname(m_w$coef$Beta.multinom.PL), unname(m_dup$coef$Beta.multinom.PL),
               tolerance = 1e-8)
  expect_equal(unname(m_w$coef$Beta.PL), unname(m_dup$coef$Beta.PL), tolerance = 1e-8)
  expect_equal(m_w$coef$Alpha.NB, m_dup$coef$Alpha.NB, tolerance = 1e-8)
  expect_equal(m_w$log.lik, m_dup$log.lik, tolerance = 1e-8)
})

test_that("evinb() with integer weights equals a fit on the row-duplicated data", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(6)
  w_int <- sample(1:2, nrow(genevzinb2), replace = TRUE)
  d <- genevzinb2
  d$w <- w_int
  ctrl <- tight_ctrl_w(max.diff.par = 1e-8, max.no.em.steps = 5000)

  m_w <- suppressMessages(suppressWarnings(evinb(
    y ~ x1 + x2 + x3, data = d, weights = w,
    bootstrap = FALSE, verbose = FALSE, control = ctrl
  )))

  idx <- rep(seq_len(nrow(d)), w_int)
  d_dup <- d[idx, ]
  m_dup <- suppressMessages(suppressWarnings(evinb(
    y ~ x1 + x2 + x3, data = d_dup,
    bootstrap = FALSE, verbose = FALSE, control = ctrl
  )))

  # evinb's fixed-ZC block reaches max.diff.par's stopping criterion from a
  # slightly different iteration path for the weighted vs. row-duplicated
  # data (same fixed point, different floating-point accumulation order
  # along the way -- same phenomenon as A.3's trajectory sensitivity,
  # review §1), so this needs a looser tolerance than the direct
  # log_lik_fun() row-duplication identity above (which matches to 1e-13).
  expect_equal(unname(m_w$coef$Beta.NB), unname(m_dup$coef$Beta.NB), tolerance = 1e-6)
  expect_equal(m_w$log.lik, m_dup$log.lik, tolerance = 1e-6)
})

test_that("nobs() is sum(weights), and AIC/BIC match a hand computation", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(7)
  d <- genevzinb2
  d$w <- runif(nrow(d), 0.5, 3)  # analytic (non-integer) weights

  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d, weights = w,
    bootstrap = FALSE, verbose = FALSE, control = tight_ctrl_w()
  )))

  expect_equal(nobs(m), sum(d$w))
  expect_equal(glance(m)$sum_weights, sum(d$w))
  expect_equal(glance(m)$nobs, nrow(d))

  npar <- length(m$par.all)
  expect_equal(m$AIC, 2 * npar - 2 * m$log.lik)
  expect_equal(m$BIC, log(sum(d$w)) * npar - 2 * m$log.lik)
})

test_that("weights must be positive and finite", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$w0 <- 0
  d$wneg <- -1
  d$wna <- NA_real_
  d$winf <- Inf

  expect_error(
    evzinb(y ~ x1, data = d, weights = w0, bootstrap = FALSE, verbose = FALSE),
    "positive and finite"
  )
  expect_error(
    evzinb(y ~ x1, data = d, weights = wneg, bootstrap = FALSE, verbose = FALSE),
    "positive and finite"
  )
  expect_error(
    evzinb(y ~ x1, data = d, weights = winf, bootstrap = FALSE, verbose = FALSE),
    "positive and finite"
  )
})

test_that("weights accepts a bare column name, a string, and a numeric vector", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$wt <- rep(1, nrow(d))

  m_bare <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1, data = d, weights = wt, bootstrap = FALSE, verbose = FALSE,
    control = tight_ctrl_w()
  )))
  m_string <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1, data = d, weights = "wt", bootstrap = FALSE, verbose = FALSE,
    control = tight_ctrl_w()
  )))
  m_vector <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1, data = d, weights = rep(1, nrow(d)), bootstrap = FALSE, verbose = FALSE,
    control = tight_ctrl_w()
  )))
  m0 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1, data = d, bootstrap = FALSE, verbose = FALSE, control = tight_ctrl_w()
  )))

  expect_equal(m_bare$coef$Beta.NB, m0$coef$Beta.NB)
  expect_equal(m_string$coef$Beta.NB, m0$coef$Beta.NB)
  expect_equal(m_vector$coef$Beta.NB, m0$coef$Beta.NB)
})

test_that("unweighted fits are unchanged (weights = 1 by default)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_true(all(m$weights == 1))
  expect_equal(nobs(m), nrow(m$data$x.nb))
  expect_equal(glance(m)$sum_weights, glance(m)$nobs)

  mi <- fit_evinb_fast(bootstrap = FALSE)
  expect_true(all(mi$weights == 1))
  expect_equal(nobs(mi), nrow(mi$data$x.nb))
})

test_that("bootstrap replicates index weights by boot_id, not by position (round9 D.2)", {
  # Same pattern as the D.1 permuted-offset test: a fit with the weight
  # column permuted to different rows, refit with the same boot_seed (so the
  # *same* boot_id draws are used), must give different bootstrap
  # coefficients if the weight is genuinely carried through and indexed
  # alongside the resampled rows.
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  set.seed(8)
  d$w <- runif(nrow(d), 0.5, 3)

  m1 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d, weights = w,
    n_bootstraps = 2, boot_seed = 11, multicore = FALSE, verbose = FALSE,
    control = tight_ctrl_w()
  )))

  d2 <- d
  d2$w <- sample(d$w)
  m2 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d2, weights = w,
    n_bootstraps = 2, boot_seed = 11, multicore = FALSE, verbose = FALSE,
    control = tight_ctrl_w()
  )))

  expect_identical(m1$bootstraps[[1]]$boot_id, m2$bootstraps[[1]]$boot_id)
  c1 <- unname(m1$bootstraps[[1]]$coef$Beta.NB)
  c2 <- unname(m2$bootstraps[[1]]$coef$Beta.NB)
  expect_false(isTRUE(all.equal(c1, c2)))
})

test_that("update() reuses a weights column, not a raw weights vector", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$w <- rep(c(1, 2), length.out = nrow(d))

  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1, data = d, weights = w, bootstrap = FALSE, verbose = FALSE,
    control = tight_ctrl_w()
  )))
  expect_identical(m$weights_col, "w")

  m2 <- suppressMessages(suppressWarnings(update(m, formula_nb. = ~ . + x2)))
  expect_true(all(m2$weights == d$w))
})
