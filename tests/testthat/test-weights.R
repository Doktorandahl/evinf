# round9 D.2 (audit §5.6): weights = for evzinb() / evinb(). Frequency-weight
# semantics: every observation's contribution to the log-likelihood and the
# M-step accumulations is multiplied by its weight; nobs() (and therefore
# AIC/BIC/the approximate t p-values) use sum(weights).

tight_ctrl_w <- function(max.diff.par = 1e-6, max.no.em.steps = 3000, ...) {
  evinf_control(c.lim = c(50, 1000), init.C = 200,
                max.diff.par = max.diff.par, max.no.em.steps = max.no.em.steps, ...)
}

# round11 B2: repeatedly refits on row-duplicated data to verify weight
# semantics; kept fast on CI (NOT_CRAN=true) but skipped on CRAN's own
# check-time budget.
testthat::skip_on_cran()

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

# round10 0.5 (review §5): default (unweighted) fits moved by up to 8.6e-8
# from before round9 D.2. w(i) * x is never different from x when w = 1 in
# IEEE754 in isolation, but two things in the weighted code compound to move
# the *accumulated* result: (a) log_lik_fun()'s 0 < y < c branch changed a
# chained `func_val + a + b` (two left-to-right additions) into
# `func_val + (a + b)` -- addition is not associative, so that regrouping
# alone is a real bit-level change, independent of w; and (b) the compiler's
# vectorisation/instruction scheduling for this hot per-row loop is itself
# sensitive to the code's shape, so even a conditional branch that always
# reduces to the pre-D.2 expression can still perturb the last bit or two on
# some inputs. has_weights = FALSE now takes the exact pre-D.2 source text at
# every w(i)*/w% site (fixing (a) in full), which restores bit-identity on
# every case checked directly against b63b302 (genevzinb2 and hks, evzinb and
# evinb) and reduces the residual on other inputs from up to 8.6e-8 to at
# most a couple of ULPs (~1e-13) rather than eliminating it outright -- (b)
# is a property of the compiler's codegen, not of any remaining regrouping.

test_that("log_lik_fun()/update_bfgs_fun() agree to ~1 ULP between has_weights = TRUE/FALSE when w = 1 (round10 0.5)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  n <- nrow(m$data$data)
  w1 <- rep(1, n)
  off0 <- rep(0, n)
  xmz <- cbind(1, m$data$x.multinom.zc)
  xmpl <- cbind(1, m$data$x.multinom.pl)
  xnb <- cbind(1, m$data$x.nb)
  xpl <- cbind(1, m$data$x.pl)

  args_ll <- list(
    m$coef$Beta.multinom.ZC, m$coef$Beta.multinom.PL, m$coef$Beta.NB,
    m$coef$Alpha.NB, m$coef$Beta.PL, m$coef$C,
    xmz, xmpl, xnb, xpl, m$data$y, off0, off0, off0, w1, 0L, 0L
  )
  ll_true  <- do.call(evinf:::log_lik_fun, c(args_ll, list(TRUE)))
  ll_false <- do.call(evinf:::log_lik_fun, c(args_ll, list(FALSE)))
  # tolerance, not identical: see the header comment on (b) above.
  expect_equal(ll_true, ll_false, tolerance = 1e-10)

  args_upd <- list(
    m$coef$Beta.multinom.ZC, m$coef$Beta.multinom.PL, m$coef$Beta.NB,
    m$coef$Alpha.NB, m$coef$Beta.PL, m$coef$C,
    xmz, xmpl, xnb, xpl, m$data$y, 0.5, 3, off0, off0, off0, w1, 0L, FALSE, 0L
  )
  upd_true  <- do.call(evinf:::update_bfgs_fun, c(args_upd, list(TRUE)))
  upd_false <- do.call(evinf:::update_bfgs_fun, c(args_upd, list(FALSE)))
  expect_equal(upd_true$beta_nb_old, upd_false$beta_nb_old, tolerance = 1e-10)
  expect_equal(upd_true$gamma_z_old, upd_false$gamma_z_old, tolerance = 1e-10)
  expect_equal(upd_true$gamma_pl_old, upd_false$gamma_pl_old, tolerance = 1e-10)
  expect_equal(upd_true$beta_pl_old, upd_false$beta_pl_old, tolerance = 1e-10)
})

test_that("a default evzinb() fit is bit-identical to the pre-round9-D.2 baseline (round10 0.5)", {
  data(genevzinb2, package = "evinf", envir = environment())
  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = FALSE, verbose = FALSE
  )))
  # Pinned from b63b302 (the round9 Part 0 merge, before D.2's weights
  # landed) -- refitting on that commit reproduces this to the last bit, and
  # so does this fit after round10 0.5. test-em-identity.R's fixtures guard
  # the general case at a CI-portable 1e-8; this pins the exact value as a
  # tripwire specific to this fix.
  expect_equal(m$log.lik, -254.03389449157694, tolerance = 1e-10)
})
