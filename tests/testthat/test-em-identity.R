# Numerical-identity tests for the R/a019.R -> R/em_*.R refactor (round 3, Part B).
#
# The baselines in tests/testthat/fixtures/em_baseline_*.rds were captured from
# the pre-refactor code by tests/testthat/fixtures/make_baselines.R. The refactor
# must reproduce them exactly.

ctrl_id <- function() evinf::evinf_control(c.lim = c(50, 1000), init.C = 200)

fit_z_id <- function(formula, data, ...) {
  suppressMessages(suppressWarnings(
    evinf::evzinb(formula, data = data, control = ctrl_id(), verbose = FALSE, ...)
  ))
}
fit_i_id <- function(formula, data, ...) {
  suppressMessages(suppressWarnings(
    evinf::evinb(formula, data = data, control = ctrl_id(), verbose = FALSE, ...)
  ))
}

capture_fit_id <- function(m) {
  list(
    beta_nb = m$coef$Beta.NB, beta_zc = m$coef$Beta.multinom.ZC,
    beta_pl_mult = m$coef$Beta.multinom.PL, beta_pl = m$coef$Beta.PL,
    alpha_nb = m$coef$Alpha.NB, c_ev = m$coef$C,
    log_lik = m$log.lik, aic = m$AIC, bic = m$BIC,
    props = m$props, resp = m$resp, fitted = m$fitted,
    converge = m$converge, n_above_c = m$n_above_c,
    loglik_recomputed = isTRUE(m$loglik_recomputed), c_profile = m$c_profile
  )
}

expect_fit_equal <- function(actual, baseline_file, tol) {
  base <- readRDS(test_path("fixtures", baseline_file))
  expect_equal(actual$beta_nb,      base$beta_nb,      tolerance = tol)
  expect_equal(actual$beta_zc,      base$beta_zc,      tolerance = tol)
  expect_equal(actual$beta_pl_mult, base$beta_pl_mult, tolerance = tol)
  expect_equal(actual$beta_pl,      base$beta_pl,      tolerance = tol)
  expect_equal(actual$alpha_nb,     base$alpha_nb,     tolerance = tol)
  expect_equal(actual$c_ev,         base$c_ev,         tolerance = tol)
  expect_equal(actual$log_lik,      base$log_lik,      tolerance = tol)
  expect_equal(actual$aic,          base$aic,          tolerance = tol)
  expect_equal(actual$bic,          base$bic,          tolerance = tol)
  expect_equal(unname(actual$props), unname(base$props), tolerance = tol)
  expect_equal(unname(actual$resp),  unname(base$resp),  tolerance = tol)
  expect_equal(actual$fitted,       base$fitted,       tolerance = tol)
  expect_identical(actual$converge, base$converge)
  expect_identical(actual$n_above_c, base$n_above_c)
  expect_identical(actual$loglik_recomputed, base$loglik_recomputed)
  expect_equal(as.data.frame(actual$c_profile), as.data.frame(base$c_profile),
               tolerance = tol)
}

data(genevzinb2, package = "evinf", envir = environment())
gf_id <- local({
  set.seed(42)
  d <- genevzinb2
  d$g <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  d
})

test_that("evzinb() estimation is identical after the em_*.R refactor", {
  m <- fit_z_id(y ~ x1 + x2 + x3, genevzinb2, bootstrap = FALSE)
  expect_fit_equal(capture_fit_id(m), "em_baseline_evzinb.rds", 1e-8)
})

test_that("evinb() estimation is identical after the em_*.R refactor", {
  m <- fit_i_id(y ~ x1 + x2 + x3, genevzinb2, bootstrap = FALSE)
  expect_fit_equal(capture_fit_id(m), "em_baseline_evinb.rds", 1e-8)
})

test_that("estimation is identical with a factor + in-formula transform", {
  f3 <- y ~ x1 + log(abs(x2) + 1) + g
  expect_fit_equal(capture_fit_id(fit_z_id(f3, gf_id, bootstrap = FALSE)),
                   "em_baseline_evzinb_factor.rds", 1e-8)
  expect_fit_equal(capture_fit_id(fit_i_id(f3, gf_id, bootstrap = FALSE)),
                   "em_baseline_evinb_factor.rds", 1e-8)
})

test_that("bootstrap coefficient matrices are identical (n = 3, seed 123)", {
  base <- readRDS(test_path("fixtures", "em_baseline_boot.rds"))
  mzb <- fit_z_id(y ~ x1 + x2 + x3, genevzinb2, bootstrap = TRUE,
                  n_bootstraps = 3, boot_seed = 123, multicore = FALSE)
  mib <- fit_i_id(y ~ x1 + x2 + x3, genevzinb2, bootstrap = TRUE,
                  n_bootstraps = 3, boot_seed = 123, multicore = FALSE)
  # compare every replicate (the fixture has all 3); a degenerate one would
  # otherwise be dropped by coefficient_extractor()'s default.
  expect_equal(as.data.frame(suppressWarnings(
                 coefficient_extractor(mzb, "all", exclude_degenerate = FALSE))),
               as.data.frame(base$evzinb), tolerance = 1e-8)
  expect_equal(as.data.frame(suppressWarnings(
                 coefficient_extractor(mib, "all", exclude_degenerate = FALSE))),
               as.data.frame(base$evinb), tolerance = 1e-8)
})

test_that("hks estimation is identical after the em_*.R refactor", {
  skip_on_cran()
  data(hks, package = "evinf", envir = environment())
  f_hks <- osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
    lntpop + brv_AllLag_log + osvAllLagDum + incomp
  f_hks_pareto <- ~ troopLag_log + epdur_log + brv_AllLag_log + lntpop

  expect_fit_equal(
    capture_fit_id(fit_z_id(f_hks, hks, formula_pareto = f_hks_pareto,
                            bootstrap = FALSE)),
    "em_baseline_evzinb_hks.rds", 1e-6
  )
  expect_fit_equal(
    capture_fit_id(fit_i_id(f_hks, hks, formula_pareto = f_hks_pareto,
                            bootstrap = FALSE)),
    "em_baseline_evinb_hks.rds", 1e-6
  )
})
