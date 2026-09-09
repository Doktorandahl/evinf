# audit 4.1 - standard S3 generics

test_that("coef() flattens with <component>_<term> names plus alpha_nb, c_ev", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  cf <- coef(m)
  expect_true(all(c("count_(Intercept)", "count_x1", "zero_x1", "evi_x1",
                    "pareto_x1", "alpha_nb", "c_ev") %in% names(cf)))
  expect_equal(unname(coef(m, "count")), unname(m$coef$Beta.NB))

  mi <- fit_evinb_fast(n_bootstraps = 6)
  expect_false(any(grepl("^zero_", names(coef(mi)))))
})

test_that("vcov() is the bootstrap covariance and errors without bootstraps", {
  m <- fit_evzinb_fast(n_bootstraps = 8)
  V <- vcov(m)
  expect_equal(rownames(V), names(coef(m)))
  expect_true(isSymmetric(unname(V)))

  m0 <- fit_evzinb_fast(bootstrap = FALSE)
  expect_error(vcov(m0), "bootstrap")
})

test_that("confint() returns percentile and approx intervals", {
  m <- fit_evzinb_fast(n_bootstraps = 10)
  ci <- confint(m)
  expect_equal(nrow(ci), length(coef(m)))
  expect_equal(ncol(ci), 2L)
  ci2 <- confint(m, parm = c("count_x1", "c_ev"), type = "approx")
  expect_equal(rownames(ci2), c("count_x1", "c_ev"))
})

test_that("logLik()/AIC()/BIC()/nobs() agree with the stored values", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  mi <- fit_evinb_fast(n_bootstraps = 5)
  expect_equal(as.numeric(logLik(m)), m$log.lik)
  expect_equal(attr(logLik(m), "df"), length(m$par.all))
  expect_equal(AIC(m), m$AIC)
  expect_equal(BIC(m), m$BIC)
  expect_equal(AIC(mi), mi$AIC)
  expect_equal(nobs(m), nrow(m$data$x.nb))
})

test_that("formula()/terms()/model.frame() dispatch per component", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_s3_class(formula(m, "zero"), "formula")
  expect_s3_class(terms(m, "pareto"), "terms")
  expect_identical(model.frame(m), m$data$data)
})

test_that("fitted() and residuals() work", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_equal(fitted(m), predict(m, type = "harmonic"))
  expect_equal(residuals(m), m$data$y - fitted(m))
  rq <- residuals(m, type = "quantile", seed = 1)
  expect_length(rq, nobs(m))
  expect_true(all(is.finite(rq)))
  # reproducible
  expect_equal(rq, residuals(m, type = "quantile", seed = 1))
})

test_that("simulate() returns a tidy data frame", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  s <- simulate(m, nsim = 3, seed = 1)
  expect_named(s, c("sim_1", "sim_2", "sim_3"))
  expect_equal(nrow(s), nobs(m))
  expect_equal(s, simulate(m, nsim = 3, seed = 1))
})

test_that("update() edits per-component formulas and other arguments", {
  m <- fit_evzinb_fast(n_bootstraps = 3)
  m2 <- suppressWarnings(suppressMessages(
    update(m, formula_pareto = . ~ . - x3, evaluate = TRUE)
  ))
  expect_equal(length(attr(stats::terms(formula(m2, "pareto")), "term.labels")), 2L)
  m3 <- suppressWarnings(suppressMessages(update(m, bootstrap = FALSE)))
  expect_null(m3$bootstraps)
})

test_that("modelsummary works with a grouped shape (skip if absent)", {
  skip_if_not_installed("modelsummary")
  m <- fit_evzinb_fast(n_bootstraps = 6)
  # As with nnet::multinom, modelsummary needs the y.level grouping in `shape`.
  expect_no_error(
    modelsummary::modelsummary(m, shape = term + y.level ~ model,
                               output = "data.frame")
  )
})
