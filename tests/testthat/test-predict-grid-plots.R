# audit 4.6 - predict_grid() and the plot methods.

test_that("predict_grid() returns a long tibble with the documented columns", {
  m <- fit_evzinb_fast(n_bootstraps = 5)

  g <- predict_grid(m, "x1", type = "harmonic", n = 12)
  expect_s3_class(g, "tbl_df")
  expect_identical(names(g),
                   c("variable", "value", "type", "estimate", "conf.low", "conf.high"))
  expect_equal(nrow(g), 12)
  expect_true(all(g$variable == "x1"))
  expect_true(all(is.finite(g$estimate)))
  expect_true(all(is.na(g$conf.low)))
})

test_that("predict_grid() holds the other covariates fixed at one value", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  # The held-fixed covariates must be constant across the grid: predictions for
  # a flat grid in x1 should trace a single smooth curve, and re-running with a
  # different `fixed` rule should move it.
  g_mean <- predict_grid(m, "x1", type = "harmonic", n = 8, fixed = "mean")
  g_med  <- predict_grid(m, "x1", type = "harmonic", n = 8, fixed = "median")
  expect_equal(nrow(g_mean), 8)
  expect_false(isTRUE(all.equal(g_mean$estimate, g_med$estimate)))

  # `at` pins a covariate explicitly.
  g_at <- predict_grid(m, "x1", type = "harmonic", n = 8, at = list(x2 = 0))
  expect_equal(nrow(g_at), 8)
})

test_that("predict_grid(type = 'states') pivots to one row per state", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  g <- predict_grid(m, "x1", type = "states", n = 6)
  expect_setequal(unique(g$type), c("pr_zero", "pr_count", "pr_evi"))
  expect_equal(nrow(g), 18)
  expect_message(predict_grid(m, "x1", type = "states", n = 4, confint = TRUE),
                 "not available")
})

test_that("predict_grid() adds bootstrap confidence intervals on request", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  g <- predict_grid(m, "x1", type = "harmonic", n = 6, confint = TRUE)
  expect_false(any(is.na(g$conf.low)))
  expect_false(any(is.na(g$conf.high)))
  expect_true(all(g$conf.low <= g$conf.high))
})

test_that("predict_grid() rejects an unknown variable", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_error(predict_grid(m, "not_a_var", type = "harmonic"), "must be one of")
})

test_that("each plot type returns a ggplot", {
  skip_if_not_installed("ggplot2")
  m <- fit_evzinb_fast(n_bootstraps = 5)

  expect_s3_class(plot(m, type = "states", variable = "x1"), "ggplot")
  expect_s3_class(plot(m, type = "prediction", variable = "x1"), "ggplot")
  expect_s3_class(plot(m, type = "coefficients"), "ggplot")
  expect_s3_class(plot(m, type = "ppc"), "ggplot")
  expect_s3_class(plot(m, type = "ppc_quantiles"), "ggplot")
})

test_that("ppc_quantiles works for evinb too", {
  skip_if_not_installed("ggplot2")
  m <- fit_evinb_fast(n_bootstraps = 5)
  expect_s3_class(plot(m, type = "ppc_quantiles"), "ggplot")
})

test_that("plot() needs a variable for states / prediction and bootstraps for coefficients", {
  skip_if_not_installed("ggplot2")
  m0 <- fit_evzinb_fast(bootstrap = FALSE)
  expect_error(plot(m0, type = "states"), "variable")
  expect_error(plot(m0, type = "prediction"), "variable")
  expect_error(plot(m0, type = "coefficients"), "bootstrap = TRUE")
})
