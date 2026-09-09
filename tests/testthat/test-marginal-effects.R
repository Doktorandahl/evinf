# audit 4.8 - marginal_effects() and marginaleffects compatibility.

test_that("marginal_effects() returns the documented tibble with finite estimates", {
  m <- fit_evzinb_fast(n_bootstraps = 8)
  me <- suppressWarnings(marginal_effects(m, variables = "x1"))

  expect_s3_class(me, "tbl_df")
  expect_identical(names(me),
                   c("variable", "contrast", "type", "estimate",
                     "std.error", "conf.low", "conf.high"))
  expect_equal(nrow(me), 1)
  expect_true(is.finite(me$estimate))
  expect_identical(me$contrast, "dydx")
})

test_that("marginal_effects(type = 'states') gives one row per state", {
  m <- fit_evzinb_fast(n_bootstraps = 8)
  me <- suppressWarnings(marginal_effects(m, variables = "x1", type = "states"))
  expect_equal(nrow(me), 3)
  expect_setequal(me$type, c("pr_zero", "pr_count", "pr_evi"))
  expect_true(all(is.finite(me$estimate)))
})

test_that("marginal_effects() contrasts a factor against its reference level", {
  d <- genevzinb2_factor()
  m <- suppressMessages(evinf::evzinb(
    y ~ x1 + g, data = d, control = .fast_control(),
    bootstrap = TRUE, n_bootstraps = 8, multicore = FALSE,
    boot_seed = 123, verbose = FALSE
  ))
  me <- suppressWarnings(marginal_effects(m, variables = "g"))
  expect_setequal(me$contrast, c("b - a", "c - a"))
  expect_true(all(is.finite(me$estimate)))
})

test_that("marginal_effects() needs bootstraps", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_error(marginal_effects(m, variables = "x1"), "bootstrap = TRUE")
})

test_that("avg_slopes() agrees with marginal_effects() for a numeric covariate", {
  skip_on_cran()
  skip_if_not_installed("marginaleffects")
  m <- fit_evzinb_fast(n_bootstraps = 8)

  me <- suppressWarnings(marginal_effects(m, variables = "x1"))
  s <- suppressWarnings(
    marginaleffects::avg_slopes(m, variables = "x1", newdata = m$data$data)
  )
  expect_equal(unname(s$estimate[1]), me$estimate[1], tolerance = 0.1)
})
