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

test_that("insight::get_data() and avg_slopes() work without newdata (issue 3.3)", {
  skip_on_cran()
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("insight")
  m <- fit_evzinb_fast(n_bootstraps = 8)

  expect_equal(nrow(insight::get_data(m)), nrow(m$data$data))

  w <- character(0)
  s <- withCallingHandlers(
    marginaleffects::avg_slopes(m, variables = "x1"),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(w, 0L)
  expect_true(is.finite(s$estimate[1]))
})

test_that("marginal_effects(type = 'quantile') is non-zero (N1)", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  me <- suppressWarnings(marginal_effects(m, variables = "x1", type = "quantile",
                                          quantile = 0.9, method = "derivative"))
  expect_equal(nrow(me), 1L)
  expect_true(is.finite(me$estimate))
  expect_true(abs(me$estimate) > 1e-6)
})

test_that("quantiles_from_evzinb(round = FALSE) returns non-integer values", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  nd <- m$data$data[1:20, ]
  q_round <- evinf:::quantiles_from_evzinb(m, 0.9, newdata = nd, round = TRUE)
  q_cont  <- evinf:::quantiles_from_evzinb(m, 0.9, newdata = nd, round = FALSE)
  expect_equal(q_round, round(q_cont))
  expect_gt(sum(abs(q_cont - round(q_cont))), 0)  # at least some non-integers
})

test_that("marginal_effects(method =) is implemented for numeric covariates (N2)", {
  m <- fit_evzinb_fast(n_bootstraps = 6)

  d <- suppressWarnings(marginal_effects(m, variables = "x1", type = "harmonic",
                                         method = "derivative"))
  diff1 <- suppressWarnings(marginal_effects(m, variables = "x1", type = "harmonic",
                                             method = "difference", delta = 1))
  expect_true(is.finite(d$estimate) && is.finite(diff1$estimate))
  expect_identical(sign(d$estimate), sign(diff1$estimate))

  # type = "quantile" defaults to method = "difference"
  q_def <- suppressWarnings(marginal_effects(m, variables = "x1",
                                             type = "quantile", quantile = 0.9))
  q_der <- suppressWarnings(marginal_effects(m, variables = "x1",
                                             type = "quantile", quantile = 0.9,
                                             method = "derivative"))
  expect_true(is.finite(q_def$estimate) && abs(q_def$estimate) > 1e-6)
  expect_identical(sign(q_def$estimate), sign(q_der$estimate))

  # delta scales the difference roughly linearly for a smooth type
  diff2 <- suppressWarnings(marginal_effects(m, variables = "x1", type = "harmonic",
                                             method = "difference", delta = 2))
  expect_gt(abs(diff2$estimate), abs(diff1$estimate))
})

test_that("marginal_effects(type = 'quantile', n_max =) subsamples with a message (N8)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  expect_message(
    me <- marginal_effects(m, variables = "x1", type = "quantile",
                           quantile = 0.9, n_max = 10),
    "10-row subsample"
  )
  expect_true(is.finite(me$estimate))

  # reproducible: two runs give the same estimate
  me2 <- suppressMessages(marginal_effects(m, variables = "x1", type = "quantile",
                                           quantile = 0.9, n_max = 10))
  expect_equal(me$estimate, me2$estimate)

  # n_max = Inf uses all rows (no message)
  expect_no_message(
    marginal_effects(m, variables = "x1", type = "quantile", quantile = 0.9,
                     n_max = Inf)
  )
})
