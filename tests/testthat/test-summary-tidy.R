test_that("summary(p_value = 'both') returns both p-value columns (audit 1.3)", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  s <- suppressWarnings(summary(m, p_value = "both"))
  expect_true(all(c("bootstrap_p", "approx_p") %in% names(s$coefficients$count)))
})

test_that("summary()/tidy() work without bootstraps (audit 1.4)", {
  m0 <- fit_evzinb_fast(bootstrap = FALSE)

  expect_message(s <- summary(m0), "bootstrap")
  expect_false("se" %in% names(s$coefficients$count))
  expect_true(is.na(s$n_failed_bootstraps))
  expect_error(summary(m0, coef = "bootstrapped_mean"), "bootstrap = TRUE")

  expect_message(td <- tidy(m0), "bootstrap")
  expect_setequal(names(td), c("y.level", "term", "estimate"))
  expect_error(tidy(m0, coef_type = "bootstrap_mean"), "bootstrap = TRUE")
})

test_that("tidy() default component is 'all' (audit 2.8)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  td <- suppressWarnings(tidy(m))
  expect_true("y.level" %in% names(td))
  expect_equal(levels(td$y.level), c("zero", "evi", "count", "pareto"))
})

test_that("tidy(confint=) produces conf.low/conf.high (audit 2.9)", {
  m <- fit_evzinb_fast(n_bootstraps = 8)

  tb <- suppressWarnings(tidy(m, component = "count", confint = "bootstrapped", conf_level = 0.9))
  expect_true(all(c("conf.low", "conf.high") %in% names(tb)))
  expect_true(all(tb$conf.low <= tb$conf.high))

  ta <- suppressWarnings(tidy(m, component = "count", confint = "approx", conf_level = 0.9))
  nobs <- nrow(m$data$x.nb); npar <- length(m$par.all)
  crit <- qt(0.95, df = nobs - npar)
  expect_equal(ta$conf.low, ta$estimate - crit * ta$std.error)
})

test_that("summary()/tidy() emit no join messages (audit 2.14)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  expect_no_message(summary(m))
  expect_no_message(suppressWarnings(tidy(m, component = "count")))
})
