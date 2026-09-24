# round11 B2: broad summary()/tidy() coverage across bootstrap option
# combinations; kept fast on CI (NOT_CRAN=true) but skipped on CRAN's own
# check-time budget.
testthat::skip_on_cran()

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

test_that("summary(standard_error = FALSE) works for every approx_t_value/p_value combination (audit0.10 §1.2)", {
  fits <- list(evzinb = fit_evzinb_fast(n_bootstraps = 6), evinb = fit_evinb_fast(n_bootstraps = 6))
  for (fit in fits) {
    for (se in c(TRUE, FALSE)) {
      for (at in c(TRUE, FALSE)) {
        for (pv in c("none", "approx", "bootstrapped", "both")) {
          info <- sprintf("class=%s se=%s at=%s pv=%s", class(fit)[1], se, at, pv)
          s <- suppressWarnings(summary(fit, standard_error = se, approx_t_value = at, p_value = pv))
          final_se <- se || pv %in% c("approx", "both")
          expect_equal(unname("se" %in% names(s$coefficients$count)), final_se, info = info)
          expect_equal(unname("approx_t" %in% names(s$coefficients$count)), final_se && at, info = info)
          expect_no_error(print(s))
        }
      }
    }
  }
})

test_that("tidy(standard_error = FALSE) works for every approx_t_value/p_value/confint combination (audit0.10 §1.2)", {
  fits <- list(evzinb = fit_evzinb_fast(n_bootstraps = 6), evinb = fit_evinb_fast(n_bootstraps = 6))
  for (fit in fits) {
    for (se in c(TRUE, FALSE)) {
      for (at in c(TRUE, FALSE)) {
        for (pv in c("none", "approx", "bootstrapped")) {
          for (ci in c("none", "approx", "bootstrapped")) {
            info <- sprintf("class=%s se=%s at=%s pv=%s ci=%s", class(fit)[1], se, at, pv, ci)
            td <- suppressWarnings(tidy(fit, component = "count", standard_error = se,
                                        approx_t_value = at, p_value = pv, confint = ci))
            final_se <- se || pv == "approx" || ci == "approx"
            expect_equal(unname("std.error" %in% names(td)), final_se, info = info)
            expect_equal(unname("statistic" %in% names(td)), final_se && at, info = info)
          }
        }
      }
    }
  }
})

test_that("bootstrap_p_value_calculator() floors at 1/B instead of returning 0 (audit0.10 §1.11)", {
  x_all_positive <- rep(1, 5)  # every draw on the same side -> raw p would be 0
  expect_equal(evinf:::bootstrap_p_value_calculator(x_all_positive, estimate = 1), 1 / 5)
  expect_equal(evinf:::bootstrap_p_value_calculator(x_all_positive, estimate = 1, symmetric = FALSE), 1 / 5)

  x_mixed <- c(-1, 1, 1, 1, 1)  # 1/5 crossed zero -> raw p = 2/5, above the floor
  expect_equal(evinf:::bootstrap_p_value_calculator(x_mixed, estimate = 1), 2 / 5)
})

test_that("print(summary()) shows '< 1/B' instead of '<2e-16' for a floored bootstrap p-value (audit0.10 §1.11)", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  # Force every usable replicate's x1 (count component) to the same sign, so
  # its bootstrapped p-value is genuinely floored at 1/n_bootstraps_used.
  for (i in seq_along(m$bootstraps)) {
    if (!inherits(m$bootstraps[[i]], "try-error")) {
      m$bootstraps[[i]]$coef$Beta.NB["x1"] <- abs(m$bootstraps[[i]]$coef$Beta.NB["x1"]) + 1
    }
  }
  s <- suppressWarnings(summary(m, p_value = "bootstrapped"))
  out <- capture.output(print(s))
  # the row we forced to one sign (count component's x1) shows a
  # "<"-prefixed (floored) p-value, not the raw exact-zero value.
  x1_line <- out[grepl("^x1\\s", out)][1]
  expect_match(x1_line, "<")
  expect_false(any(grepl("2e-16", out)))
})

test_that("summary() works when exactly one bootstrap replicate is usable (round11)", {
  # Found via R CMD check --run-donttest on ?evzinb's hks example:
  # purrr::reduce(rbind) over a length-1 list returns the bare colMeans()
  # vector unchanged (never calling rbind()), so as_tibble(.name_repair =
  # ~prop_names) read it as 3 rows of 1 column instead of 1 row of 3 columns
  # and errored ("Repaired names have length 3 instead of length 1").
  m <- fit_evzinb_fast(n_bootstraps = 2)
  m$bootstraps[[1]] <- try(stop("boom"), silent = TRUE)
  expect_no_error(s <- suppressWarnings(summary(m)))
  expect_true(all(c("zero", "count", "evi") %in% s$component_proportions$state))
})
