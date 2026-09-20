# audit0.10 §1.13 (D.5) - quantiles_from_evinb() gets the same near-zero
# pareto-alpha clamp+warning as quantiles_from_evzinb(), via a shared helper
# (renamed evinf_clamp_alpha_pl() and generalised to harmonic_calc(),
# explog_calc() and em_fitted_values() too in round9 0.1, review §2).

test_that("evinf_clamp_alpha_pl() clamps values below the floor and warns with the count (round9 0.1)", {
  expect_warning(clamped <- evinf:::evinf_clamp_alpha_pl(c(0.001, 0.5)),
                 "1 fitted Pareto alpha value")
  expect_equal(clamped, c(0.01, 0.5))
  expect_no_warning(evinf:::evinf_clamp_alpha_pl(c(0.5, 0.9)))

  expect_warning(
    clamped2 <- evinf:::evinf_clamp_alpha_pl(c(0.001, 0.002, 0.5), floor = 0.01),
    "2 fitted Pareto alpha values"
  )
  expect_equal(clamped2, c(0.01, 0.01, 0.5))

  expect_warning(
    evinf:::evinf_clamp_alpha_pl(c(0.05, 0.5), floor = 0.1, context = "test context"),
    "test context"
  )
})

test_that("quantiles_from_evinb() clamps near-zero pareto alpha and warns, matching quantiles_from_evzinb() (audit0.10 §1.13, D.5)", {
  m <- fit_evinb_fast(bootstrap = FALSE)
  n <- nobs(m)
  testthat::local_mocked_bindings(
    fitted_alpha_from_evzinb = function(object, newdata = NULL, return_data = FALSE) {
      tibble::tibble(pareto_alpha = rep(0.001, n))
    },
    .package = "evinf"
  )
  expect_warning(q <- evinf:::quantiles_from_evinb(m, quantile = 0.5),
                 "below the floor")
  expect_true(all(is.finite(q)))
})
