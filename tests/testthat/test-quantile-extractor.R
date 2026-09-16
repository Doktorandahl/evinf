# audit0.10 §1.13 (D.5) - quantiles_from_evinb() gets the same near-zero
# pareto-alpha clamp+warning as quantiles_from_evzinb(), via a shared helper.

test_that("evinf_clamp_pareto_alpha() clamps values below 1e-02 and warns", {
  expect_warning(clamped <- evinf:::evinf_clamp_pareto_alpha(c(0.001, 0.5)),
                 "below 1e-02")
  expect_equal(clamped, c(0.01, 0.5))
  expect_no_warning(evinf:::evinf_clamp_pareto_alpha(c(0.5, 0.9)))
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
                 "below 1e-02")
  expect_true(all(is.finite(q)))
})
