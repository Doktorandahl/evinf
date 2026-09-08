test_that("compare_models() fits the bootstrapped ZINB with the full formula (audit 1.2)", {
  m <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + x2, formula_zi = ~x1,
    data = { data(genevzinb2, package = "evinf", envir = environment()); genevzinb2 },
    bootstrap = TRUE, n_bootstraps = 5, multicore = FALSE, boot_seed = 123, verbose = FALSE
  )))
  comp <- suppressWarnings(suppressMessages(compare_models(m)))

  ok <- comp$zinb$bootstraps[!vapply(comp$zinb$bootstraps, inherits, logical(1), "try-error")]
  expect_gt(length(ok), 0)
  for (b in ok) {
    expect_equal(length(b$coefficients$zero), 2L)  # (Intercept) + x1
    expect_equal(length(b$coefficients$count), 3L) # (Intercept) + x1 + x2
  }
})

test_that("compare_models() supports evinb and renames the first slot (audit 1.6)", {
  mi <- fit_evinb_fast(n_bootstraps = 5)

  expect_message(comp <- compare_models(mi), "zinb_comparison")
  expect_identical(comp$model, comp$evzinb)
  expect_false("zinb" %in% names(comp))
  expect_error(compare_models(mi, zinb_comparison = TRUE), "not available")
})

test_that("predict.nbboot()/predict.zinbboot() are registered (audit 2.5)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  comp <- suppressWarnings(suppressMessages(compare_models(m)))
  expect_true(is.numeric(predict(comp$nb, pred = "bootstrap_median")))
  expect_true(is.numeric(predict(comp$zinb, pred = "bootstrap_median")))
})
