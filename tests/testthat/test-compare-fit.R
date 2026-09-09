# audit 4.7 - compare_fit(), oob_evaluation() and the evzinbcomp tidiers.

make_comp <- function() {
  m <- fit_evzinb_fast(n_bootstraps = 8)
  suppressWarnings(suppressMessages(compare_models(m)))
}

test_that("compare_fit() returns the documented shape and class", {
  comp <- make_comp()
  cf <- suppressWarnings(compare_fit(comp, metrics = c("aic", "bic")))

  expect_s3_class(cf, "evinf_compare_fit")
  expect_s3_class(cf, "tbl_df")
  expect_identical(names(cf),
                   c("model", "metric", "median_difference",
                     "prop_evinf_better", "n_pairs"))
  expect_setequal(unique(cf$metric), c("aic", "bic"))
  expect_setequal(unique(cf$model), c("nb", "zinb"))
  expect_true(all(cf$prop_evinf_better >= 0 & cf$prop_evinf_better <= 1))
  expect_true(all(cf$n_pairs > 0))

  expect_output(print(cf), "compared - evinf")
})

test_that("compare_fit() can add out-of-bag RMSE / RMSLE", {
  comp <- make_comp()
  cf <- suppressWarnings(compare_fit(comp, metrics = c("rmse", "rmsle")))
  expect_setequal(unique(cf$metric), c("rmse", "rmsle"))
  expect_true(all(is.finite(cf$median_difference)))
})

test_that("oob_evaluation() on an evzinbcomp gives one column per model", {
  comp <- make_comp()
  oe <- suppressWarnings(oob_evaluation(comp, metric = "rmse"))
  expect_s3_class(oe, "tbl_df")
  expect_true(all(c("evinf", "nb", "zinb") %in% names(oe)))
  expect_equal(nrow(oe), 8)
})

test_that("tidy() / glance() on an evzinbcomp stack the models with a model column", {
  comp <- make_comp()
  td <- suppressWarnings(generics::tidy(comp))
  gl <- suppressWarnings(generics::glance(comp))
  expect_true("model" %in% names(td))
  expect_true("model" %in% names(gl))
  expect_true(all(c("model", "nb", "zinb") %in% td$model))
})

test_that("plot.evzinbcomp() returns a ggplot", {
  skip_if_not_installed("ggplot2")
  comp <- make_comp()
  expect_s3_class(plot(comp, metrics = c("aic", "bic")), "ggplot")
})

test_that("compare_fit() rejects a non-evzinbcomp object", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_error(compare_fit(m), "evzinbcomp")
})
