# audit 4.7 - compare_fit(), oob_evaluation() and the evzinbcomp tidiers.

make_comp <- function() {
  m <- fit_evzinb_fast(n_bootstraps = 8)
  suppressWarnings(suppressMessages(compare_models(m)))
}

# round11 B2: dominated by bootstrap-heavy compare_models()/compare_fit()
# fits; kept fast on CI (NOT_CRAN=true) but skipped on CRAN's own
# check-time budget.
testthat::skip_on_cran()

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

  expect_output(print(cf), "evinf - compared")
})

test_that("compare_fit() sign: evinf clearly winning on AIC gives a negative median and prop_evinf_better near 1 (audit0.10 §1.6)", {
  comp <- make_comp()
  # evzinb() adds two multinomial-logit blocks (zero-inflation, extreme-value)
  # on top of the NB mean, so on data actually generated from an EVZINB
  # process (genevzinb2) it comfortably beats the plain NB on AIC.
  cf <- suppressWarnings(compare_fit(comp, metrics = "aic"))
  nb_row <- cf[cf$model == "nb", ]
  expect_lt(nb_row$median_difference, 0)
  expect_gt(nb_row$prop_evinf_better, 0.9)
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

test_that("compare_fit() returns NA AIC/BIC for razorised/winsorised slots (round8 0.5, review §6)", {
  m <- fit_evzinb_fast(n_bootstraps = 8)
  comp <- suppressWarnings(suppressMessages(
    compare_models(m, winsorize = TRUE, razorize = TRUE, cutoff_value = 20)
  ))
  cf <- suppressWarnings(compare_fit(comp, metrics = c("aic", "bic", "rmse", "rmsle")))

  modified <- cf[grepl("_razor$|_winsor$", cf$model), ]
  expect_true(all(is.na(modified$median_difference[modified$metric %in% c("aic", "bic")])))
  expect_true(all(is.na(modified$prop_evinf_better[modified$metric %in% c("aic", "bic")])))
  # RMSE/RMSLE stay comparable (OOB error is always against the raw outcome).
  expect_true(all(is.finite(modified$median_difference[modified$metric %in% c("rmse", "rmsle")])))

  unmodified <- cf[!grepl("_razor$|_winsor$", cf$model), ]
  expect_true(all(is.finite(unmodified$median_difference)))

  expect_output(print(cf), "not comparable to the")

  # round9 0.6 (review §6): n_pairs stays numeric ("as computed") on the
  # object itself, but print() blanks it for the not-comparable rows instead
  # of showing a count next to an NA metric. A comparable row (rmse) still
  # ends in its numeric n_pairs; a not-comparable one (aic) ends right after
  # the NA/NA columns, with nothing trailing.
  expect_true(all(modified$n_pairs[modified$metric %in% c("aic", "bic")] > 0))
  printed <- capture.output(print(cf))
  razor_aic_line <- grep("^\\s*nb_razor\\s+aic\\s", printed, value = TRUE)
  razor_rmse_line <- grep("^\\s*nb_razor\\s+rmse\\s", printed, value = TRUE)
  expect_length(razor_aic_line, 1L)
  expect_length(razor_rmse_line, 1L)
  expect_match(razor_aic_line, "NA\\s+NA\\s*$")
  expect_match(razor_rmse_line, "[0-9]\\s*$")
})

test_that("oob_evaluation() on an evzinbcomp masks every column at the same positions (round8 0.4, review §5)", {
  comp <- make_comp()
  # Force one evinf bootstrap replicate to look degenerate (without needing a
  # genuinely degenerate fit) so its oob error comes back NA -- the nb/zinb
  # columns must then be NA at that same position too.
  comp$model$bootstraps[[1]]$degenerate <- TRUE

  oe <- suppressWarnings(oob_evaluation(comp, metric = "rmse"))
  expect_true(all(is.na(oe[1, ])))

  n_nonmissing <- vapply(oe, function(col) sum(!is.na(col)), integer(1))
  expect_length(unique(n_nonmissing), 1L)

  excluded <- attr(oe, "excluded")
  expect_true(is.logical(excluded))
  expect_length(excluded, nrow(oe))
  expect_true(excluded[1])
  expect_identical(which(excluded), which(apply(is.na(oe), 1, any)))
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
