# Coverage for the nbboot / zinbboot tidy, predict and coefficient_extractor
# methods reached through compare_models(), and predict.evinb() bootstrap paths.

cached_comp <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      m <- fit_evzinb_fast(n_bootstraps = 8)
      cache <<- suppressWarnings(suppressMessages(compare_models(m)))
    }
    cache
  }
})

test_that("tidy.nbboot() / tidy.zinbboot() honour their options", {
  comp <- cached_comp()

  for (ct in c("original", "bootstrap_mean", "bootstrap_median")) {
    tn <- suppressWarnings(generics::tidy(comp$nb, coef_type = ct))
    expect_s3_class(tn, "tbl_df")
    expect_true("estimate" %in% names(tn))
  }
  tz <- suppressWarnings(generics::tidy(comp$zinb, component = "all",
                                       p_value = "approx",
                                       confint = "bootstrapped"))
  expect_true(all(c("conf.low", "conf.high") %in% names(tz)))

  tz2 <- suppressWarnings(generics::tidy(comp$zinb, component = "count",
                                        confint = "approx", p_value = "none"))
  expect_true("term" %in% names(tz2))
})

test_that("coefficient_extractor() works for every supported class", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  mi <- fit_evinb_fast(n_bootstraps = 6)
  comp <- cached_comp()

  expect_s3_class(coefficient_extractor(m, "all"), "tbl_df")
  expect_s3_class(coefficient_extractor(mi, "all"), "tbl_df")
  expect_s3_class(coefficient_extractor(comp$zinb, "all"), "tbl_df")
  expect_s3_class(coefficient_extractor(comp$nb), "tbl_df")
})

test_that("predict.zinbboot() / predict.nbboot() cover pred, confint, quantile", {
  comp <- cached_comp()
  n <- nrow(comp$model$data$data)

  for (ty in c("predicted", "counts", "zi", "count_state", "states")) {
    p <- suppressWarnings(predict(comp$zinb, type = ty))
    if (ty == "states") expect_equal(nrow(p), n) else
      expect_length(as.numeric(unlist(p)), n)
  }
  qz <- suppressWarnings(predict(comp$zinb, type = "quantile", quantile = 0.9))
  expect_length(as.numeric(unlist(qz)), n)

  cz <- suppressWarnings(predict(comp$zinb, type = "predicted", confint = TRUE))
  expect_true(all(c("ci_lb", "ci_ub") %in% names(cz)))

  # regression: predict.zinbboot(type = "counts", pred = "original") used to
  # error with "$ operator is invalid for atomic vectors".
  cc <- suppressWarnings(predict(comp$zinb, type = "counts"))
  expect_length(as.numeric(unlist(cc)), n)
  cca <- suppressWarnings(predict(comp$zinb, type = "counts", confint = TRUE))
  expect_true(all(c("count", "ci_lb", "ci_ub") %in% names(cca)))
  expect_equal(nrow(cca), n)

  cn <- suppressWarnings(predict(comp$nb, type = "predicted", confint = TRUE))
  expect_true(all(c("ci_lb", "ci_ub") %in% names(cn)))
  qn <- suppressWarnings(predict(comp$nb, type = "quantile", quantile = 0.5))
  expect_length(as.numeric(unlist(qn)), n)
})

test_that("predict.evinb() bootstrap summaries and CIs work", {
  m <- fit_evinb_fast(n_bootstraps = 8)
  n <- nrow(m$data$data)

  for (pr in c("bootstrap_median", "bootstrap_mean")) {
    p <- suppressWarnings(predict(m, type = "harmonic", pred = pr))
    expect_length(as.numeric(unlist(p)), n)
  }
  ci <- suppressWarnings(predict(m, type = "harmonic", confint = TRUE))
  expect_true(all(c("ci_lb", "ci_ub") %in% names(ci)))
  ciq <- suppressWarnings(predict(m, type = "quantile", quantile = 0.75,
                                  confint = TRUE))
  expect_true(all(c("ci_lb", "ci_ub") %in% names(ciq)))
})
