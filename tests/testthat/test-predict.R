# Broad coverage of predict.evzinb() / predict.evinb() and the *boot predict
# methods reached through compare_models().

evzinb_types <- c("harmonic", "explog", "counts", "pareto_alpha", "zi", "evinf",
                  "count_state", "states", "all")
evinb_types  <- c("harmonic", "explog", "counts", "pareto_alpha", "evinf",
                  "count_state", "states", "all")

test_that("predict.evzinb() returns a value for every prediction type", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  n <- nrow(m$data$data)
  for (ty in evzinb_types) {
    p <- suppressWarnings(predict(m, type = ty))
    if (ty %in% c("states", "all")) {
      expect_equal(nrow(p), n, info = ty)
    } else {
      expect_length(as.numeric(unlist(p)), n)
    }
  }
  q <- suppressWarnings(predict(m, type = "quantile", quantile = 0.9))
  expect_length(as.numeric(q), n)
})

test_that("predict.evinb() returns a value for every prediction type", {
  m <- fit_evinb_fast(n_bootstraps = 5)
  n <- nrow(m$data$data)
  for (ty in evinb_types) {
    p <- suppressWarnings(predict(m, type = ty))
    if (ty %in% c("states", "all")) {
      expect_equal(nrow(p), n, info = ty)
    } else {
      expect_length(as.numeric(unlist(p)), n)
    }
  }
})

test_that("predict.evzinb() bootstrap summaries and confidence intervals work", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  n <- nrow(m$data$data)

  for (pr in c("bootstrap_median", "bootstrap_mean")) {
    p <- suppressWarnings(predict(m, type = "harmonic", pred = pr))
    expect_length(as.numeric(unlist(p)), n)
  }

  ci <- suppressWarnings(predict(m, type = "harmonic", confint = TRUE))
  expect_true(all(c("ci_lb", "ci_ub") %in% names(ci)))
  expect_equal(nrow(ci), n)

  ciq <- suppressWarnings(predict(m, type = "quantile", quantile = 0.8,
                                 confint = TRUE))
  expect_true(all(c("ci_lb", "ci_ub") %in% names(ciq)))

  rb <- suppressWarnings(predict(m, type = "harmonic", confint = TRUE,
                                 return_bootstraps = TRUE))
  expect_type(rb, "list")
})

test_that("predict.evzinb() accepts newdata and errors on a bad quantile call", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  nd <- m$data$data[1:10, ]
  expect_equal(nrow(suppressWarnings(predict(m, newdata = nd, type = "states"))), 10)
  expect_length(suppressWarnings(predict(m, newdata = nd, type = "harmonic")), 10)
  expect_error(suppressWarnings(predict(m, type = "quantile")), "quantile")
  expect_error(predict(m, type = "states", confint = TRUE), "vector outputs")
})

test_that("predict on the nb / zinb bootstrap objects works via compare_models()", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  comp <- suppressWarnings(suppressMessages(compare_models(m)))

  for (ty in c("original", "bootstrap_median", "bootstrap_mean")) {
    pnb <- suppressWarnings(predict(comp$nb, pred = ty))
    pzi <- suppressWarnings(predict(comp$zinb, pred = ty))
    expect_length(as.numeric(unlist(pnb)), nrow(m$data$data))
    expect_length(as.numeric(unlist(pzi)), nrow(m$data$data))
  }
  cinb <- suppressWarnings(predict(comp$nb, confint = TRUE))
  expect_true(all(c("ci_lb", "ci_ub") %in% names(cinb)))
})

test_that("compare_models() winsorize / razorize branches run", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  comp <- suppressWarnings(suppressMessages(
    compare_models(m, winsorize = TRUE, razorize = TRUE, cutoff_value = 30)
  ))
  expect_true(all(c("nb", "zinb", "nb_winsor", "zinb_winsor",
                    "nb_razor", "zinb_razor") %in% names(comp)))
  expect_s3_class(suppressWarnings(generics::glance(comp)), "tbl_df")
})
