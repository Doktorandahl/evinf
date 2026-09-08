test_that("lr_test() reports 2 * (logLik_full - logLik_restricted) (audit 1.1)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  res <- suppressWarnings(lr_test(m, "x1"))

  expect_s3_class(res, "tbl_df")
  expect_equal(res$statistic, 2 * (res$loglik_full - res$loglik_restricted))
  expect_equal(res$prob, pchisq(res$statistic, res$df, lower.tail = FALSE))
  # restricting x1 removes one column from each of the four components
  expect_equal(res$df, 4L)
})

test_that("lr_test(bootstrap = TRUE) doubles the bootstrap statistics too (audit 1.1)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  res <- suppressWarnings(lr_test(m, "x1", bootstrap = TRUE))

  expect_named(res, c("results", "boot_results"))
  bs <- res$boot_results[["x1"]]
  ok <- !is.na(bs$statistic)
  expect_equal(bs$statistic[ok], 2 * (bs$ll_full[ok] - bs$ll_reduced[ok]))
})

test_that("lr_test() works on evinb objects (audit 1.5)", {
  mi <- fit_evinb_fast(bootstrap = FALSE)
  res <- suppressWarnings(lr_test(mi, "x1"))
  expect_s3_class(res, "tbl_df")
  # x1 is in nb, evi and pareto (no zero-inflation component)
  expect_equal(res$df, 3L)
})

test_that("lr_test() df counts dropped design columns, not terms (audit 1.7)", {
  d <- genevzinb2_factor()
  m <- suppressWarnings(suppressMessages(evzinb(y ~ x1 + g, formula_zi = ~x1,
                                               formula_evi = ~x1, formula_pareto = ~x1,
                                               data = d, bootstrap = FALSE, verbose = FALSE)))
  res <- suppressWarnings(lr_test(m, "g"))
  # g has three levels -> two dummy columns, only in the count component
  expect_equal(res$df, 2L)
})
