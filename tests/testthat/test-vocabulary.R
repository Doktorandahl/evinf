test_that("component aliases are accepted with a deprecation warning (audit 2.7)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)

  expect_warning(a <- coefficient_extractor(m, "nb"), "deprecated")
  b <- coefficient_extractor(m, "count")
  expect_identical(a, b)

  expect_warning(tidy(m, component = "zi"), "deprecated")
})

test_that("canonical names are used across tidy()/summary()/coefficient_extractor() (audit 2.7)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)

  expect_named(summary(m)$coefficients, c("count", "zero", "evi", "pareto"))
  expect_equal(sort(unique(coefficient_extractor(m, "all")$.component)),
               c("count", "evi", "pareto", "zero"))
})

test_that("predict(type=) accepts the new aliases and renames state columns (audit 2.7)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)

  expect_equal(predict(m, type = "zero"), predict(m, type = "zi"))
  expect_equal(predict(m, type = "evi"), predict(m, type = "evinf"))

  st <- predict(m, type = "states")
  expect_true(all(c("pr_zero", "pr_count", "pr_evi") %in% names(st)))
  expect_equal(st$pr_zero, st$pr_zc)          # deprecated duplicate kept
  expect_equal(st$pr_evi, st$pr_pareto)
})

test_that("fitted object exposes both new and deprecated component names (audit 2.7)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_equal(m$fitted$prob_evi, m$fitted$prob_pareto)
  expect_equal(colnames(m$props), c("zero", "count", "evi"))
})
