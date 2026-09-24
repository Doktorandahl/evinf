# round11 A5 (review §5): predict(m, type = "harmonic", nonsense_arg = TRUE)
# was accepted silently -- `...` was captured but never checked. This is
# exactly why A4's clamp_alpha_pl looked implemented when it wasn't: a typo'd
# or misplaced argument produced no error, no warning, and a default-argument
# result. evinf_validate_predict_args() (R/predict_distribution.R) closes
# both halves of that gap: a name that is not a known predict() argument at
# all errors, and a known argument that does nothing for the requested `type`
# warns.

test_that("predict() errors on an unknown argument, naming the valid ones for that type", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_error(
    predict(m, type = "harmonic", nonsense_arg = TRUE),
    "unknown argument.*nonsense_arg.*type = .harmonic."
  )
  # the error lists at least the universal arguments
  expect_error(predict(m, type = "harmonic", nonsense_arg = TRUE), "conf_level")
})

test_that("predict() errors on a typo'd type-specific argument name", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_error(
    predict(m, type = "exceedance", treshold = c(10, 20)),
    "unknown argument.*treshold"
  )
  expect_error(
    predict(m, type = "draws", ndraws = 100),
    "unknown argument.*ndraws"
  )
})

test_that("predict() warns (does not error) on a known-but-irrelevant argument", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_warning(
    p <- predict(m, type = "harmonic", threshold = c(10, 20)),
    "threshold.*ignored.*type = .harmonic."
  )
  # the irrelevant argument had no effect on the result
  expect_equal(unname(p), unname(suppressWarnings(predict(m, type = "harmonic"))))

  expect_warning(
    predict(m, type = "distribution", n_draws = 5),
    "n_draws.*ignored"
  )
})

test_that("predict() accepts every type's own arguments without warning", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  expect_no_warning(predict(m, type = "quantile", quantile = 0.5))
  expect_no_warning(predict(m, type = "distribution", support = 0:10, format = "matrix"))
  expect_no_warning(predict(m, type = "exceedance", threshold = c(10, 20)))
  expect_no_warning(predict(m, type = "draws", n_draws = 5, seed = 1))
  expect_no_warning(predict(m, type = "explog", clamp_alpha_pl = TRUE))
  # universal arguments never warn regardless of type
  expect_no_warning(predict(m, type = "harmonic", pred = "bootstrap_median"))
  expect_no_warning(predict(m, type = "counts", confint = TRUE, conf_level = 0.9))
})

test_that("predict.evinb() validates arguments the same way as predict.evzinb()", {
  m <- fit_evinb_fast(bootstrap = FALSE)
  expect_error(predict(m, type = "harmonic", nonsense_arg = TRUE), "unknown argument")
  expect_warning(predict(m, type = "harmonic", threshold = 10), "ignored")
})

test_that("fitted() forwards clamp_alpha_pl to predict() (round11 A4)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_equal(
    unname(suppressMessages(fitted(m, type = "explog", clamp_alpha_pl = TRUE))),
    unname(suppressMessages(predict(m, type = "explog", clamp_alpha_pl = TRUE)))
  )
})

test_that("predict_grid() and marginal_effects() take no unforwarded `...`, so nothing to guard (round11 A5)", {
  # Both build their predict()/evinf_ame_one() argument lists from their own
  # explicit formals -- no `...` -- so a typo'd argument to either is already
  # a base R "unused argument" error, not a silently-swallowed one.
  expect_null(formals(predict_grid)[["..."]])
  expect_null(formals(marginal_effects)[["..."]])
})
