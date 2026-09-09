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

  suppressWarnings(expect_message(comp <- compare_models(mi), "zinb_comparison"))
  expect_identical(comp$model, comp$evzinb)
  expect_false("zinb" %in% names(comp))
  expect_error(suppressWarnings(compare_models(mi, zinb_comparison = TRUE)),
               "not available")
})

test_that("predict.nbboot()/predict.zinbboot() are registered (audit 2.5)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  comp <- suppressWarnings(suppressMessages(compare_models(m)))
  expect_true(is.numeric(predict(comp$nb, pred = "bootstrap_median")))
  expect_true(is.numeric(predict(comp$zinb, pred = "bootstrap_median")))
})

test_that("predict.*boot(pred = 'original') uses the full-sample model, not a stray index", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  comp <- suppressWarnings(suppressMessages(compare_models(m)))
  nd <- { data(genevzinb2, package = "evinf", envir = environment()); genevzinb2[1:5, ] }
  for (slot in c("nb", "zinb")) {
    expect_length(predict(comp[[slot]], pred = "original"), 100L)
    expect_length(predict(comp[[slot]], newdata = nd, pred = "original"), 5L)
  }
})

test_that("compare_models() boot_refit_family() refactor is output-identical (issue 4.1)", {
  # Baseline captured from the pre-refactor code (twelve foreach blocks) for
  # evzinb(y ~ x1 + x2 + x3, genevzinb2, n_bootstraps = 5, boot_seed = 123,
  #        control = evinf_control(c.lim = c(50, 1000), init.C = 200)).
  base_rmsle_nb   <- c(2.5423001991, 2.6974881320, 2.6921071076, 2.5755132916, 2.6253084900)
  base_rmsle_zinb <- c(2.4376849746, 2.5878614556, 2.7054591325, 2.6991174880, 2.8190669965)
  base_rmse_nb    <- c(125.4124041023, 127.1601900074, 115.6458855064, 99.6642756981, 87.0973165485)
  base_rmse_zinb  <- c(119.8545366841, 108.9729871326, 147.1069088120, 103.2540030443, 133.0767432541)

  m <- fit_evzinb_fast(n_bootstraps = 5)
  comp <- suppressWarnings(suppressMessages(compare_models(m)))

  oe_rmsle <- suppressWarnings(oob_evaluation(comp, metric = "rmsle"))
  oe_rmse  <- suppressWarnings(oob_evaluation(comp, metric = "rmse"))

  expect_equal(unname(oe_rmsle$nb),   base_rmsle_nb,   tolerance = 1e-8)
  expect_equal(unname(oe_rmsle$zinb), base_rmsle_zinb, tolerance = 1e-8)
  expect_equal(unname(oe_rmse$nb),    base_rmse_nb,    tolerance = 1e-6)
  expect_equal(unname(oe_rmse$zinb),  base_rmse_zinb,  tolerance = 1e-6)
})

test_that("compare_models() winsorize + razorize still name all six slots (issue 4.1)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  comp <- suppressWarnings(suppressMessages(
    compare_models(m, winsorize = TRUE, razorize = TRUE, cutoff_value = 30)
  ))
  expect_setequal(
    setdiff(names(comp), c("model", "evzinb")),
    c("nb", "zinb", "nb_winsor", "zinb_winsor", "nb_razor", "zinb_razor")
  )
  for (s in setdiff(names(comp), c("model", "evzinb"))) {
    expect_s3_class(comp[[s]], if (grepl("zinb", s)) "zinbboot" else "nbboot")
    expect_length(comp[[s]]$bootstraps, 5L)
  }
})
