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
  # Baseline captured for evzinb(y ~ x1 + x2 + x3, genevzinb2, n_bootstraps = 5,
  # boot_seed = 123, control = evinf_control(c.lim = c(50, 1000), init.C = 200)).
  # Regenerated for round 5 Part A (furrr per-element L'Ecuyer streams replace
  # %dorng%, so the bootstrap resamples - and hence these pins - changed).
  base_rmsle_nb   <- c(2.7702147669, 2.4477240167, 2.7208083737, 2.5454411241, 2.5069424574)
  base_rmsle_zinb <- c(2.6426746324, 2.5363944513, 2.6628125688, 2.6642800734, 2.5294352959)
  base_rmse_nb    <- c(118.2469492256, 103.4495487558, 100.7679581460, 117.0706907652, 98.4691091791)
  base_rmse_zinb  <- c(111.5976544278, 98.5212097513, 102.5366179051, 116.3133542114, 110.7324410291)

  m <- fit_evzinb_fast(n_bootstraps = 5)
  comp <- suppressWarnings(suppressMessages(compare_models(m)))

  oe_rmsle <- suppressWarnings(oob_evaluation(comp, metric = "rmsle"))
  oe_rmse  <- suppressWarnings(oob_evaluation(comp, metric = "rmse"))

  expect_equal(unname(oe_rmsle$nb),   base_rmsle_nb,   tolerance = 1e-8)
  expect_equal(unname(oe_rmsle$zinb), base_rmsle_zinb, tolerance = 1e-8)
  expect_equal(unname(oe_rmse$nb),    base_rmse_nb,    tolerance = 1e-6)
  expect_equal(unname(oe_rmse$zinb),  base_rmse_zinb,  tolerance = 1e-6)
})

test_that("boot_refit_one() remaps a resample onto data_razor's own rows (audit0.10 §1.5)", {
  set.seed(1)
  full_n <- 100
  keep <- sort(sample(seq_len(full_n), 60))
  boot_id <- sample(seq_len(full_n), full_n, replace = TRUE)
  expected_local <- match(boot_id[boot_id %in% keep], keep)

  data_razor <- data.frame(y = seq_along(keep), x1 = rnorm(length(keep)))
  sp <- list(type = "nb", data = data_razor,
             formulas = list(formula_nb = y ~ x1),
             has_init_theta = FALSE, init_theta = NULL,
             keep = keep, y_orig = data_razor$y)

  testthat::local_mocked_bindings(
    inner_nb = function(bootstrap, data, formulas, init_theta, y_orig) bootstrap$boot_id,
    .package = "evinf"
  )
  got <- evinf:::boot_refit_one(sp, boot_id)
  expect_equal(got, expected_local)
  expect_true(all(got >= 1 & got <= length(keep)))
  expect_false(anyNA(got))
})

test_that("compare_models(razorize = TRUE) refits on exactly the retained rows (audit0.10 §1.5)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  boot_id <- m$bootstraps[[1]]$boot_id
  y_full <- m$data$data$y
  cutoff_value <- 20
  keep <- which(y_full < sort(y_full, decreasing = TRUE)[cutoff_value])
  data_razor <- m$data$data[keep, ]

  local_id <- match(boot_id[boot_id %in% keep], keep)
  expected_fit <- suppressWarnings(MASS::glm.nb(y ~ x1 + x2 + x3, data = data_razor[local_id, ]))

  comp <- suppressWarnings(suppressMessages(
    compare_models(m, nb_comparison = TRUE, zinb_comparison = FALSE,
                   razorize = TRUE, cutoff_value = cutoff_value, multicore = FALSE)
  ))
  got_fit <- comp$nb_razor$bootstraps[[1]]

  expect_false(inherits(got_fit, "try-error"))
  expect_false(anyNA(coef(got_fit)))
  expect_equal(unname(coef(got_fit)), unname(coef(expected_fit)), tolerance = 1e-6)
})

test_that("compare_models(winsorize = TRUE) computes OOB error against the raw outcome (audit0.10 §1.5)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  comp <- suppressWarnings(suppressMessages(
    compare_models(m, nb_comparison = TRUE, zinb_comparison = FALSE,
                   winsorize = TRUE, cutoff_value = 20, multicore = FALSE)
  ))
  b <- comp$nb_winsor$bootstraps[[1]]
  expect_false(inherits(b, "try-error"))

  boot_id <- m$bootstraps[[1]]$boot_id
  y_raw <- m$data$data$y
  dv_raw <- y_raw[-boot_id]
  # The bug computed this against the winsorised outcome; assert the raw and
  # winsorised outcomes actually differ in the OOB set here, so this test is
  # not accidentally trivial.
  y_winsor_val <- sort(y_raw, decreasing = TRUE)[20]
  expect_true(any(dv_raw > y_winsor_val))

  expected_rmse <- sqrt(mean((dv_raw - b$oob_predictions)^2))
  expect_equal(b$oob_rmse, expected_rmse, tolerance = 1e-10)
})

test_that("boot_refit_one()'s try()s are silent on failure (audit0.10 §1.5)", {
  data(genevzinb2, package = "evinf", envir = environment())
  sp <- list(type = "nb", data = genevzinb2,
             formulas = list(formula_nb = y ~ nonexistent_column),
             has_init_theta = FALSE, init_theta = NULL, keep = NULL,
             y_orig = genevzinb2$y)
  boot_id <- seq_len(nrow(genevzinb2))

  msgs <- capture.output(
    fit <- evinf:::boot_refit_one(sp, boot_id),
    type = "message"
  )
  expect_true(inherits(fit, "try-error"))
  expect_length(msgs, 0)
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
