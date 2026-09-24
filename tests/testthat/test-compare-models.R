# round11 B2: dominated by bootstrap-heavy compare_models() fits; kept fast
# on CI (NOT_CRAN=true) but skipped on CRAN's own check-time budget.
testthat::skip_on_cran()

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

test_that("compare_models() boot_refit_family()'s OOB metrics match a fresh in-session refit on the same boot_ids (audit0.10 §2.2, E.2)", {
  # Was pinned to hard-coded literals (issue 4.1's refactor-identity check),
  # which drifted in the 4th significant figure under pscl 1.5.9. Compares
  # against a fresh call to the same underlying inner_nb()/inner_zinb()
  # helpers instead, so the test tracks whatever pscl/MASS actually produce
  # rather than a frozen snapshot -- this still catches a boot_refit_family()
  # regression (wrong boot_id/spec paired up, wrong ordering) without being
  # sensitive to point-release numerical drift in either dependency.
  m <- fit_evzinb_fast(n_bootstraps = 5)
  comp <- suppressWarnings(suppressMessages(compare_models(m)))

  dv_f <- all.vars(m$formulas$formula_nb)[1]
  deparse_rhs <- function(f) {
    if (is.null(f)) return("1")
    rhs <- if (length(f) == 3L) f[[3]] else f[[2]]
    paste(deparse(rhs), collapse = " ")
  }
  f_zinb <- stats::as.formula(paste(dv_f, "~", deparse_rhs(m$formulas$formula_nb),
                                    "|", deparse_rhs(m$formulas$formula_zi)))
  y_orig <- dplyr::pull(m$data$data[dv_f])

  # init_theta must be *omitted*, not passed as NULL: MASS::glm.nb()'s
  # init.theta has no default and checks missing(init.theta) internally, so
  # an explicit NULL takes the wrong branch and glm.nb() errors. This
  # mirrors boot_refit_one()'s own has_init_theta branching -- compare_models()
  # was called above with its default init_theta = NULL, i.e. "no init_theta".
  fresh_nb <- fresh_zinb <- vector("list", length(m$bootstraps))
  for (k in seq_along(m$bootstraps)) {
    bid <- m$bootstraps[[k]]$boot_id
    fresh_nb[[k]]   <- suppressWarnings(evinf:::inner_nb(
      list(boot_id = bid), m$data$data, m$formulas, y_orig = y_orig))
    fresh_zinb[[k]] <- suppressWarnings(evinf:::inner_zinb(
      list(boot_id = bid), m$data$data, m$formulas, f_zinb, y_orig))
  }

  oe_rmsle <- suppressWarnings(oob_evaluation(comp, metric = "rmsle"))
  oe_rmse  <- suppressWarnings(oob_evaluation(comp, metric = "rmse"))

  extract <- function(fits, field) {
    vapply(fits, function(f) {
      if (inherits(f, "try-error") || is.null(f[[field]])) NA_real_ else f[[field]]
    }, numeric(1))
  }

  expect_equal(unname(oe_rmsle$nb),   extract(fresh_nb, "oob_rmsle"))
  expect_equal(unname(oe_rmsle$zinb), extract(fresh_zinb, "oob_rmsle"))
  expect_equal(unname(oe_rmse$nb),    extract(fresh_nb, "oob_rmse"))
  expect_equal(unname(oe_rmse$zinb),  extract(fresh_zinb, "oob_rmse"))
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
    inner_nb = function(bootstrap, data, formulas, init_theta, y_orig, weights_col = NULL) bootstrap$boot_id,
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

test_that("inner_nb()/inner_zinb() fail cleanly on an empty boot_id (round8 0.7, review §7)", {
  # data[-integer(0), ] returns zero rows, not all of them -- a resample where
  # every drawn row failed to survive razorising. Guard against it explicitly
  # rather than silently comparing OOB error on 0 rows.
  data(genevzinb2, package = "evinf", envir = environment())
  b_empty <- list(boot_id = integer(0))

  r_nb <- evinf:::inner_nb(b_empty, genevzinb2,
                           formulas = list(formula_nb = y ~ x1),
                           init_theta = NULL, y_orig = genevzinb2$y)
  expect_true(inherits(r_nb, "try-error"))
  expect_match(conditionMessage(attr(r_nb, "condition")), "boot_id is empty")

  f_zinb <- y ~ x1 | x1
  r_zinb <- evinf:::inner_zinb(b_empty, genevzinb2,
                               formulas = list(formula_nb = y ~ x1),
                               f_zinb = f_zinb, y_orig = genevzinb2$y)
  expect_true(inherits(r_zinb, "try-error"))
  expect_match(conditionMessage(attr(r_zinb, "condition")), "boot_id is empty")
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

# round10 0.2 (review §2): every competitor fit (nb / zinb / poisson / zip,
# including winsorised/razorised variants and every bootstrap refit) was
# fitted unweighted, so a weighted evzinb() fit's AIC/BIC comparisons were
# between a weighted and an unweighted likelihood. nobs.glm() itself counts
# nonzero-weight rows (base R's convention), not sum(weights), so these
# tests check the weights that actually reached each fit via weights()/
# $weights and via reproducing a fit on the row-duplicated data.

test_that("compare_models() competitor fits are weighted (round10 0.2)", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(31)
  w_int <- sample(1:3, nrow(genevzinb2), replace = TRUE)
  d <- genevzinb2
  d$w <- w_int
  ctrl <- evinf_control(c.lim = c(50, 1000), init.C = 200)

  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d, weights = w,
    n_bootstraps = 3, boot_seed = 1, multicore = FALSE, verbose = FALSE, control = ctrl
  )))
  comp <- suppressWarnings(suppressMessages(
    compare_models(m, poisson_comparison = TRUE, zip_comparison = TRUE)
  ))

  expect_equal(sum(stats::weights(comp$nb$full_run)), sum(w_int))
  expect_equal(sum(comp$zinb$full_run$weights), sum(w_int))
  expect_equal(sum(stats::weights(comp$poisson$full_run)), sum(w_int))
  expect_equal(sum(comp$zip$full_run$weights), sum(w_int))

  # The NB baseline -- a well-behaved IRLS fit -- reproduces a fit on the
  # row-duplicated data essentially exactly. (zeroinfl()'s BFGS optimiser can
  # land on a different local optimum for the same weighted-vs-duplicated
  # likelihood depending on its default starting values, so that comparison
  # is not made here -- see the commit message for a worked example.)
  idx <- rep(seq_len(nrow(d)), w_int)
  d_dup <- d[idx, ]
  m_dup <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d_dup,
    n_bootstraps = 3, boot_seed = 1, multicore = FALSE, verbose = FALSE, control = ctrl
  )))
  comp_dup <- suppressWarnings(suppressMessages(compare_models(m_dup)))
  expect_equal(unname(coef(comp$nb$full_run)), unname(coef(comp_dup$nb$full_run)),
               tolerance = 1e-6)
  expect_equal(as.numeric(logLik(comp$nb$full_run)), as.numeric(logLik(comp_dup$nb$full_run)),
               tolerance = 1e-6)
})

test_that("compare_models() carries weights into winsorised/razorised competitor fits (round10 0.2)", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(32)
  w_int <- sample(1:3, nrow(genevzinb2), replace = TRUE)
  d <- genevzinb2
  d$w <- w_int
  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d, weights = w,
    n_bootstraps = 3, boot_seed = 1, multicore = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200)
  )))
  comp <- suppressWarnings(suppressMessages(
    compare_models(m, winsorize = TRUE, razorize = TRUE, cutoff_value = 5)
  ))
  expect_equal(sum(stats::weights(comp$nb_winsor$full_run)), sum(w_int))
  keep <- which(d$y < sort(d$y, decreasing = TRUE)[5])
  expect_equal(sum(stats::weights(comp$nb_razor$full_run)), sum(w_int[keep]))
})

test_that("compare_models() carries weights into bootstrap competitor refits (round10 0.2)", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(33)
  w_int <- sample(1:3, nrow(genevzinb2), replace = TRUE)
  d <- genevzinb2
  d$w <- w_int
  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d, weights = w,
    n_bootstraps = 3, boot_seed = 1, multicore = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200)
  )))
  comp <- suppressWarnings(suppressMessages(compare_models(m)))

  b1 <- m$bootstraps[[1]]
  data_ib <- d[b1$boot_id, ]
  hand <- MASS::glm.nb(y ~ x1 + x2 + x3, data = data_ib, weights = w)
  expect_equal(unname(comp$nb$bootstraps[[1]]$fit_stats["logLik"]),
               as.numeric(logLik(hand)), tolerance = 1e-6)
})
