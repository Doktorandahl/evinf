# round9 0.1 (review §2): a fitted Pareto shape (alpha_pl) that collapses
# toward 0 makes exp(1/alpha_pl)/1/alpha_pl-based quantities (harmonic/explog
# predictions, the continuous mixture quantile, $fitted's tail summaries)
# silently return Inf or an astronomically large finite number.
# evinf_clamp_alpha_pl() itself is exercised in test-quantile-extractor.R,
# which also predates and motivated this shared helper.

test_that("predict(type = 'harmonic'/'explog') clamp near-zero alpha_pl and warn naming the count", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  n <- nobs(m)
  testthat::local_mocked_bindings(
    fitted_alpha_from_evzinb = function(object, newdata = NULL, return_data = FALSE) {
      tibble::tibble(pareto_alpha = c(0.001, rep(0.5, n - 1)))
    },
    .package = "evinf"
  )
  expect_warning(h <- predict(m, type = "harmonic"), "1 fitted Pareto alpha value")
  expect_true(all(is.finite(h)))
  expect_warning(e <- predict(m, type = "explog"), "1 fitted Pareto alpha value")
  expect_true(all(is.finite(e)))
})

test_that("evinf_control(alpha_pl_floor = ) is honoured end-to-end", {
  m <- fit_evzinb_fast(bootstrap = FALSE, control = .fast_control(alpha_pl_floor = 0.2))
  expect_equal(m$control$alpha_pl_floor, 0.2)

  n <- nobs(m)
  testthat::local_mocked_bindings(
    fitted_alpha_from_evzinb = function(object, newdata = NULL, return_data = FALSE) {
      tibble::tibble(pareto_alpha = rep(0.1, n))
    },
    .package = "evinf"
  )
  expect_warning(h <- predict(m, type = "harmonic"), "below the floor")

  prbs <- suppressWarnings(evinf:::prob_from_evzinb(m))
  cnts <- evinf:::counts_from_evzinb(m)
  expected <- prbs$pr_count * cnts$count +
    prbs$pr_pareto * m$coef$C * (1 + 0.2) / 0.2
  expect_equal(unname(h), expected, tolerance = 1e-8)
})

test_that("glance() reports the unclamped min_alpha_pl and print() notes a collapse", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_gte(glance(m)$min_alpha_pl, m$control$alpha_pl_floor)
  expect_false(grepl("collapsed", paste(utils::capture.output(print(m)), collapse = "\n")))

  m2 <- m
  m2$fitted$alpha.pl[1] <- 1e-20
  expect_equal(glance(m2)$min_alpha_pl, 1e-20)
  out <- paste(utils::capture.output(print(m2)), collapse = "\n")
  expect_true(grepl("collapsed", out))
  expect_true(grepl("alpha_pl_floor = 0.01", out))
})

test_that("harmonic prediction stays finite when the fitted Pareto shape overflows exp() to Inf", {
  # round9 0.1 follow-up, found via CI: on some platforms
  # genevzinb2_factor()'s y ~ x1 + g fit (test-model-matrix.R:30) gives an
  # extreme Beta.PL coefficient for a sparse factor level, pushing one row's
  # Pareto-shape linear predictor past ~709.78 -- alpha_pl = exp(that)
  # overflows to literal Inf. This is the *opposite* of the near-zero
  # collapse the floor above guards against, and harmonic_calc()'s old
  # (1 + alpha) / alpha formula computed Inf/Inf = NaN there. Rewritten as
  # 1/alpha + 1, which is exact for any finite alpha and gives the
  # mathematically correct limit (1, i.e. the extreme-value contribution
  # -> C) as alpha -> Inf -- no ceiling/clamp needed for this failure mode.
  m <- fit_evzinb_fast(bootstrap = FALSE)
  n <- nobs(m)
  testthat::local_mocked_bindings(
    fitted_alpha_from_evzinb = function(object, newdata = NULL, return_data = FALSE) {
      tibble::tibble(pareto_alpha = c(Inf, 1e300, exp(709), rep(1, n - 3)))
    },
    .package = "evinf"
  )
  expect_no_warning(h <- predict(m, type = "harmonic"))
  expect_true(all(is.finite(h)))

  prbs <- suppressWarnings(evinf:::prob_from_evzinb(m))
  cnts <- evinf:::counts_from_evzinb(m)
  expect_equal(unname(h[1]),
               unname(prbs$pr_count[1] * cnts$count[1] + prbs$pr_pareto[1] * m$coef$C))
})

test_that("hks with a four-covariate specification yields finite predictions despite a collapsing alpha_pl (round9 0.1, review §2)", {
  skip_on_cran()
  data(hks, package = "evinf", envir = environment())
  f_hks4 <- osvAll ~ troopLag + policeLag + militaryobserversLag + epduration
  f_hks_pareto <- ~ troopLag_log + epdur_log + brv_AllLag_log + lntpop

  m <- suppressMessages(suppressWarnings(evzinb(
    f_hks4, data = hks, formula_pareto = f_hks_pareto,
    bootstrap = FALSE, verbose = FALSE
  )))

  expect_true(is.finite(glance(m)$min_alpha_pl))
  for (ty in c("harmonic", "explog", "counts", "pareto_alpha")) {
    expect_true(all(is.finite(predict(m, type = ty))), info = ty)
  }
  expect_true(all(is.finite(m$fitted$y.hat.pl_exp.E.logy)))
  expect_true(all(is.finite(m$fitted$y.hat.pl_median)))
})
