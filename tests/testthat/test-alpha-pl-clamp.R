# round9 0.1 (review §2): a fitted Pareto shape (alpha_pl) that collapses
# toward 0 makes exp(1/alpha_pl)/1/alpha_pl-based quantities (harmonic/explog
# predictions, the continuous mixture quantile, $fitted's tail summaries)
# silently return Inf or an astronomically large finite number.
# evinf_clamp_alpha_pl() itself is exercised in test-quantile-extractor.R,
# which also predates and motivated this shared helper.

# round11 B2: includes hks fits and a bootstrap confint check; kept fast on
# CI (NOT_CRAN=true) but skipped on CRAN's own check-time budget.
testthat::skip_on_cran()

test_that("predict(type = 'harmonic') clamps near-zero alpha_pl (via alpha_pl_floor) and warns naming the count", {
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
})

# round11 A4: type = "explog" is exempt from alpha_pl_floor (unlike harmonic),
# so clamp_alpha_pl = FALSE (the default) uses alpha_pl as-is -- at 0.001 this
# genuinely overflows to Inf, and the warning names clamp_alpha_pl as the fix.
test_that("predict(type = 'explog') is unclamped by default and may be Inf; clamp_alpha_pl fixes it (round11 A4)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  n <- nobs(m)
  testthat::local_mocked_bindings(
    fitted_alpha_from_evzinb = function(object, newdata = NULL, return_data = FALSE) {
      tibble::tibble(pareto_alpha = c(0.001, rep(0.5, n - 1)))
    },
    .package = "evinf"
  )
  expect_warning(e <- predict(m, type = "explog"),
                "1 fitted Pareto alpha value.*clamp_alpha_pl = TRUE")
  expect_true(is.infinite(e[1]))
  expect_equal(attr(e, "clamp_alpha_pl"), FALSE)

  expect_no_warning(e_clamped <- suppressMessages(
    predict(m, type = "explog", clamp_alpha_pl = TRUE)))
  expect_true(all(is.finite(e_clamped)))
  expect_equal(attr(e_clamped, "clamp_alpha_pl"), 0.1)
})

# round10 0.6 / round11 A4: explog is exempt from alpha_pl_floor entirely
# (unlike harmonic) -- warn well below where exp(1/alpha_pl) is still finite
# but not sane, regardless of what alpha_pl_floor is set to.

test_that("predict(type = 'explog') warns below 0.1 regardless of alpha_pl_floor (round11 A4)", {
  m <- fit_evzinb_fast(bootstrap = FALSE, control = .fast_control(alpha_pl_floor = 0.001))
  n <- nobs(m)
  testthat::local_mocked_bindings(
    fitted_alpha_from_evzinb = function(object, newdata = NULL, return_data = FALSE) {
      tibble::tibble(pareto_alpha = c(0.05, rep(0.5, n - 1)))
    },
    .package = "evinf"
  )
  # 0.05 is above the (looser) 0.001 floor, so no clamp warning fires -- only
  # the "effectively undefined" one, since explog is unusable there regardless.
  expect_warning(e <- predict(m, type = "explog"), "1 fitted Pareto alpha value.*below 0\\.1")
  expect_true(all(is.finite(e)))
})

test_that("predict(type = 'explog') does not warn when every alpha_pl is >= 0.1", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  n <- nobs(m)
  testthat::local_mocked_bindings(
    fitted_alpha_from_evzinb = function(object, newdata = NULL, return_data = FALSE) {
      tibble::tibble(pareto_alpha = rep(0.5, n))
    },
    .package = "evinf"
  )
  expect_no_warning(predict(m, type = "explog"))
})

test_that("the legacy $fitted$y.hat.pl_* fields are gone (round10 0.6)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  removed <- c("y.hat.pl_exp.E.logy", "y.hat.pl_E.inv.y", "y.hat.pl_median", "y.hat.pl_mean")
  expect_false(any(removed %in% names(m$fitted)))
  # fitted() still works, routed through predict() (unaffected by the removal).
  expect_equal(unname(fitted(m)), unname(predict(m, type = "harmonic")))
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

# round11 A4: the exact reproduction from the review -- a plain hks fit whose
# fitted alpha_pl genuinely collapses (min ~= 4e-32, not a mocked value).
test_that("clamp_alpha_pl fixes the hks explog blowup (round11 A4)", {
  skip_on_cran()
  data(hks, package = "evinf", envir = environment())
  f_hks4 <- osvAll ~ troopLag + policeLag + militaryobserversLag + brv_AllLag_log

  m <- suppressMessages(evzinb(
    f_hks4, data = hks, n_bootstraps = 2, multicore = FALSE,
    bootstrap = TRUE, verbose = FALSE
  ))
  expect_lt(glance(m)$min_alpha_pl, 0.1)

  expect_warning(p_unclamped <- predict(m, type = "explog"), "clamp_alpha_pl")
  expect_equal(attr(p_unclamped, "clamp_alpha_pl"), FALSE)

  expect_no_warning(p_clamped <- suppressMessages(
    predict(m, type = "explog", clamp_alpha_pl = TRUE)))
  expect_lte(max(p_clamped), m$coef$C * exp(10))
  expect_equal(attr(p_clamped, "clamp_alpha_pl"), 0.1)

  # a positive number clamps there instead of at 0.1
  p_numeric <- suppressMessages(predict(m, type = "explog", clamp_alpha_pl = 0.5))
  expect_equal(attr(p_numeric, "clamp_alpha_pl"), 0.5)
  expect_lte(max(p_numeric), m$coef$C * exp(1 / 0.5))

  # the harmonic path is unaffected by clamp_alpha_pl (A5 warns that the
  # argument is irrelevant for this type, but the values are unchanged)
  h1 <- suppressWarnings(predict(m, type = "harmonic"))
  expect_warning(h2 <- predict(m, type = "harmonic", clamp_alpha_pl = TRUE),
                "ignored")
  expect_equal(h1, h2)
})

test_that("confint = TRUE clamps every bootstrap replicate alike (round11 A4)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  n <- nobs(m)
  # every replicate's alpha_pl (mocking ignores which object is passed) is
  # forced below 0.1, so an unclamped explog would be enormous everywhere.
  testthat::local_mocked_bindings(
    fitted_alpha_from_evzinb = function(object, newdata = NULL, return_data = FALSE) {
      tibble::tibble(pareto_alpha = rep(0.001, n))
    },
    .package = "evinf"
  )
  ci <- suppressMessages(suppressWarnings(predict(
    m, type = "explog", clamp_alpha_pl = TRUE, confint = TRUE,
    return_bootstraps = TRUE
  )))
  # every replicate's alpha_pl was forced to 0.001 (unclamped: exp(1000) =
  # Inf); clamped to 0.1, explog is finite and bounded by that replicate's
  # own C * exp(10).
  boot_max <- vapply(ci$bootstraps, function(b) max(b$explog), numeric(1))
  expect_true(all(is.finite(boot_max)))
  boot_C <- vapply(evinf_usable_bootstraps(m, TRUE), function(b) b$coef$C, numeric(1))
  expect_true(all(boot_max <= boot_C * exp(10) * 1.01))
  expect_equal(attr(ci$ci, "clamp_alpha_pl"), 0.1)
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
})
