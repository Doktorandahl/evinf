# round12 B1 (dev/review_round11.md §4): a hurdle competitor in
# compare_models(), added alongside the existing NB/ZINB/Poisson/ZIP
# baselines so the table can also answer "what does the extreme-value state
# buy over the like-for-like hurdle", not just over the standard toolkit
# model.

fit_evzinb_hurdle_fast <- function(...) {
  data(genevzinb2, package = "evinf", envir = environment())
  suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2,
    family = evinf_family(zero = "hurdle"), n_bootstraps = 4, boot_seed = 1,
    multicore = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200),
    ...
  )))
}

fit_evzinb_hurdle_poisson_fast <- function(...) {
  data(genevzinb2, package = "evinf", envir = environment())
  suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2,
    family = evinf_family(count = "poisson", zero = "hurdle"),
    n_bootstraps = 4, boot_seed = 1, multicore = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200),
    ...
  )))
}

# round11 B2: this file is dominated by bootstrap-heavy compare_models()
# fits; kept fast on CI (NOT_CRAN=true) but skipped on CRAN's own check to
# stay within its check-time budget.
testthat::skip_on_cran()

test_that("hurdle_comparison defaults to TRUE only for a hurdle-family fit", {
  m_hurdle <- fit_evzinb_hurdle_fast()
  cmp_hurdle <- suppressWarnings(compare_models(m_hurdle, multicore = FALSE))
  expect_true("hurdle" %in% names(cmp_hurdle))

  m_mix <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 4)
  cmp_mix <- suppressWarnings(compare_models(m_mix, multicore = FALSE))
  expect_false("hurdle" %in% names(cmp_mix))
})

test_that("hurdle_comparison is user-overridable either way", {
  m_mix <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 4)
  cmp <- suppressWarnings(compare_models(m_mix, hurdle_comparison = TRUE, multicore = FALSE))
  expect_true("hurdle" %in% names(cmp))

  m_hurdle <- fit_evzinb_hurdle_fast()
  cmp2 <- suppressWarnings(compare_models(m_hurdle, hurdle_comparison = FALSE, multicore = FALSE))
  expect_false("hurdle" %in% names(cmp2))
})

test_that("compare_models() guards hurdle_comparison for evinb the same way as zinb_comparison", {
  mi <- fit_evinb_fast(n_bootstraps = 4)
  suppressWarnings(expect_message(
    comp <- compare_models(mi, nb_comparison = TRUE), "hurdle_comparison"
  ))
  expect_false("hurdle" %in% names(comp))
  expect_error(
    suppressWarnings(compare_models(mi, hurdle_comparison = TRUE)),
    "not available"
  )
})

test_that("the hurdle competitor's count dist matches the fitted model's own", {
  m_nb <- fit_evzinb_hurdle_fast()
  cmp_nb <- suppressWarnings(compare_models(m_nb, multicore = FALSE))
  expect_equal(cmp_nb$hurdle$full_run$dist$count, "negbin")

  m_pois <- fit_evzinb_hurdle_poisson_fast()
  cmp_pois <- suppressWarnings(compare_models(m_pois, multicore = FALSE))
  expect_equal(cmp_pois$hurdle$full_run$dist$count, "poisson")
})

test_that("tidy()/glance()/predict()/coefficient_extractor() are registered for hurdleboot", {
  m <- fit_evzinb_hurdle_fast()
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))

  td <- tidy(comp$hurdle)
  expect_true(all(c("count", "zero") %in% td$y.level))
  expect_equal(nrow(td), 8L)  # 4 count + 4 zero terms

  gl <- glance(comp$hurdle)
  expect_equal(nrow(gl), 1L)
  expect_true(is.numeric(gl$alpha) && !is.na(gl$alpha))  # negbin count -> a dispersion

  ce <- coefficient_extractor(comp$hurdle)
  expect_true(all(c("count", "zero") %in% ce$.component))

  expect_true(is.numeric(predict(comp$hurdle, type = "predicted")))
  expect_true(is.numeric(predict(comp$hurdle, type = "counts")))
  expect_true(is.numeric(predict(comp$hurdle, type = "quantile", quantile = 0.5)))
  expect_true(is.numeric(predict(comp$hurdle, pred = "bootstrap_median")))
  states <- predict(comp$hurdle, type = "states")
  expect_true(all(c("pr_zc", "pr_count") %in% names(states)))
  expect_equal(states$pr_zc + states$pr_count, rep(1, nrow(states)), tolerance = 1e-8)
})

# round13 C1: pscl::hurdle()'s own zero-hurdle coefficients model P(Y > 0);
# tidy()/coefficient_extractor() on $hurdle now negate the zero component
# onto evinf's own P(zero state) scale, so a modelsummary()-style table's
# zero block reads directly against every other competitor's (see C.4's
# amended NEWS bullet). The count component, and predict() (already
# normalised in round12 B1), are untouched.
test_that("tidy(cmp$hurdle, component = 'zero') is the exact negative of pscl::hurdle()'s own coefficients (round13 C1)", {
  m <- fit_evzinb_hurdle_fast()
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))
  fit <- comp$hurdle$full_run
  pscl_coef <- coef(fit)

  zero_td <- tidy(comp$hurdle, component = "zero")
  zero_pscl <- unname(pscl_coef[grep("^zero_", names(pscl_coef))])
  expect_equal(unname(zero_td$estimate), -zero_pscl)

  count_td <- tidy(comp$hurdle, component = "count")
  count_pscl <- unname(pscl_coef[grep("^count_", names(pscl_coef))])
  expect_equal(unname(count_td$estimate), count_pscl)

  ce <- coefficient_extractor(comp$hurdle, component = "zero")
  raw_zero_boot <- purrr::map(comp$hurdle$bootstraps, "coefficients") %>% purrr::map("zero") %>% dplyr::bind_rows()
  expect_equal(as.data.frame(ce), as.data.frame(-raw_zero_boot))

  # CI bounds are the swapped negatives of the un-negated (raw pscl-scale) draws.
  zero_ci <- tidy(comp$hurdle, component = "zero", confint = "bootstrapped")
  qs <- c(0.025, 0.975)
  raw_low  <- vapply(raw_zero_boot, stats::quantile, numeric(1), probs = qs[1])
  raw_high <- vapply(raw_zero_boot, stats::quantile, numeric(1), probs = qs[2])
  expect_equal(zero_ci$conf.low, unname(-raw_high))
  expect_equal(zero_ci$conf.high, unname(-raw_low))
})

test_that("predict.hurdleboot(type = 'zi') is P(zero state), matching predict.zinbboot()'s direction", {
  # round12 B2: pscl::hurdle()'s own zero-hurdle component models P(Y > 0) --
  # the opposite of predict(type = 'zi')'s documented meaning everywhere else
  # in the package -- so prob_from_hurdle() must flip it back, not reuse
  # prob_from_znb()'s raw column labels unchanged.
  m <- fit_evzinb_hurdle_fast()
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))
  hrd <- comp$hurdle$full_run

  p_y0 <- unname(predict(hrd, type = "prob")[, 1])
  expect_equal(predict(comp$hurdle, type = "zi"), p_y0, tolerance = 1e-8)
  expect_equal(predict(comp$hurdle, type = "count_state"), 1 - p_y0, tolerance = 1e-8)
})

test_that("predict.*boot(pred = 'original') on hurdle uses the full-sample model, not a stray index", {
  m <- fit_evzinb_hurdle_fast()
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))
  nd <- { data(genevzinb2, package = "evinf", envir = environment()); genevzinb2[1:5, ] }
  expect_length(predict(comp$hurdle, pred = "original"), 100L)
  expect_length(predict(comp$hurdle, newdata = nd, pred = "original"), 5L)
})

test_that("tidy()/glance() include the hurdle row generically", {
  m <- fit_evzinb_hurdle_fast()
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))

  td <- tidy(comp)
  expect_true("hurdle" %in% td$model)

  gl <- glance(comp)
  expect_true("hurdle" %in% gl$model)
})

test_that("compare_fit()/oob_evaluation() work generically with the new hurdle slot", {
  m <- fit_evzinb_hurdle_fast()
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))

  cf <- compare_fit(comp)
  expect_true("hurdle" %in% cf$model)

  oob <- oob_evaluation(comp)
  expect_true("hurdle" %in% names(oob))
})

test_that("winsorize/razorize build a hurdle variant too", {
  m <- fit_evzinb_hurdle_fast()
  comp <- suppressWarnings(compare_models(
    m, winsorize = TRUE, razorize = TRUE, cutoff_value = 5, multicore = FALSE
  ))
  expect_true(all(c("hurdle_winsor", "hurdle_razor") %in% names(comp)))
})

test_that("a weighted hurdle fit's competitor is weighted, matching a pscl::hurdle() fit on row-duplicated data", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(41)
  w_int <- sample(1:3, nrow(genevzinb2), replace = TRUE)
  d <- genevzinb2
  d$w <- w_int
  ctrl <- evinf_control(c.lim = c(50, 1000), init.C = 200)

  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d, weights = w,
    family = evinf_family(zero = "hurdle"),
    n_bootstraps = 3, boot_seed = 1, multicore = FALSE, verbose = FALSE, control = ctrl
  )))
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))

  expect_equal(sum(comp$hurdle$full_run$weights), sum(w_int))

  idx <- rep(seq_len(nrow(d)), w_int)
  d_dup <- d[idx, ]
  hand <- pscl::hurdle(y ~ x1 + x2 + x3 | x1 + x2 + x3, data = d_dup,
                       dist = "negbin", zero.dist = "binomial")
  expect_equal(unname(coef(comp$hurdle$full_run)), unname(coef(hand)), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(comp$hurdle$full_run)), as.numeric(logLik(hand)),
               tolerance = 1e-4)
})

# round13 B1/B3: every fixture above uses y ~ x1 + x2 + x3 with no separate
# formula_zi, so the count and zero RHSs are always identical -- exactly why
# znb_formula_extractor()'s list(nb=, zi=) branch (hit only when they
# differ) went unexercised and its `[,-1]` intercept-stripping bug shipped
# alongside predict.hurdleboot() itself. This fixture is asymmetric.
test_that("predict() on an asymmetric-formula hurdle competitor works for every type (round13 B1/B3)", {
  data(genevzinb2, package = "evinf", envir = environment())
  m <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + x2, formula_zi = ~x1, data = genevzinb2,
    family = evinf_family(zero = "hurdle"),
    control = evinf_control(c.lim = c(50, 1000), init.C = 200),
    n_bootstraps = 4, boot_seed = 1, multicore = FALSE, verbose = FALSE
  )))
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))
  fit <- comp$hurdle$full_run

  for (pr in c("original", "bootstrap_median")) {
    expect_true(all(is.finite(predict(comp$hurdle, type = "predicted", pred = pr))))
    expect_true(all(is.finite(predict(comp$hurdle, type = "counts", pred = pr))))
    expect_true(all(is.finite(predict(comp$hurdle, type = "zi", pred = pr))))
    expect_true(all(is.finite(predict(comp$hurdle, type = "count_state", pred = pr))))
    expect_true(all(is.finite(unlist(predict(comp$hurdle, type = "states", pred = pr)))))
    expect_true(all(is.finite(predict(comp$hurdle, type = "quantile", quantile = 0.5, pred = pr))))
    expect_true(all(is.finite(unlist(predict(comp$hurdle, type = "all", quantile = 0.5, pred = pr)))))
  }

  states <- predict(comp$hurdle, type = "states")
  expect_equal(states$pr_zc + states$pr_count, rep(1, nrow(states)), tolerance = 1e-8)

  # pscl::hurdle()'s own zero-hurdle coefficients model P(Y > 0); prob_from_hurdle()
  # flips them back to P(zero state), so the hand-computed check is negated
  # relative to the raw plogis(), unlike the ZINB/ZIP case (round12 B1).
  hand_pr_y0 <- as.numeric(stats::plogis(model.matrix(~x1, genevzinb2) %*% fit$coefficients$zero))
  expect_equal(states$pr_zc, 1 - hand_pr_y0, tolerance = 1e-8)

  expect_equal(
    as.numeric(predict(comp$hurdle, type = "counts")),
    as.numeric(predict(fit, type = "count", newdata = genevzinb2))
  )

  # brute-force inversion of F(y) = p0 + (1 - p0)(F_c(y) - f_c(0))/(1 - f_c(0))
  cnts <- as.numeric(predict(comp$hurdle, type = "counts"))
  q0 <- stats::dnbinom(0, mu = cnts, size = fit$theta)
  for (qq in c(0.05, 0.5, 0.9)) {
    target <- q0 + (1 - q0) * pmax(qq - states$pr_zc, 0) / (1 - states$pr_zc)
    want <- ifelse(qq - states$pr_zc > 0, stats::qnbinom(target, mu = cnts, size = fit$theta), 0)
    expect_equal(as.numeric(predict(comp$hurdle, type = "quantile", quantile = qq)), want)
  }
})

test_that("an offset model's hurdle competitor keeps the offset", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(42)
  d <- genevzinb2
  d$ex <- runif(nrow(d), 0.5, 2)

  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + offset(log(ex)), data = d,
    family = evinf_family(zero = "hurdle"),
    n_bootstraps = 3, boot_seed = 1, multicore = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200)
  )))
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))

  expect_true(!is.null(attr(stats::terms(comp$hurdle$full_run$terms$count), "offset")))
  # formula_zi (and so the hurdle competitor's zero part) inherits formula_nb
  # with any offset() stripped, not carried along -- an offset applies only
  # where written explicitly (evinf_component_formula(), R/design_matrix.R).
  hand <- pscl::hurdle(y ~ x1 + x2 + offset(log(ex)) | x1 + x2,
                       data = d, dist = "negbin", zero.dist = "binomial")
  expect_equal(unname(coef(comp$hurdle$full_run)), unname(coef(hand)), tolerance = 1e-6)
})
