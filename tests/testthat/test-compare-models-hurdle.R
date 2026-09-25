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
