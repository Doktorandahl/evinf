# round9 E.3 (audit §5.5): Poisson / ZIP competitor baselines in
# compare_models(), added for a family = "poisson" evzinb()/evinb() fit
# alongside the existing nb / zinb baselines.

fit_evzinb_poisson_fast <- function(...) {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$y[1:8] <- d$y[1:8] + 40L
  suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2, data = d, family = "poisson", n_bootstraps = 4, boot_seed = 1,
    multicore = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200),
    ...
  )))
}

test_that("poisson_comparison/zip_comparison default to TRUE only for a Poisson-family fit", {
  m_pois <- fit_evzinb_poisson_fast()
  cmp_pois <- suppressWarnings(compare_models(m_pois, multicore = FALSE))
  expect_true(all(c("poisson", "zip") %in% names(cmp_pois)))

  m_nb <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 4)
  cmp_nb <- suppressWarnings(compare_models(m_nb, multicore = FALSE))
  expect_false(any(c("poisson", "zip") %in% names(cmp_nb)))
})

test_that("poisson_comparison/zip_comparison are user-overridable either way", {
  m_nb <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 4)
  cmp <- suppressWarnings(compare_models(
    m_nb, poisson_comparison = TRUE, zip_comparison = TRUE, multicore = FALSE
  ))
  expect_true(all(c("poisson", "zip") %in% names(cmp)))

  m_pois <- fit_evzinb_poisson_fast()
  cmp2 <- suppressWarnings(compare_models(
    m_pois, poisson_comparison = FALSE, zip_comparison = FALSE, multicore = FALSE
  ))
  expect_false(any(c("poisson", "zip") %in% names(cmp2)))
})

test_that("compare_models() guards zip_comparison for evinb the same way as zinb_comparison", {
  mi <- fit_evinb_fast(n_bootstraps = 4)
  suppressWarnings(expect_message(
    comp <- compare_models(mi, poisson_comparison = TRUE), "zip_comparison"
  ))
  expect_false("zip" %in% names(comp))
  expect_true("poisson" %in% names(comp))
  expect_error(
    suppressWarnings(compare_models(mi, zip_comparison = TRUE)),
    "not available"
  )
})

test_that("predict.poissonboot()/predict.zipboot() are registered and functional", {
  m <- fit_evzinb_poisson_fast()
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))

  expect_true(is.numeric(predict(comp$poisson, type = "predicted")))
  expect_true(is.numeric(predict(comp$poisson, type = "quantile", quantile = 0.9)))
  expect_true(is.numeric(predict(comp$poisson, pred = "bootstrap_median")))

  expect_true(is.numeric(predict(comp$zip, type = "predicted")))
  expect_true(is.numeric(predict(comp$zip, type = "quantile", quantile = 0.9)))
  z_states <- predict(comp$zip, type = "states")
  expect_true(all(c("pr_zc", "pr_count") %in% names(z_states)))
  expect_true(is.numeric(predict(comp$zip, pred = "bootstrap_median")))
})

test_that("predict.*boot(pred = 'original') on poisson/zip uses the full-sample model, not a stray index", {
  m <- fit_evzinb_poisson_fast()
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))
  nd <- { data(genevzinb2, package = "evinf", envir = environment()); genevzinb2[1:5, ] }
  for (slot in c("poisson", "zip")) {
    expect_length(predict(comp[[slot]], pred = "original"), 100L)
    expect_length(predict(comp[[slot]], newdata = nd, pred = "original"), 5L)
  }
})

test_that("tidy()/glance() include poisson/zip rows generically, alpha = NA (no dispersion)", {
  m <- fit_evzinb_poisson_fast()
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))

  td <- tidy(comp)
  expect_true(all(c("poisson", "zip") %in% td$model))

  gl <- glance(comp)
  expect_true(all(c("poisson", "zip") %in% gl$model))
  expect_true(all(is.na(gl$alpha[gl$model %in% c("poisson", "zip")])))
  # glance.nbboot()'s npar is glm.nb()$rank, which (like glm()$rank) already
  # excludes the dispersion parameter -- so nb and poisson agree here. zinb's
  # npar explicitly adds +1 for theta (nrow(vcov) + 1), so zip (no +1) is one
  # fewer.
  expect_equal(gl$npar[gl$model == "poisson"], gl$npar[gl$model == "nb"])
  expect_equal(gl$npar[gl$model == "zip"], gl$npar[gl$model == "zinb"] - 1)
})

test_that("compare_fit()/oob_evaluation() work generically with the new poisson/zip slots", {
  m <- fit_evzinb_poisson_fast()
  comp <- suppressWarnings(compare_models(m, multicore = FALSE))

  cf <- compare_fit(comp)
  expect_true(all(c("poisson", "zip") %in% cf$model))

  oob <- oob_evaluation(comp)
  expect_true(all(c("poisson", "zip") %in% names(oob)))
})

test_that("winsorize/razorize build poisson/zip variants too", {
  m <- fit_evzinb_poisson_fast()
  comp <- suppressWarnings(compare_models(
    m, winsorize = TRUE, razorize = TRUE, cutoff_value = 5, multicore = FALSE
  ))
  expect_true(all(c("poisson_winsor", "zip_winsor", "poisson_razor", "zip_razor") %in% names(comp)))
})
