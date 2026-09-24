test_that("lr_test() reports 2 * (logLik_full - logLik_restricted) (audit 1.1)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  res <- suppressWarnings(lr_test(m, "x1"))

  expect_s3_class(res, "tbl_df")
  expect_equal(res$statistic, 2 * (res$loglik_full - res$loglik_restricted))
  expect_equal(res$prob, pchisq(res$statistic, res$df, lower.tail = FALSE))
  # restricting x1 removes one column from each of the four components
  expect_equal(res$df, 4L)
})

test_that("lr_test(bootstrap = TRUE) doubles the bootstrap statistics too (audit 1.1)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  res <- suppressWarnings(lr_test(m, "x1", bootstrap = TRUE))

  expect_named(res, c("results", "boot_results"))
  bs <- res$boot_results[["x1"]]
  ok <- !is.na(bs$statistic)
  expect_equal(bs$statistic[ok], 2 * (bs$ll_full[ok] - bs$ll_reduced[ok]))
})

test_that("lr_test() works on evinb objects (audit 1.5)", {
  mi <- fit_evinb_fast(bootstrap = FALSE)
  res <- suppressWarnings(lr_test(mi, "x1"))
  expect_s3_class(res, "tbl_df")
  # x1 is in nb, evi and pareto (no zero-inflation component)
  expect_equal(res$df, 3L)
})

test_that("lr_test() df counts dropped design columns, not terms (audit 1.7)", {
  d <- genevzinb2_factor()
  m <- suppressWarnings(suppressMessages(evzinb(y ~ x1 + g, formula_zi = ~x1,
                                               formula_evi = ~x1, formula_pareto = ~x1,
                                               data = d, bootstrap = FALSE, verbose = FALSE)))
  res <- suppressWarnings(lr_test(m, "g"))
  # g has three levels -> two dummy columns, only in the count component
  expect_equal(res$df, 2L)
})

test_that("lr_test(bootstrap = TRUE) survives a try-error replicate (audit0.10 §1.1)", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  m$bootstraps[[2]] <- try(stop("boom"), silent = TRUE)

  res <- suppressWarnings(lr_test(m, "x1", bootstrap = TRUE))

  bs <- res$boot_results[["x1"]]
  # one replicate was dropped up front by evinf_usable_bootstraps(); the
  # remaining 5 must all line up (no length mismatch / no spurious NAs from
  # misalignment between ll_full and ll_reduced).
  expect_equal(nrow(bs), 5L)
  expect_equal(res$results$n_bootstraps_used, 5L)
  expect_true(all(!is.na(bs$ll_full)))
})

test_that("lr_test() refits emit no deprecation warning (audit0.10 §1.7)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  expect_no_warning(lr_test(m, "x1"))
  expect_no_warning(lr_test(m, "x1", bootstrap = TRUE))
})

test_that("lr_test(bootstrap = TRUE) excludes degenerate replicates by default (audit0.10 §1.1)", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  m$bootstraps[[3]]$degenerate <- TRUE
  m$bootstraps[[3]]$degenerate_reason <- "injected for test"

  res_default <- suppressWarnings(lr_test(m, "x1", bootstrap = TRUE))
  expect_equal(res_default$results$n_bootstraps_used, 5L)
  expect_equal(nrow(res_default$boot_results[["x1"]]), 5L)

  res_kept <- suppressWarnings(
    lr_test(m, "x1", bootstrap = TRUE, exclude_degenerate = FALSE))
  expect_equal(res_kept$results$n_bootstraps_used, 6L)
  expect_equal(nrow(res_kept$boot_results[["x1"]]), 6L)
})

# round10 0.1 (review §1): lr_refit_restricted() dropped weights, family and
# offset() terms, so a weighted or offset fit's LR statistic was silently
# wrong, and a Poisson or hurdle fit's restricted refit errored or compared
# against the wrong nested model.

test_that("a weighted fit's LR statistic equals the row-duplicated one (round10 0.1)", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(21)
  w_int <- sample(1:3, nrow(genevzinb2), replace = TRUE)
  d <- genevzinb2
  d$w <- w_int
  ctrl <- evinf_control(c.lim = c(50, 1000), init.C = 200,
                        max.diff.par = 1e-6, max.no.em.steps = 3000)

  m_w <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d, weights = w,
    bootstrap = FALSE, verbose = FALSE, control = ctrl
  )))
  res_w <- suppressWarnings(lr_test(m_w, "x1"))
  expect_gte(res_w$statistic, 0)

  idx <- rep(seq_len(nrow(d)), w_int)
  d_dup <- d[idx, ]
  m_dup <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = d_dup,
    bootstrap = FALSE, verbose = FALSE, control = ctrl
  )))
  res_dup <- suppressWarnings(lr_test(m_dup, "x1"))

  expect_equal(res_w$statistic, res_dup$statistic, tolerance = 1e-6)
})

test_that("an offset model's LR statistic keeps the offset in the restricted refit (round10 0.1)", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(22)
  d <- genevzinb2
  d$ex <- runif(nrow(d), 0.5, 2)
  ctrl <- evinf_control(c.lim = c(50, 1000), init.C = 200,
                        max.diff.par = 1e-6, max.no.em.steps = 3000)

  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + offset(log(ex)), data = d,
    bootstrap = FALSE, verbose = FALSE, control = ctrl
  )))
  res <- suppressWarnings(lr_test(m, "x2"))

  # Ground truth: the same warm-started control lr_test() uses internally
  # (restricted_control()), but with the reduced formula typed by hand
  # instead of produced by formula_var_remover() -- so this is an
  # independent check that the offset survived the reduction.
  lr_obj <- structure(list(control = m$control, coef = m$coef, family = m$family),
                      class = class(m))
  reduced_formulas <- list(nb = y ~ x1 + offset(log(ex)), zi = y ~ x1,
                           evinf = y ~ x1, pareto = y ~ x1)
  ctrl_ws <- restricted_control(lr_obj, reduced_formulas, m$data$data)
  m_hand <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + offset(log(ex)), data = d,
    bootstrap = FALSE, verbose = FALSE, control = ctrl_ws
  )))

  expect_equal(res$df, 4L)
  expect_equal(res$loglik_restricted, m_hand$log.lik, tolerance = 1e-6)
  expect_equal(res$statistic, 2 * (m$log.lik - m_hand$log.lik), tolerance = 1e-6)
})

test_that("lr_test() runs on a Poisson-count fit, with and without bootstrap (round10 0.1)", {
  m <- suppressWarnings(fit_evzinb_fast(family = "poisson", n_bootstraps = 5))
  res <- suppressWarnings(lr_test(m, "x1"))
  expect_s3_class(res, "tbl_df")
  expect_true(is.finite(res$statistic))

  res_boot <- suppressWarnings(lr_test(m, "x1", bootstrap = TRUE))
  expect_named(res_boot, c("results", "boot_results"))
})

test_that("lr_test()'s restricted refit is itself a hurdle model (round10 0.1)", {
  m <- fit_evzinb_fast(family = evinf_family(zero = "hurdle"), n_bootstraps = 5)
  res <- suppressWarnings(lr_test(m, "x1"))
  expect_true(is.finite(res$statistic))

  reduced <- formula_var_remover(m$formulas, "x1", m$data$data)$formulas
  lr_obj <- structure(list(control = m$control, coef = m$coef, family = m$family,
                           weights = m$weights), class = class(m))
  refit <- lr_refit_restricted(reduced, m$data$data, lr_obj, m$data$data,
                               weights = lr_obj$weights)
  expect_equal(refit$family$zero, "hurdle")
  expect_equal(res$loglik_restricted, refit$log.lik, tolerance = 1e-8)

  res_boot <- suppressWarnings(lr_test(m, "x1", bootstrap = TRUE))
  expect_named(res_boot, c("results", "boot_results"))
})

test_that("formula_var_remover() keeps a no-intercept specification (round10 0.1)", {
  data(genevzinb2, package = "evinf", envir = environment())
  formulas <- list(
    formula_nb = y ~ x1 + x2 - 1,
    formula_zi = NULL,
    formula_evi = y ~ x1 + x2,
    formula_pareto = y ~ x1 + x2
  )
  out <- formula_var_remover(formulas, "x2", genevzinb2)
  expect_equal(attr(stats::terms(out$formulas$nb), "intercept"), 0L)
  expect_false("x2" %in% all.vars(out$formulas$nb))
})
