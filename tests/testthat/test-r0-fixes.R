# Phase-1 review fixes (dev/review_round1.md, items R0.1-R0.8).

test_that("block column may have NAs where the formula variables do not (R0.1)", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  set.seed(1)
  d$grp <- sample(1:20, nrow(d), replace = TRUE)
  d$grp[3] <- NA_integer_

  m <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + x2 + x3, data = d, block = "grp",
    n_bootstraps = 4, boot_seed = 1, multicore = FALSE, verbose = FALSE
  )))
  expect_equal(length(m$block_vec), nrow(m$data$data))
  expect_equal(nrow(m$data$data), sum(!is.na(d$grp)))
  expect_length(m$bootstraps, 4L)

  expect_error(evzinb(y ~ x1, data = d, block = 1:5, bootstrap = FALSE),
               "single string")
})

test_that("lr_test() refits use the full model's control settings (R0.2)", {
  m <- suppressMessages(evzinb(
    y ~ x1 + x2 + x3,
    data = { data(genevzinb2, package = "evinf", envir = environment()); genevzinb2 },
    n_bootstraps = 5, boot_seed = 1, verbose = FALSE,
    control = evinf_control(c.lim = c(20, 500))
  ))
  reduced <- formula_var_remover(m$formulas, "x1", m$data$data)$formulas
  args <- restricted_fit_args(m, reduced, m$data$data)
  expect_equal(args$c.lim, c(20, 500))
  expect_equal(args$init.C, as.numeric(m$coef$C))
  expect_equal(args$init.Alpha.NB, as.numeric(m$coef$Alpha.NB))

  res <- suppressWarnings(suppressMessages(lr_test(m, "x1")))
  expect_gte(res$statistic, 0)
})

test_that("log-likelihood recompute is silent and recorded (R0.3)", {
  m <- expect_no_warning(suppressMessages(fit_evzinb_fast(bootstrap = FALSE)))
  expect_false(isTRUE(m$loglik_recomputed))

  # bootstraps never warn
  expect_no_warning(suppressMessages(fit_evzinb_fast(n_bootstraps = 4)))
})

test_that("multicore = NULL leaves the user's future plan untouched (R0.4)", {
  oplan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(oplan), add = TRUE)
  plan_before <- future::plan()

  # multicore = NULL: evinf must not touch the plan the user set.
  evinf:::evinf_with_plan(NULL, NULL, invisible(NULL))
  expect_identical(class(future::plan()), class(plan_before))

  # multicore = TRUE / FALSE set a temporary plan but restore it on exit.
  f <- function(mc) evinf:::evinf_with_plan(mc, 2, future::plan())
  inner <- f(FALSE)
  expect_identical(class(future::plan()), class(plan_before))
  expect_true(inherits(inner, "sequential"))
})

test_that("predict() type='all' and confint use canonical state names first (R0.5)", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  nm_all <- names(predict(m, type = "all", quantile = 0.9))
  expect_true(all(c("pr_zero", "pr_count", "pr_evi") %in% nm_all))
  expect_lt(match("pr_zero", nm_all), match("pr_zc", nm_all))

  ci_zi <- suppressWarnings(predict(m, type = "zi", confint = TRUE))
  expect_equal(names(ci_zi)[1], "pr_zero")
  ci_ev <- suppressWarnings(predict(m, type = "evinf", confint = TRUE))
  expect_equal(names(ci_ev)[1], "pr_evi")
})

test_that("design matrix rejects non-finite columns (R0.8)", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$z <- c(0, d$x1[-1])          # a zero so log(z) = -Inf, row not otherwise NA
  expect_error(
    suppressWarnings(suppressMessages(
      evzinb(y ~ log(z), data = d, bootstrap = FALSE, verbose = FALSE)
    )),
    "[Nn]on-finite"
  )
})

test_that("fitted-value computation is unchanged by vectorisation (4.13)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  # values captured from the pre-vectorisation loop
  expect_equal(unname(head(m$fitted$y.hat.pl_median, 3)),
               c(42.1708, 64.5629, 324.8574), tolerance = 1e-3)
  expect_equal(m$log.lik, -251.2682021, tolerance = 1e-6)
})

test_that("print(summary()) reports observations at or above C_EV (4.2)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  expect_output(print(summary(m)), "Observations at or above C_EV")
})
