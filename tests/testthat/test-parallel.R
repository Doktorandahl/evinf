# Round 5 Part A: the future/furrr backend.
#
# furrr::furrr_options(seed = <int>) hands each element (bootstrap / refit) its
# own L'Ecuyer stream, so the *resampling* is identical under plan(sequential)
# and plan(multisession) for the same boot_seed: every $boot_id matches
# bit-for-bit, and the full-sample fit (which is not parallelised) is identical.
# The per-bootstrap fits then agree to within cross-process floating-point noise
# in the iterative optimisers (Accelerate / OpenBLAS are not bit-reproducible
# across processes); the seeds below are chosen so no resample sits on a
# likelihood ridge where that noise would flip the EM into another basin, so a
# tight tolerance still holds. Degenerate resamples are the subject of Part C.
#
# The multisession halves are skip_on_cran() (spinning up workers is expensive);
# there is deliberately no skip_on_os() - CI runs them on every platform.

ctl <- .fast_control()

coefs_num <- function(m) {
  ce <- suppressWarnings(coefficient_extractor(m, "all"))
  as.matrix(ce[vapply(ce, is.numeric, logical(1))])
}
boot_ids_of <- function(m) lapply(m$bootstraps, function(b)
  if (inherits(b, "try-error")) NULL else b$boot_id)

# Run `code` under a plan, restoring the previous plan as soon as it returns.
seq_run <- function(code) {
  op <- future::plan("sequential"); on.exit(future::plan(op), add = TRUE)
  force(code)
}
ms_run <- function(code, workers = 2) {
  op <- future::plan(future::multisession, workers = workers)
  on.exit(future::plan(op), add = TRUE)
  force(code)
}

fit_z <- function(...) suppressWarnings(suppressMessages(evinf::evzinb(
  y ~ x1 + x2 + x3, data = genevzinb2, control = ctl, verbose = FALSE, ...)))
fit_i <- function(...) suppressWarnings(suppressMessages(evinf::evinb(
  y ~ x1 + x2 + x3, data = genevzinb2, control = ctl, verbose = FALSE, ...)))

data(genevzinb2, package = "evinf", envir = environment())

test_that("evzinb()/evinb() bootstraps: sequential and multisession resample identically", {
  seq_z <- seq_run(fit_z(bootstrap = TRUE, n_bootstraps = 4, boot_seed = 202, multicore = NULL))
  seq_i <- seq_run(fit_i(bootstrap = TRUE, n_bootstraps = 4, boot_seed = 202, multicore = NULL))

  skip_on_cran()
  ms_z <- ms_run(fit_z(bootstrap = TRUE, n_bootstraps = 4, boot_seed = 202, multicore = NULL))
  ms_i <- ms_run(fit_i(bootstrap = TRUE, n_bootstraps = 4, boot_seed = 202, multicore = NULL))

  # the documented Reproducibility contract: boot_id is fully determined by
  # boot_seed and independent of the backend AND the number of workers.
  expect_identical(boot_ids_of(ms_z), boot_ids_of(seq_z))
  expect_identical(boot_ids_of(ms_i), boot_ids_of(seq_i))
  ms_z3 <- ms_run(fit_z(bootstrap = TRUE, n_bootstraps = 4, boot_seed = 202,
                        multicore = NULL), workers = 3)
  expect_identical(boot_ids_of(ms_z3), boot_ids_of(seq_z))

  expect_equal(ms_z$coef, seq_z$coef, tolerance = 1e-12)   # full fit is not parallel
  expect_equal(ms_i$coef, seq_i$coef, tolerance = 1e-12)
  expect_equal(coefs_num(ms_z), coefs_num(seq_z), tolerance = 1e-6)
  expect_equal(coefs_num(ms_i), coefs_num(seq_i), tolerance = 1e-6)
})

test_that("add_bootstraps(): sequential and multisession resample identically", {
  base <- seq_run(fit_z(bootstrap = TRUE, n_bootstraps = 3, boot_seed = 202, multicore = NULL))
  seq_ab <- seq_run(suppressWarnings(add_bootstraps(base, 3, boot_seed = 303, multicore = NULL)))

  skip_on_cran()
  ms_ab <- ms_run(suppressWarnings(add_bootstraps(base, 3, boot_seed = 303, multicore = NULL)))
  expect_identical(boot_ids_of(ms_ab), boot_ids_of(seq_ab))
  expect_equal(coefs_num(ms_ab), coefs_num(seq_ab), tolerance = 1e-6)
})

test_that("lr_test(bootstrap = TRUE): sequential ~= multisession", {
  m <- seq_run(fit_z(bootstrap = TRUE, n_bootstraps = 4, boot_seed = 31, multicore = NULL))
  seq_lr <- seq_run(suppressWarnings(suppressMessages(
    lr_test(m, "x1", bootstrap = TRUE, multicore = NULL))))

  skip_on_cran()
  ms_lr <- ms_run(suppressWarnings(suppressMessages(
      lr_test(m, "x1", bootstrap = TRUE, multicore = NULL))))
  expect_equal(ms_lr$results$statistic, seq_lr$results$statistic, tolerance = 1e-4)
  expect_equal(ms_lr$boot_results$x1$statistic,
               seq_lr$boot_results$x1$statistic, tolerance = 1e-3)
})

test_that("compare_models(): sequential ~= multisession", {
  m <- seq_run(fit_z(bootstrap = TRUE, n_bootstraps = 5, boot_seed = 41, multicore = NULL))
  seq_cm <- seq_run(suppressWarnings(suppressMessages(
    compare_models(m, multicore = NULL))))
  seq_oob <- suppressWarnings(oob_evaluation(seq_cm, metric = "rmse")$nb)

  skip_on_cran()
  ms_cm <- ms_run(suppressWarnings(suppressMessages(compare_models(m, multicore = NULL))))
  ms_oob <- suppressWarnings(oob_evaluation(ms_cm, metric = "rmse")$nb)
  expect_equal(ms_oob, seq_oob, tolerance = 1e-6)
})

test_that("predict(pred = 'bootstrap_median', confint = TRUE): sequential ~= multisession", {
  m <- seq_run(fit_z(bootstrap = TRUE, n_bootstraps = 5, boot_seed = 51, multicore = NULL))
  seq_p <- seq_run(suppressWarnings(
    predict(m, type = "harmonic", pred = "bootstrap_median",
            confint = TRUE, multicore = NULL)))

  skip_on_cran()
  ms_p <- ms_run(suppressWarnings(predict(m, type = "harmonic", pred = "bootstrap_median",
                             confint = TRUE, multicore = NULL)))
  expect_equal(as.data.frame(ms_p), as.data.frame(seq_p), tolerance = 1e-4)
})

test_that("marginal_effects(): sequential ~= multisession", {
  m <- seq_run(fit_z(bootstrap = TRUE, n_bootstraps = 5, boot_seed = 61, multicore = NULL))
  seq_me <- seq_run(suppressWarnings(
    marginal_effects(m, variables = "x1", type = "states", multicore = NULL)))

  skip_on_cran()
  ms_me <- ms_run(suppressWarnings(
      marginal_effects(m, variables = "x1", type = "states", multicore = NULL)))
  expect_equal(as.data.frame(ms_me), as.data.frame(seq_me), tolerance = 1e-4)
})

test_that("multicore = TRUE/FALSE match the corresponding explicit plans, user plan restored", {
  skip_on_cran()
  user_plan <- future::plan(future::multisession, workers = 2)
  withr::defer(future::plan(user_plan))

  # multicore = FALSE == plan(sequential)
  mc_false <- fit_z(bootstrap = TRUE, n_bootstraps = 3, boot_seed = 71, multicore = FALSE)
  seq_expl <- seq_run(fit_z(bootstrap = TRUE, n_bootstraps = 3, boot_seed = 71, multicore = NULL))
  expect_identical(boot_ids_of(mc_false), boot_ids_of(seq_expl))
  expect_equal(coefs_num(mc_false), coefs_num(seq_expl), tolerance = 1e-3)

  # multicore = TRUE == plan(multisession)
  mc_true <- fit_z(bootstrap = TRUE, n_bootstraps = 3, boot_seed = 72,
                   multicore = TRUE, ncores = 2)
  ms_expl <- ms_run(fit_z(bootstrap = TRUE, n_bootstraps = 3, boot_seed = 72, multicore = NULL))
  expect_identical(boot_ids_of(mc_true), boot_ids_of(ms_expl))
  expect_equal(coefs_num(mc_true), coefs_num(ms_expl), tolerance = 1e-3)

  # the user's plan is back
  expect_true(inherits(future::plan(), "multisession"))
})
