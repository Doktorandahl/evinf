# issue 4.3 - targeted coverage for the thin internal files
# (parallel_backend.R, progress.R, extra_distributions.R) and the
# bootstrap_mean / bootstrap_median prediction paths.

test_that("nbinomdist2() validates its arguments", {
  expect_error(evinf:::nbinomdist2(size = 10), "prob or mu")
  expect_error(evinf:::nbinomdist2(prob = "a"), "numeric")
  expect_error(evinf:::nbinomdist2(prob = 1.5), "\\[0,1\\]")
  expect_error(evinf:::nbinomdist2(size = 2.5, prob = 0.3), "integer")
  expect_error(evinf:::nbinomdist2(mu = "a"), "numeric")
  d <- evinf:::nbinomdist2(mu = 3, size = 2)
  expect_s3_class(d, "nbinomdist")
})

test_that("evinf_progress_run() runs its progressr path", {
  old <- progressr::handlers("void")
  on.exit(progressr::handlers(old), add = TRUE)

  # verbose = TRUE forces the progressr::with_progress() branch
  expect_true(evinf:::evinf_progress_active(TRUE))
  out <- evinf:::evinf_progress_run(3, TRUE, function(p) {
    for (i in 1:3) p()
    "done"
  })
  expect_identical(out, "done")
})

test_that("bootstrap_mean / bootstrap_median predictions run for every type", {
  m  <- fit_evzinb_fast(n_bootstraps = 6)
  mi <- fit_evinb_fast(n_bootstraps = 6)
  n <- nrow(m$data$data)

  for (pr in c("bootstrap_mean", "bootstrap_median")) {
    for (ty in c("harmonic", "explog", "counts", "pareto_alpha")) {
      pz <- suppressWarnings(predict(m, type = ty, pred = pr))
      pv <- suppressWarnings(predict(mi, type = ty, pred = pr))
      expect_length(as.numeric(unlist(pz)), n)
      expect_length(as.numeric(unlist(pv)), n)
    }
    expect_length(
      as.numeric(unlist(suppressWarnings(
        predict(m, type = "quantile", quantile = 0.9, pred = pr)))),
      n
    )
  }
})

# Q2: exercise the multicore (PSOCK on Windows CI, fork elsewhere) bootstrap
# path and assert it is identical to the sequential path for the same boot_seed
# -- the %dorng% stream is backend-independent, so results must match to 1e-10.
# skip_on_cran() only (was skip_on_os("windows")) so R-CMD-check.yaml runs it.
coef_num <- function(m) {
  ce <- suppressWarnings(coefficient_extractor(m, "all"))
  as.matrix(ce[vapply(ce, is.numeric, logical(1))])
}

test_that("multicore bootstrap == sequential for the same boot_seed (Q2)", {
  skip_on_cran()
  skip_if_not_installed("doParallel")

  data(genevzinb2, package = "evinf", envir = environment())
  ctl <- .fast_control()
  fitz <- function(mc) suppressWarnings(suppressMessages(evinf::evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, control = ctl, bootstrap = TRUE,
    n_bootstraps = 4, boot_seed = 99, multicore = mc, ncores = 2,
    verbose = FALSE)))
  fiti <- function(mc) suppressWarnings(suppressMessages(evinf::evinb(
    y ~ x1 + x2 + x3, data = genevzinb2, control = ctl, bootstrap = TRUE,
    n_bootstraps = 4, boot_seed = 99, multicore = mc, ncores = 2,
    verbose = FALSE)))

  mz_seq <- fitz(FALSE); mz_par <- fitz(TRUE)
  expect_equal(coef_num(mz_par), coef_num(mz_seq), tolerance = 1e-10)

  mi_seq <- fiti(FALSE); mi_par <- fiti(TRUE)
  expect_equal(coef_num(mi_par), coef_num(mi_seq), tolerance = 1e-10)

  # add_bootstraps(): the new batch must match across backends
  ab_seq <- suppressWarnings(add_bootstraps(mz_seq, 3, boot_seed = 7,
                                            multicore = FALSE))
  ab_par <- suppressWarnings(add_bootstraps(mz_seq, 3, boot_seed = 7,
                                            multicore = TRUE, ncores = 2))
  expect_equal(coef_num(ab_par), coef_num(ab_seq), tolerance = 1e-10)

  # predict(pred = "bootstrap_median")
  p_seq <- suppressWarnings(predict(mz_seq, type = "harmonic",
                                    pred = "bootstrap_median", multicore = FALSE))
  p_par <- suppressWarnings(predict(mz_seq, type = "harmonic",
                                    pred = "bootstrap_median",
                                    multicore = TRUE, ncores = 2))
  expect_equal(as.numeric(p_par), as.numeric(p_seq), tolerance = 1e-10)

  # lr_test(bootstrap = TRUE)
  lr_seq <- suppressWarnings(suppressMessages(
    lr_test(mz_seq, "x1", bootstrap = TRUE, multicore = FALSE)))
  lr_par <- suppressWarnings(suppressMessages(
    lr_test(mz_seq, "x1", bootstrap = TRUE, multicore = TRUE, ncores = 2)))
  expect_equal(lr_par$results$chisq_mean, lr_seq$results$chisq_mean,
               tolerance = 1e-10)

  # compare_models()
  cm_seq <- suppressWarnings(suppressMessages(
    compare_models(mz_seq, multicore = FALSE)))
  cm_par <- suppressWarnings(suppressMessages(
    compare_models(mz_seq, multicore = TRUE, ncores = 2)))
  expect_equal(
    suppressWarnings(oob_evaluation(cm_par, metric = "rmse")$nb),
    suppressWarnings(oob_evaluation(cm_seq, metric = "rmse")$nb),
    tolerance = 1e-10
  )
})

test_that("a user-registered parallel backend survives multicore = FALSE", {
  skip_on_cran()
  skip_if_not_installed("doParallel")

  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
    foreach::registerDoSEQ()
  }, add = TRUE)
  be <- evinf:::evinf_setup_backend(multicore = FALSE)
  expect_true(be$user_backend)
  expect_identical(be$stop(), invisible(NULL))   # no-op: user's backend survives
  expect_true(foreach::getDoParRegistered())
})
