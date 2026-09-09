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

test_that("the multicore backend path works end to end", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("doParallel")

  data(genevzinb2, package = "evinf", envir = environment())
  m <- suppressWarnings(suppressMessages(evinf::evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, control = .fast_control(),
    bootstrap = TRUE, n_bootstraps = 4, multicore = TRUE, ncores = 2,
    boot_seed = 1, verbose = FALSE
  )))
  expect_length(m$bootstraps, 4L)

  p <- suppressWarnings(predict(m, type = "harmonic", confint = TRUE,
                                multicore = TRUE, ncores = 2))
  expect_true(all(c("ci_lb", "ci_ub") %in% names(p)))

  comp <- suppressWarnings(suppressMessages(
    compare_models(m, multicore = TRUE, ncores = 2)
  ))
  expect_length(comp$nb$bootstraps, 4L)

  # a user-registered backend is left untouched when multicore = FALSE
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
