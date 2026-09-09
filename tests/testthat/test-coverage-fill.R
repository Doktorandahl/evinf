# issue 4.3 - targeted coverage for the thin internal files
# (parallel.R, extra_distributions.R) and the
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

test_that("evinf_pmap() runs its progressr path and is seed-stable", {
  old <- progressr::handlers("void")
  on.exit(progressr::handlers(old), add = TRUE)

  out <- evinf:::evinf_pmap(1:3, function(i) i^2, seed = 1L, label = "x",
                            verbose = TRUE)
  expect_identical(out, list(1, 4, 9))

  # per-element L'Ecuyer streams -> identical regardless of chunk_size
  a <- evinf:::evinf_pmap(1:5, function(i) stats::runif(1), seed = 42L)
  b <- evinf:::evinf_pmap(1:5, function(i) stats::runif(1), seed = 42L,
                          chunk_size = 5L)
  expect_equal(unlist(a), unlist(b))

  expect_error(evinf:::evinf_pmap(1:3, function(i) i, seed = c(1, 2)),
               "single integer")
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

# furrr::furrr_options(seed = <int>) hands each bootstrap its own L'Ecuyer
# stream, so the resample indices ($boot_id) are bit-identical under
# multicore = TRUE and FALSE. The fitted coefficients then agree to within
# cross-process floating-point noise in the EM optimiser (which a handful of
# near-degenerate resamples can amplify), so they are compared at a loose
# tolerance. (The broader plan/multicore matrix lives in test-parallel.R.)
coef_num <- function(m) {
  ce <- suppressWarnings(coefficient_extractor(m, "all"))
  as.matrix(ce[vapply(ce, is.numeric, logical(1))])
}
boot_ids_of <- function(m) lapply(m$bootstraps, function(b)
  if (inherits(b, "try-error")) NULL else b$boot_id)

test_that("multicore = TRUE bootstrap == multicore = FALSE for the same boot_seed", {
  skip_on_cran()

  data(genevzinb2, package = "evinf", envir = environment())
  ctl <- .fast_control()
  fitz <- function(mc) suppressWarnings(suppressMessages(evinf::evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, control = ctl, bootstrap = TRUE,
    n_bootstraps = 4, boot_seed = 202, multicore = mc, ncores = 2,
    verbose = FALSE)))

  mz_seq <- fitz(FALSE); mz_par <- fitz(TRUE)
  expect_identical(boot_ids_of(mz_par), boot_ids_of(mz_seq))
  expect_equal(coef_num(mz_par), coef_num(mz_seq), tolerance = 1e-6)

  # add_bootstraps(): the new batch must match across backends
  ab_seq <- suppressWarnings(add_bootstraps(mz_seq, 3, boot_seed = 303,
                                            multicore = FALSE))
  ab_par <- suppressWarnings(add_bootstraps(mz_seq, 3, boot_seed = 303,
                                            multicore = TRUE, ncores = 2))
  expect_identical(boot_ids_of(ab_par), boot_ids_of(ab_seq))
  expect_equal(coef_num(ab_par), coef_num(ab_seq), tolerance = 1e-6)

  # compare_models(): the out-of-bag metric aggregates over all bootstraps,
  # so cross-process noise averages down.
  cm_seq <- suppressWarnings(suppressMessages(
    compare_models(mz_seq, multicore = FALSE)))
  cm_par <- suppressWarnings(suppressMessages(
    compare_models(mz_seq, multicore = TRUE, ncores = 2)))
  expect_equal(
    suppressWarnings(oob_evaluation(cm_par, metric = "rmse")$nb),
    suppressWarnings(oob_evaluation(cm_seq, metric = "rmse")$nb),
    tolerance = 1e-6
  )
})
