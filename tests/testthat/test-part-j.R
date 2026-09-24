# round10 Part J (audit §5.11): engineering -- runtime estimate (J.1),
# chunking (J.2), evinf_bench_data() (J.3).

# --- J.3: evinf_bench_data() -------------------------------------------------

# round11 B2: benchmark-data generation and chunking/runtime-estimate tests;
# kept fast on CI (NOT_CRAN=true) but skipped on CRAN's own check-time budget.
testthat::skip_on_cran()

test_that("evinf_bench_data() is reproducible, correctly shaped, and leaves the caller's RNG state untouched", {
  set.seed(123)
  pre <- .Random.seed
  d1 <- evinf_bench_data(200, seed = 1)
  expect_identical(.Random.seed, pre)

  d2 <- evinf_bench_data(200, seed = 1)
  expect_identical(d1, d2)

  expect_equal(nrow(d1), 200L)
  expect_equal(names(d1), c("y", "x1", "x2", "x3"))
  expect_true(all(d1$y >= 0 & d1$y == round(d1$y)))
  expect_true(mean(d1$y == 0) > 0 && mean(d1$y == 0) < 1)  # a genuine mixture, not degenerate

  d_noseed <- evinf_bench_data(50)
  expect_equal(nrow(d_noseed), 50L)
})

test_that("evinf_bench_data() produces data evzinb() fits successfully, at small and large n", {
  d <- evinf_bench_data(300, seed = 7)
  m <- suppressMessages(suppressWarnings(
    evzinb(y ~ x1 + x2 + x3, data = d, bootstrap = FALSE, verbose = FALSE)
  ))
  expect_true(m$converge)

  d_big <- evinf_bench_data(5000, seed = 8)
  m_big <- suppressMessages(suppressWarnings(
    evzinb(y ~ x1 + x2 + x3, data = d_big, bootstrap = FALSE, verbose = FALSE)
  ))
  expect_true(m_big$converge)
})

# --- J.1: bootstrap runtime projection --------------------------------------

test_that("evinf_estimate_boot_runtime() emits only when verbose = TRUE or the projection exceeds 60s", {
  fake_bootrun <- function(spec, blk, tv) {
    Sys.sleep(0.001)
    list(ok = TRUE)
  }

  # n_bootstraps <= 4: the probe would cover every replicate, so there is
  # nothing left to project -- no message regardless of verbose.
  expect_no_message(
    evinf:::evinf_estimate_boot_runtime(fake_bootrun, list(), NULL, NULL, 4, 1L, verbose = TRUE)
  )
  expect_no_message(
    evinf:::evinf_estimate_boot_runtime(fake_bootrun, list(), NULL, NULL, 2, 1L, verbose = FALSE)
  )

  # verbose = TRUE: always messages, even for a fast (well under 60s) projection.
  expect_message(
    evinf:::evinf_estimate_boot_runtime(fake_bootrun, list(), NULL, NULL, 8, 1L, verbose = TRUE),
    "projected bootstrap runtime"
  )

  # verbose = FALSE, fast probe -> projection well under 60s: silent.
  expect_no_message(
    evinf:::evinf_estimate_boot_runtime(fake_bootrun, list(), NULL, NULL, 8, 1L, verbose = FALSE)
  )

  # verbose = FALSE, but a large n_bootstraps projects > 60s from the same
  # per-replicate time: messages anyway.
  expect_message(
    evinf:::evinf_estimate_boot_runtime(fake_bootrun, list(), NULL, NULL, 200000, 1L, verbose = FALSE),
    "projected bootstrap runtime"
  )
})

test_that("evinf_estimate_boot_runtime()'s probe does not change the real dispatch's bootstrap output", {
  # The probe (R/parallel.R's own comment) deliberately recomputes positions
  # 1..k under the *same* boot_seed as the real dispatch, precisely so the
  # real dispatch below stays untouched -- this is an end-to-end regression
  # guard for that: two full evzinb() bootstrap runs with the same seed give
  # identical boot_id sequences (the actual resampling, which is what the
  # real dispatch's reproducibility promise is about).
  data(genevzinb2, package = "evinf", envir = environment())
  ctrl <- evinf_control(c.lim = c(50, 1000), init.C = 200)
  m1 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = TRUE, n_bootstraps = 8,
    multicore = FALSE, boot_seed = 999, control = ctrl, verbose = FALSE
  )))
  m2 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = TRUE, n_bootstraps = 8,
    multicore = FALSE, boot_seed = 999, control = ctrl, verbose = FALSE
  )))
  same <- mapply(function(a, b) isTRUE(all.equal(a$boot_id, b$boot_id)),
                 m1$bootstraps, m2$bootstraps)
  expect_true(all(same))
})

# --- J.2: chunking -----------------------------------------------------------

test_that("evinf_pmap()'s default chunk_size doesn't change results, and neither does evinf_control(chunk_size = )", {
  future::plan("sequential")
  r1 <- evinf:::evinf_pmap(1:20, function(i) runif(2), seed = 7, chunk_size = 1L)
  r5 <- evinf:::evinf_pmap(1:20, function(i) runif(2), seed = 7, chunk_size = 5L)
  rAuto <- evinf:::evinf_pmap(1:20, function(i) runif(2), seed = 7)
  expect_identical(r1, r5)
  expect_identical(r1, rAuto)
})

test_that("evinf_control(chunk_size = )'s default is NULL (\"auto\"), and it validates an explicit override", {
  expect_null(evinf_control()$chunk_size)
  expect_identical(evinf_control(chunk_size = 3L)$chunk_size, 3L)
  expect_error(evinf_control(chunk_size = -1), "positive integer")
  expect_error(evinf_control(chunk_size = 1.5), "positive integer")
})

test_that("evzinb()'s bootstrap replicates are identical across chunk_size = 1, 3 and the auto default (round10 J.2)", {
  data(genevzinb2, package = "evinf", envir = environment())
  run_with <- function(cs) {
    suppressMessages(suppressWarnings(evzinb(
      y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = TRUE, n_bootstraps = 12,
      multicore = FALSE, boot_seed = 555,
      control = evinf_control(c.lim = c(50, 1000), init.C = 200, chunk_size = cs),
      verbose = FALSE
    )))
  }
  m1 <- run_with(1L)
  m3 <- run_with(3L)
  m_auto <- run_with(NULL)
  ids1 <- lapply(m1$bootstraps, `[[`, "boot_id")
  ids3 <- lapply(m3$bootstraps, `[[`, "boot_id")
  ids_auto <- lapply(m_auto$bootstraps, `[[`, "boot_id")
  expect_identical(ids1, ids3)
  expect_identical(ids1, ids_auto)
})

test_that("evzinb()/evinb() bootstrap dispatch: no message for a short run, and a message under verbose = TRUE for n_bootstraps > 4", {
  data(genevzinb2, package = "evinf", envir = environment())
  ctrl <- evinf_control(c.lim = c(50, 1000), init.C = 200)

  expect_no_message(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = TRUE, n_bootstraps = 8,
    multicore = FALSE, boot_seed = 1, control = ctrl, verbose = FALSE
  )))

  expect_message(
    suppressWarnings(evinb(
      y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = TRUE, n_bootstraps = 8,
      multicore = FALSE, boot_seed = 1, control = ctrl, verbose = TRUE
    )),
    "projected bootstrap runtime"
  )
})
