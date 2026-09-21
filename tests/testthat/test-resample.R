# round9 F (audit §5.7): evinf_resample_ids() and bootstrap_scheme = for
# panel/time-series bootstrapping.

test_that("iid scheme returns n indices in 1:n, ignoring block/time", {
  ids <- evinf:::evinf_resample_ids(20, "iid")
  expect_length(ids, 20L)
  expect_true(all(ids %in% 1:20))
})

test_that("cluster scheme draws whole units together and requires block_vec", {
  block <- rep(1:4, each = 5)
  ids <- evinf:::evinf_resample_ids(20, "cluster", block_vec = block)
  expect_length(ids, 20L)
  # every drawn unit contributes all 5 of its own rows, contiguously in the
  # output (evinf_resample_ids() concatenates whole `which(block == u)` runs).
  drawn_units <- split(ids, rep(seq_len(4), each = 5))
  for (u in drawn_units) {
    expect_length(unique(block[u]), 1L)
    expect_length(u, 5L)
  }
  expect_error(evinf:::evinf_resample_ids(20, "cluster"), "requires block_vec")
})

test_that("moving_block scheme respects unit boundaries and forms contiguous runs", {
  block <- rep(1:3, each = 10)
  time <- rep(1:10, times = 3)
  ids <- suppressMessages(evinf:::evinf_resample_ids(
    30, "moving_block", block_vec = block, time_vec = time, block_length = 3
  ))
  expect_length(ids, 30L)
  expect_true(all(ids %in% 1:30))
  # concatenated per unit, in unit order -- rows 1:10 come only from unit 1
  # (row indices == time here), etc.
  expect_true(all(ids[1:10] %in% 1:10))
  expect_true(all(ids[11:20] %in% 11:20))
  expect_true(all(ids[21:30] %in% 21:30))
  # blocks of length 3: within a unit's output, runs of consecutive integers
  # come in chunks of (up to) 3.
  d <- diff(ids[1:10])
  run_len <- rle(d == 1)$lengths[rle(d == 1)$values]
  expect_true(all(run_len <= 2))  # a run of k consecutive +1 diffs is a block of k+1
})

test_that("stationary scheme's mean block length matches L over many draws", {
  set.seed(1)
  n <- 500
  ids_lengths <- replicate(200, {
    ids <- suppressMessages(evinf:::evinf_resample_ids(
      n, "stationary", time_vec = seq_len(n), block_length = 10
    ))
    d <- diff(ids)
    breaks <- sum(d != 1 & d != -(n - 1))
    n / (breaks + 1)
  })
  expect_equal(mean(ids_lengths), 10, tolerance = 0.15)
})

test_that("moving_block/stationary require `time`", {
  expect_error(evinf:::evinf_resample_ids(20, "moving_block"), "requires `time`")
  expect_error(evinf:::evinf_resample_ids(20, "stationary"), "requires `time`")
})

test_that("a fixed seed reproduces indices exactly under every scheme", {
  block <- rep(1:3, each = 10)
  time <- rep(1:10, times = 3)
  for (scheme in c("iid", "cluster", "moving_block", "stationary")) {
    set.seed(42)
    a <- suppressMessages(evinf:::evinf_resample_ids(
      30, scheme, block_vec = if (scheme != "iid") block else NULL,
      time_vec = if (scheme %in% c("moving_block", "stationary")) time else NULL
    ))
    set.seed(42)
    b <- suppressMessages(evinf:::evinf_resample_ids(
      30, scheme, block_vec = if (scheme != "iid") block else NULL,
      time_vec = if (scheme %in% c("moving_block", "stationary")) time else NULL
    ))
    expect_identical(a, b, info = scheme)
  }
})

test_that("time must be strictly increasing within a unit -- duplicated or out-of-order errors, never silently sorted", {
  block <- rep(1:2, each = 5)
  time_dup <- c(1, 2, 3, 3, 5, 1:5)
  expect_error(
    evinf:::evinf_resample_ids(10, "moving_block", block_vec = block, time_vec = time_dup, block_length = 2),
    "not strictly increasing"
  )
  time_unsorted <- c(1, 2, 4, 3, 5, 1:5)
  expect_error(
    evinf:::evinf_resample_ids(10, "moving_block", block_vec = block, time_vec = time_unsorted, block_length = 2),
    "not strictly increasing"
  )
})

test_that("block_length = NULL defaults to ceiling(T^(1/3)) per unit, with a message", {
  block <- rep(1:2, each = c(8, 27)[1])  # placeholder, overwritten below
  block <- c(rep(1, 8), rep(2, 27))
  time <- c(1:8, 1:27)
  expect_message(
    ids <- evinf:::evinf_resample_ids(35, "moving_block", block_vec = block, time_vec = time),
    "ceiling\\(T\\^\\(1/3\\)\\)"
  )
  expect_length(ids, 35L)
})

test_that("a user-supplied block_length exceeding a unit's length errors rather than silently clamping", {
  block <- rep(1:2, each = 5)
  time <- rep(1:5, times = 2)
  expect_error(
    evinf:::evinf_resample_ids(10, "moving_block", block_vec = block, time_vec = time, block_length = 6),
    "exceeds the length"
  )
})

# --- Integration: bootstrap_scheme = via the public API --------------------

fit_panel_evzinb <- function(scheme, ...) {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  n_units <- 10
  T_per <- nrow(d) / n_units
  d$id <- rep(1:n_units, each = T_per)
  d$tm <- rep(1:T_per, times = n_units)
  suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2, data = d, block = id, time = tm, bootstrap_scheme = scheme,
    n_bootstraps = 4, boot_seed = 1, multicore = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200),
    ...
  )))
}

test_that("evzinb(bootstrap_scheme = 'moving_block'/'stationary') fits and bootstraps successfully", {
  for (scheme in c("moving_block", "stationary")) {
    m <- fit_panel_evzinb(scheme)
    expect_true(m$converge)
    expect_identical(m$bootstrap_scheme, scheme)
    expect_length(m$bootstraps, 4L)
    for (b in m$bootstraps) {
      if (inherits(b, "try-error")) next
      expect_length(b$boot_id, nrow(m$data$x.nb))
      expect_true(is.numeric(b$oob_fraction) && b$oob_fraction >= 0 && b$oob_fraction <= 1)
    }
  }
})

test_that("block = alone still defaults to the cluster scheme (unchanged pre-F behaviour)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_identical(m$bootstrap_scheme, "iid")
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$id <- rep(1:10, each = 10)
  m2 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2, data = d, block = id, bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200)
  )))
  expect_identical(m2$bootstrap_scheme, "cluster")
})

test_that("lr_test(bootstrap = TRUE), compare_models() and oob_evaluation() run under moving_block", {
  m <- fit_panel_evzinb("moving_block")
  lr <- lr_test(m, "x1", bootstrap = TRUE)
  expect_true(is.list(lr))
  cmp <- suppressWarnings(compare_models(m, multicore = FALSE))
  expect_true(all(c("nb", "zinb") %in% names(cmp)))
  oob <- oob_evaluation(cmp)
  expect_true("evinf" %in% names(oob))
})

test_that("add_bootstraps() carries the panel scheme forward without re-emitting the block_length message per replicate", {
  m <- fit_panel_evzinb("moving_block")
  msgs <- character(0)
  m2 <- withCallingHandlers(
    add_bootstraps(m, 2, boot_seed = 999),
    message = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleMessage") }
  )
  expect_length(m2$bootstraps, 6L)
  expect_false(any(grepl("block_length not supplied", msgs)))
})

test_that("update() reuses time / bootstrap_scheme / block_length", {
  m <- fit_panel_evzinb("stationary")
  m2 <- suppressMessages(suppressWarnings(update(m, formula_nb. = ~ .)))
  expect_identical(m2$bootstrap_scheme, "stationary")
  expect_identical(m2$time, "tm")
})

test_that("evinb() supports bootstrap_scheme = moving_block/stationary the same way", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$id <- rep(1:10, each = 10)
  d$tm <- rep(1:10, times = 10)
  mi <- suppressMessages(suppressWarnings(evinb(
    y ~ x1 + x2, data = d, block = id, time = tm, bootstrap_scheme = "moving_block",
    n_bootstraps = 3, boot_seed = 1, multicore = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200)
  )))
  expect_true(mi$converge)
  expect_identical(mi$bootstrap_scheme, "moving_block")
  for (b in mi$bootstraps) {
    if (inherits(b, "try-error")) next
    expect_true(is.numeric(b$oob_fraction))
  }
})
