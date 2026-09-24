# round10 Part G (audit §5.10): multiple starts, trace plots, per-bootstrap
# convergence detail, and check_evinf().

# --- G.1: n_starts ------------------------------------------------------

# round11 B2: multi-start (n_starts > 1) refits several times over; kept fast
# on CI (NOT_CRAN=true) but skipped on CRAN's own check-time budget.
testthat::skip_on_cran()

test_that("n_starts = 1 is bit-identical to a plain fit (round10 G.1)", {
  data(genevzinb2, package = "evinf", envir = environment())
  m_plain <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = FALSE, verbose = FALSE
  )))
  m_explicit <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(n_starts = 1)
  )))
  expect_identical(m_plain$log.lik, m_explicit$log.lik)
  expect_identical(unname(m_plain$coef$Beta.NB), unname(m_explicit$coef$Beta.NB))
  expect_null(m_plain$starts)
  expect_null(m_explicit$starts)
  expect_null(m_plain$start_seed)
})

test_that("n_starts = 10 recovers at least the known-good optimum on the factor example (round10 G.1)", {
  d <- genevzinb2_factor()
  m10 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + g, data = d, bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(n_starts = 10, start_jitter = 0.5)
  )))
  # round8/round9 review history: the good optimum on this example is
  # log-lik ~= -253.555; a worse accumulation path can land near -253.68 or
  # below. n_starts = 10 should reach at least the good one (often better).
  expect_gt(m10$log.lik, -253.6)
  expect_s3_class(m10$starts, "tbl_df")
  expect_equal(nrow(m10$starts), 10L)
  expect_setequal(names(m10$starts),
                  c("start", "loglik", "C", "converged", "n_em_steps",
                    "loglik_trace", "c_trace"))
  expect_true(is.numeric(m10$start_seed))
  # the best row's own loglik must equal the fit's reported log-lik
  expect_equal(max(m10$starts$loglik), m10$log.lik, tolerance = 1e-8)
})

test_that("n_starts reproduces with the same start_seed", {
  d <- genevzinb2_factor()
  fit <- function(seed) suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + g, data = d, bootstrap = FALSE, verbose = FALSE, multicore = FALSE,
    start_seed = seed, control = evinf_control(n_starts = 4, start_jitter = 0.5)
  )))
  m_a <- fit(123)
  m_b <- fit(123)
  expect_identical(m_a$starts$loglik, m_b$starts$loglik)
  expect_identical(m_a$log.lik, m_b$log.lik)
})

test_that("glance()/print() report n_starts / n_starts_at_best (round10 G.1)", {
  d <- genevzinb2_factor()
  m10 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + g, data = d, bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(n_starts = 5, start_jitter = 0.5)
  )))
  g <- glance(m10)
  expect_equal(g$n_starts, 5L)
  expect_true(g$n_starts_at_best >= 1L && g$n_starts_at_best <= 5L)
  expect_output(print(m10), "starts reached the best log-likelihood")

  m1 <- fit_evzinb_fast(bootstrap = FALSE)
  g1 <- glance(m1)
  expect_true(is.na(g1$n_starts))
  expect_true(is.na(g1$n_starts_at_best))
})

test_that("bootstraps warm-start from the winning (not the default) start (round10 G.1/D7)", {
  # The bootstrap warm-start already reads object$coef directly (unchanged
  # code path); with n_starts > 1 that's the *winning* start's coefficients
  # by construction, since object$coef is only ever set from starts_result$best.
  d <- genevzinb2_factor()
  m10 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + g, data = d, n_bootstraps = 2, boot_seed = 1, multicore = FALSE,
    verbose = FALSE, control = evinf_control(n_starts = 5, start_jitter = 0.5)
  )))
  expect_length(m10$bootstraps, 2L)
  ok <- !vapply(m10$bootstraps, inherits, logical(1), "try-error")
  expect_true(any(ok))
})

# --- G.2: plot(type = "trace") -------------------------------------------

test_that("plot(type = 'trace') works for single- and multi-start fits (round10 G.2)", {
  rlang::check_installed("ggplot2")
  m1 <- fit_evzinb_fast(bootstrap = FALSE)
  p1 <- plot(m1, type = "trace")
  expect_s3_class(p1, "ggplot")

  d <- genevzinb2_factor()
  m5 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + g, data = d, bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(n_starts = 5, start_jitter = 0.5)
  )))
  p5 <- plot(m5, type = "trace")
  expect_s3_class(p5, "ggplot")
  # multi-start data includes every start plus the "best" overlay
  expect_true(nrow(p5$data) > nrow(p1$data))
})

# --- G.3: per-bootstrap convergence detail --------------------------------

test_that("n_em_steps/c_converged/c_warmup_capped are stored on every bootstrap replicate (round10 G.3)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  for (b in m$bootstraps) {
    if (inherits(b, "try-error")) next
    expect_true(is.numeric(b$n_em_steps) && b$n_em_steps > 0)
    expect_true(is.logical(b$c_converged))
    expect_true(is.logical(b$c_warmup_capped))
  }

  mi <- fit_evinb_fast(n_bootstraps = 5)
  for (b in mi$bootstraps) {
    if (inherits(b, "try-error")) next
    expect_true(is.numeric(b$n_em_steps) && b$n_em_steps > 0)
  }
})

test_that("failed_bootstraps() surfaces not_converged replicates (round10 G.3)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  # Inject a replicate that "ran" but didn't converge, without being an
  # error or flagged degenerate.
  b1 <- m$bootstraps[[1]]
  b1$converge <- FALSE
  b1$degenerate <- FALSE
  b1$n_em_steps <- 3L
  b1$c_converged <- FALSE
  b1$c_warmup_capped <- TRUE
  m$bootstraps[[1]] <- b1

  fb <- failed_bootstraps(m)
  expect_true("not_converged" %in% fb$type)
  row <- fb[fb$type == "not_converged", ]
  expect_equal(unname(row$n_em_steps), 3L)
  expect_false(row$c_converged)
  expect_true(row$c_warmup_capped)
})

test_that("glance() reports median_boot_em_steps / n_boot_c_capped (round10 G.3)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  g <- glance(m)
  expect_true(is.numeric(g$median_boot_em_steps) && g$median_boot_em_steps > 0)
  expect_true(is.numeric(g$n_boot_c_capped) && g$n_boot_c_capped >= 0)

  m0 <- fit_evzinb_fast(bootstrap = FALSE)
  g0 <- glance(m0)
  expect_true(is.na(g0$median_boot_em_steps))
  expect_true(is.na(g0$n_boot_c_capped))
})

# --- G.4: check_evinf() ----------------------------------------------------

test_that("check_evinf() returns a classed tibble with sensible statuses (round10 G.4)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  cc <- check_evinf(m)
  expect_s3_class(cc, "evinf_check")
  expect_true(all(c("check", "status", "detail") %in% names(cc)))
  expect_true(all(cc$status %in% c("ok", "note", "warning")))
  expect_true(all(c(
    "EM converged", "C_EV profile converged (convergence phase)",
    "C_EV profile converged (warm-up phase)", "C_EV on candidate-grid boundary",
    "Fitted Pareto shape (alpha_pl) not collapsed", "Bootstrap replicates usable"
  ) %in% cc$check))
  expect_output(print(cc), "evinf diagnostic report")
})

test_that("check_evinf() adds a start-agreement row only when n_starts > 1 (round10 G.4)", {
  m1 <- fit_evzinb_fast(bootstrap = FALSE)
  expect_false("Multiple starts agree" %in% check_evinf(m1)$check)

  d <- genevzinb2_factor()
  m5 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + g, data = d, bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(n_starts = 5, start_jitter = 0.5)
  )))
  expect_true("Multiple starts agree" %in% check_evinf(m5)$check)
})

test_that("check_evinf() adds an OOB-fraction row only for the block schemes (round10 G.4)", {
  m_iid <- fit_evzinb_fast(n_bootstraps = 3)
  expect_false("Out-of-bag fraction (block scheme)" %in% check_evinf(m_iid)$check)

  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$id <- rep(1:10, each = 10)
  d$tm <- rep(1:10, times = 10)
  m_block <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2, data = d, block = id, time = tm, bootstrap_scheme = "moving_block",
    n_bootstraps = 3, boot_seed = 1, multicore = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200)
  )))
  expect_true("Out-of-bag fraction (block scheme)" %in% check_evinf(m_block)$check)
})

test_that("check_evinf() rejects a non-evinf object", {
  expect_error(check_evinf(list()), "fitted evzinb / evinb")
})

# round11 C2: the tests above only ever see a healthy fit (every status "ok"
# or "note"); trip each "warning"-capable check directly by mutating the
# fields check_evinf() reads, so the branch that actually reports "warning"
# is exercised at least once per check.
test_that("check_evinf() reports 'warning' for EM non-convergence", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  m$converge <- FALSE
  cc <- check_evinf(m)
  expect_equal(cc$status[cc$check == "EM converged"], "warning")
})

test_that("check_evinf() reports 'warning' when the C_EV profile did not converge", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  m$c_converged <- FALSE
  cc <- check_evinf(m)
  expect_equal(cc$status[cc$check == "C_EV profile converged (convergence phase)"], "warning")
})

test_that("check_evinf() reports 'warning' when C_EV sits on a candidate-grid boundary", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  m$coef$C <- min(m$c_profile$c)
  cc <- check_evinf(m)
  expect_equal(cc$status[cc$check == "C_EV on candidate-grid boundary"], "warning")
})

test_that("check_evinf() reports 'warning' when the fitted Pareto shape has collapsed", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  m$fitted$alpha.pl[1] <- 1e-20
  cc <- check_evinf(m)
  expect_equal(cc$status[cc$check == "Fitted Pareto shape (alpha_pl) not collapsed"], "warning")
})

test_that("check_evinf() reports 'note' (not 'ok') when a bootstrap replicate is unusable", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  m$bootstraps[[1]] <- try(stop("boom"), silent = TRUE)
  cc <- check_evinf(m)
  expect_equal(cc$status[cc$check == "Bootstrap replicates usable"], "note")
})
