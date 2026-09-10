# Round 5 Part C: degenerate bootstrap replicates (evinf_control(alpha_floor =)).

test_that("alpha_floor flags replicates and every summary drops them by default", {
  data(genevzinb2, package = "evinf", envir = environment())
  # an absurdly high floor forces (almost) every replicate to be degenerate,
  # regardless of the data (1e6 is a stress value, not realistic), so the plumbing is exercised deterministically.
  ctl <- evinf::evinf_control(c.lim = c(50, 1000), init.C = 200, alpha_floor = 1e6)
  m <- suppressWarnings(suppressMessages(evinf::evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, control = ctl,
    bootstrap = TRUE, n_bootstraps = 5, boot_seed = 123, multicore = FALSE,
    verbose = FALSE)))

  flags <- vapply(m$bootstraps, function(b) isTRUE(b$degenerate), logical(1))
  expect_true(all(flags))
  expect_true(all(vapply(m$bootstraps, function(b)
    grepl("alpha_floor", b$degenerate_reason), logical(1))))

  # failed_bootstraps(): type column, all "degenerate"
  fb <- failed_bootstraps(m)
  expect_named(fb, c("id", "type", "message"))
  expect_equal(nrow(fb), 5L)
  expect_true(all(fb$type == "degenerate"))

  # glance(): counted, and n_bootstraps is the usable count (0 here)
  g <- glance(m)
  expect_equal(g$n_degenerate_bootstraps, 5L)
  expect_equal(g$n_bootstraps, 0L)

  # bootstrap summaries refuse / empty by default, work with exclude_degenerate
  expect_error(suppressWarnings(vcov(m)), "usable bootstrap")
  expect_true(is.matrix(suppressWarnings(vcov(m, exclude_degenerate = FALSE))))
  ce_all  <- suppressWarnings(coefficient_extractor(m, "all"))
  ce_keep <- suppressWarnings(coefficient_extractor(m, "all", exclude_degenerate = FALSE))
  expect_equal(nrow(ce_all), 0L)
  expect_gt(nrow(ce_keep), 0L)
})

test_that("the default alpha_floor leaves well-behaved fits alone", {
  data(genevzinb2, package = "evinf", envir = environment())
  m <- suppressWarnings(suppressMessages(evinf::evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2,
    control = evinf::evinf_control(c.lim = c(50, 1000), init.C = 200),
    bootstrap = TRUE, n_bootstraps = 5, boot_seed = 123, multicore = FALSE,
    verbose = FALSE)))
  n_deg <- sum(vapply(m$bootstraps, function(b) isTRUE(b$degenerate), logical(1)))
  expect_equal(n_deg, 0L)
  expect_equal(glance(m)$n_degenerate_bootstraps, 0L)
  expect_equal(nrow(failed_bootstraps(m)), 0L)
})

test_that("a control from an older evinf (no alpha_floor) still validates", {
  ctl <- evinf::evinf_control()
  ctl$alpha_floor <- NULL
  ctl <- evinf:::validate_evinf_control(ctl)
  expect_equal(ctl$alpha_floor, 0.001)
})
