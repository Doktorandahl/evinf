# Round 5 Part C / round 6 Part A: degenerate bootstrap replicates
# (evinf_control(alpha_floor =, coef_limit =)).

# Build a real bootstrapped fit and mutate individual replicates in place, so
# each degeneracy criterion can be exercised without needing a pathological
# dataset.
fit_with_boots <- function(...) {
  data(genevzinb2, package = "evinf", envir = environment())
  suppressWarnings(suppressMessages(evinf::evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2,
    control = evinf::evinf_control(c.lim = c(50, 1000), init.C = 200),
    bootstrap = TRUE, n_bootstraps = 6, boot_seed = 123, multicore = FALSE,
    verbose = FALSE, ...)))
}

test_that("non-convergence and oversized coefficients flag a replicate (round 6)", {
  m <- fit_with_boots()
  reflag <- function(b) evinf:::evinf_flag_degenerate(
    b, m$data$x.pl[b$boot_id, , drop = FALSE], m$control)

  nc <- reflag(within(m$bootstraps[[1]], converge <- FALSE))
  expect_true(nc$degenerate)
  expect_match(nc$degenerate_reason, "did not converge")

  big <- m$bootstraps[[2]]
  big$coef$Beta.NB[["x1"]] <- 80
  big <- reflag(big)
  expect_true(big$degenerate)
  expect_match(big$degenerate_reason, "x1 = 80 in Beta.NB exceeds coef_limit")

  nab <- m$bootstraps[[3]]
  nab$coef$Alpha.NB <- Inf
  nab <- reflag(nab)
  expect_true(nab$degenerate)
  expect_match(nab$degenerate_reason, "non-finite coefficient in Alpha.NB")

  # convergence is checked before coefficients
  both <- m$bootstraps[[4]]
  both$converge <- FALSE
  both$coef$Beta.PL[[1]] <- 999
  expect_match(reflag(both)$degenerate_reason, "did not converge")
})

test_that("injected degenerate replicates flow through the accessors (round 6)", {
  m <- fit_with_boots()
  m$bootstraps[[1]]$converge <- FALSE
  m$bootstraps[[2]]$coef$Beta.NB[["x1"]] <- 80
  m$bootstraps <- lapply(m$bootstraps, function(b)
    evinf:::evinf_flag_degenerate(b, m$data$x.pl[b$boot_id, , drop = FALSE],
                                  m$control))

  fb <- failed_bootstraps(m)
  expect_setequal(fb$id[fb$type == "degenerate"], c("bootstrap_1", "bootstrap_2"))
  expect_match(fb$message[fb$id == "bootstrap_1"], "did not converge")

  expect_equal(glance(m)$n_degenerate_bootstraps, 2L)
  expect_equal(glance(m)$n_bootstraps, 4L)

  n_drop <- nrow(suppressWarnings(coefficient_extractor(m, "all")))
  n_keep <- nrow(suppressWarnings(
    coefficient_extractor(m, "all", exclude_degenerate = FALSE)))
  expect_equal(n_keep - n_drop, 2L * 4L)  # 2 replicates x 4 components

  # vcov() / confint() / tidy() all shift when the two replicates re-enter
  v_drop <- suppressWarnings(vcov(m))
  v_keep <- suppressWarnings(vcov(m, exclude_degenerate = FALSE))
  expect_false(isTRUE(all.equal(v_drop, v_keep)))
  td_keep <- suppressWarnings(tidy(m, exclude_degenerate = FALSE, standard_error = TRUE))
  td_drop <- suppressWarnings(tidy(m, standard_error = TRUE))
  expect_true(all(is.finite(td_drop$std.error)))
  expect_false(isTRUE(all.equal(td_keep$std.error, td_drop$std.error)))
})

test_that("C_EV on a candidate-grid endpoint is counted and warned, not excluded (round 6)", {
  data(genevzinb2, package = "evinf", envir = environment())
  # a candidate grid with only two observed values ({170, 173}, ~85th pctile):
  # every replicate's C_EV is then necessarily on an endpoint.
  ctl <- evinf::evinf_control(c.lim = c(170, 173), init.C = 170)
  expect_warning(
    m <- suppressMessages(evinf::evzinb(
      y ~ x1 + x2 + x3, data = genevzinb2, control = ctl,
      bootstrap = TRUE, n_bootstraps = 4, boot_seed = 1, multicore = FALSE,
      verbose = FALSE)),
    "boundary of the candidate range"
  )
  n_ok <- sum(!vapply(m$bootstraps, inherits, logical(1), "try-error"))
  expect_equal(m$n_c_on_boundary, n_ok)
  expect_equal(glance(m)$n_c_on_boundary, m$n_c_on_boundary)
  # informational, not degeneracy
  expect_equal(sum(vapply(m$bootstraps, function(b) isTRUE(b$degenerate),
                          logical(1))), 0L)
  expect_output(print(m), "boundary of the candidate range")
})


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

test_that("a control from an older evinf (no alpha_floor / coef_limit) still validates", {
  ctl <- evinf::evinf_control()
  ctl$alpha_floor <- NULL
  ctl$coef_limit <- NULL
  ctl <- evinf:::validate_evinf_control(ctl)
  expect_equal(ctl$alpha_floor, 0.001)
  expect_equal(ctl$coef_limit, 50)
})
