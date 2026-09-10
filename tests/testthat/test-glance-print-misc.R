test_that("glance() is unrounded and carries convergence / bootstrap counts (audit 2.10)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  m0 <- fit_evzinb_fast(bootstrap = FALSE)
  mi <- fit_evinb_fast(bootstrap = FALSE)

  g <- glance(m)
  expect_true(all(c("converged", "n_bootstraps", "n_failed_bootstraps",
                    "n_degenerate_bootstraps") %in% names(g)))
  expect_true(g$converged)
  expect_equal(g$n_bootstraps + g$n_failed_bootstraps +
                 g$n_degenerate_bootstraps, 5L)

  expect_true(is.na(glance(m0)$n_bootstraps))

  # evinb glance was rounded to 2 dp before
  expect_gt(abs(glance(mi)$logLik - round(glance(mi)$logLik, 2)), 0)
})

test_that("gm_evzinb has rows for the new glance columns (audit 2.10)", {
  data(gm_evzinb, package = "evinf", envir = environment())
  expect_true(all(c("converged", "n_bootstraps", "n_failed_bootstraps",
                    "n_degenerate_bootstraps", "n_c_on_boundary") %in%
                    gm_evzinb$raw))
})

test_that("summary has a print method (audit 2.6)", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  expect_output(print(summary(m)), "EVZINB model summary")
  expect_output(print(m), "C_EV")
  expect_output(print(m), "Observations at or above C_EV")
})

test_that("evinb() bootstraps are named (audit 2.1)", {
  mi <- fit_evinb_fast(n_bootstraps = 5)
  expect_equal(names(mi$bootstraps), paste0("bootstrap_", 1:5))
})

test_that("evinb par.all / AIC / BIC exclude the disabled zero-inflation parameters", {
  mi <- fit_evinb_fast(bootstrap = FALSE)        # y ~ x1 + x2 + x3
  # evi (4) + nb (4) + alpha (1) + pareto (4) + C (1) = 14; the 4 fixed ZC
  # coefficients are not counted.
  expect_equal(length(mi$par.all), 14L)
  expect_equal(glance(mi)$npar, 14L)
  expect_equal(mi$AIC, 2 * 14 - 2 * mi$log.lik)
  expect_equal(mi$BIC, log(length(mi$data$y)) * 14 - 2 * mi$log.lik)
})

test_that("pdf.pl.type is validated (audit 2.12)", {
  expect_error(evinf_control(pdf.pl.type = "bogus"))
  d <- { data(genevzinb2, package = "evinf", envir = environment()); genevzinb2 }
  expect_error(
    suppressWarnings(evinb(y ~ x1, data = d, bootstrap = FALSE, pdf.pl.type = "bogus"))
  )
})

test_that("revzinb_fit() returns one draw per row (audit 2.15)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  draw <- revzinb_fit(m)
  expect_length(draw, nrow(m$data$data))
  expect_true(all(draw >= 0))
})

test_that("fitting emits no 'sequential' backend warning with multicore = FALSE (audit 2.11)", {
  d <- { data(genevzinb2, package = "evinf", envir = environment()); genevzinb2 }
  warns <- character(0)
  withCallingHandlers(
    suppressMessages(evzinb(y ~ x1 + x2 + x3, data = d, n_bootstraps = 3,
                            multicore = FALSE, boot_seed = 1, verbose = FALSE)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("sequential|backend registered", warns)))
})
