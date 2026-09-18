# Round 4 Q1: em_fit(model = "evinb") holds the zero-inflation block fixed in
# BOTH phases (previously the warm-up ran the EVZINB EM step, which updates it).
#
# The change is numerically inert on any input evinb() can produce (the ZI
# intercept is hard-initialised at -100, so exp(-100 + x'b) ~ 0 and the spurious
# update never moved anything) - test-em-identity.R passes unchanged. These are
# characterisation tests that lock in the corrected internal behaviour.

make_xo <- function() {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- evinf:::evinf_design(y ~ x1 + x2 + x3, genevzinb2)
  list(
    y = genevzinb2$y,
    xo = list(X.multinom.ZC = d$X, X.multinom.PL = d$X, X.NB = d$X,
              X.PL = d$X, offset.nb = rep(0, nrow(d$X))),
    np = ncol(d$X) + 1L
  )
}

em_ini <- function(np, zc) {
  list(Beta.multinom.ZC = zc, Beta.multinom.PL = rep(0, np),
       Beta.NB = rep(0, np), Alpha.NB = 0.01, Beta.PL = rep(0, np), C = 200)
}

quiet_em_fit <- function(...) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(evinf:::em_fit(...)))
  out
}

test_that("em_fit(model = 'evinb') returns Beta.multinom.ZC exactly as supplied", {
  f <- make_xo()
  ctl <- evinf::evinf_control(c.lim = c(50, 1000), init.C = 200)
  zc <- c(-100, 0.3, -0.7, 1.1)
  r <- quiet_em_fit(f$y, f$xo, em_ini(f$np, zc), ctl, model = "evinb")
  expect_identical(r$par.mat$Beta.multinom.ZC, zc)
})

test_that("em_fit(model = 'evinb') is invariant to the ZI slopes", {
  f <- make_xo()
  ctl <- evinf::evinf_control(c.lim = c(50, 1000), init.C = 200)
  r0 <- quiet_em_fit(f$y, f$xo, em_ini(f$np, c(-100, 0, 0, 0)), ctl, model = "evinb")
  r1 <- quiet_em_fit(f$y, f$xo, em_ini(f$np, c(-100, 5, -5, 3)), ctl, model = "evinb")
  expect_equal(r0$par.mat$Beta.NB, r1$par.mat$Beta.NB, tolerance = 1e-10)
  expect_equal(r0$par.mat$Beta.PL, r1$par.mat$Beta.PL, tolerance = 1e-10)
  expect_identical(r0$par.mat$C, r1$par.mat$C)
})

test_that("em_fit(model = 'evzinb') still updates the ZI block", {
  f <- make_xo()
  ctl <- evinf::evinf_control(c.lim = c(50, 1000), init.C = 200)
  r <- quiet_em_fit(f$y, f$xo, em_ini(f$np, rep(0, f$np)), ctl, model = "evzinb")
  expect_gt(max(abs(r$par.mat$Beta.multinom.ZC)), 0)
})

# audit0.10 §1.4: the C_EV outer loop had no iteration cap.

test_that("em_fit() stops with converge = FALSE when the C_EV profile oscillates (audit0.10 §1.4)", {
  f <- make_xo()
  ctl <- evinf::evinf_control(c.lim = c(50, 1000), init.C = 200, max.c.iter = 5)

  flip_state <- new.env()
  flip_state$flip <- TRUE
  fake_profile <- function(y, x_obj, par, c_candidates) {
    c_hat <- if (flip_state$flip) 100 else 200
    flip_state$flip <- !flip_state$flip
    list(profile = data.frame(c = c_candidates, loglik = rep(-100, length(c_candidates))),
         c_hat = c_hat, loglik_max = -100)
  }
  testthat::local_mocked_bindings(em_profile_c = fake_profile, .package = "evinf")

  # Both phases oscillate here, so both the convergence-phase warning and the
  # round8 0.6 warm-up-phase warning fire; collect them explicitly instead of
  # letting expect_warning() only check the first and leak the second.
  collected <- character(0)
  r <- withCallingHandlers(
    quiet_em_fit(f$y, f$xo, em_ini(f$np, rep(0, f$np)), ctl, model = "evzinb"),
    warning = function(w) {
      collected <<- c(collected, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("^em_fit\\(\\): the C_EV profile did not settle", collected)))
  expect_true(any(grepl("warm-up C_EV profile did not settle", collected)))
  expect_false(r$converge)
  expect_false(r$c_converged)
  expect_true(r$c_warmup_capped)
})

test_that("em_fit() records a warm-up-only C_EV cap separately from convergence (round8 0.6, review §7)", {
  f <- make_xo()
  ctl <- evinf::evinf_control(c.lim = c(50, 1000), init.C = 200, max.c.iter = 3)

  call_state <- new.env()
  call_state$n <- 0L
  fake_profile <- function(y, x_obj, par, c_candidates) {
    call_state$n <- call_state$n + 1L
    # Oscillate for exactly the warm-up phase's calls (max.c.iter = 3), then
    # settle immediately once the convergence phase starts.
    c_hat <- if (call_state$n <= 3L) {
      if (call_state$n %% 2 == 1) 100 else 200
    } else {
      par$C
    }
    list(profile = data.frame(c = c_candidates, loglik = rep(-100, length(c_candidates))),
         c_hat = c_hat, loglik_max = -100)
  }
  testthat::local_mocked_bindings(em_profile_c = fake_profile, .package = "evinf")

  collected <- character(0)
  r <- withCallingHandlers(
    quiet_em_fit(f$y, f$xo, em_ini(f$np, rep(0, f$np)), ctl, model = "evzinb"),
    warning = function(w) {
      collected <<- c(collected, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(r$c_warmup_capped)
  expect_true(r$c_converged)
  expect_true(r$converge)
  expect_true(any(grepl("warm-up C_EV profile did not settle", collected)))
  expect_false(any(grepl("^em_fit\\(\\): the C_EV profile did not settle", collected)))
})

test_that("object$props / object$resp are computed at the final parameters, not one EM step behind (audit0.10 §1.11, D.3)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)

  # Recompute props directly at the model's own final coefficients (mirrors
  # em_fit()'s internal evinf_stable_props3() call at the returned parameters).
  ext_zc <- cbind(1, m$data$x.multinom.zc)
  ext_pl_mult <- cbind(1, m$data$x.multinom.pl)
  expected_props <- evinf:::evinf_stable_props3(
    as.numeric(ext_zc %*% m$coef$Beta.multinom.ZC),
    as.numeric(ext_pl_mult %*% m$coef$Beta.multinom.PL)
  )
  expect_equal(unname(as.matrix(m$props)), unname(expected_props), tolerance = 1e-8)

  # resp is the E-step responsibilities at that same props / final coef and
  # final C_EV, not one EM step behind them.
  expected_resp <- evinf:::evinf_responsibilities(
    m$data$y, m$fitted$mu.nb, m$coef$Alpha.NB, m$fitted$alpha.pl,
    m$coef$C, expected_props
  )
  expect_equal(unname(as.matrix(m$resp)), unname(expected_resp), tolerance = 1e-8)

  mi <- fit_evinb_fast(bootstrap = FALSE)
  ext_pl_mult_i <- cbind(1, mi$data$x.multinom.pl)
  expected_props_i <- evinf:::evinf_stable_props3(
    rep(-Inf, nrow(ext_pl_mult_i)),
    as.numeric(ext_pl_mult_i %*% mi$coef$Beta.multinom.PL)
  )
  expect_equal(unname(as.matrix(mi$props)), unname(expected_props_i[, c("count", "evi")]),
               tolerance = 1e-8)
})

test_that("em_profile_c() returns a single c_hat on an exact tie (audit0.10 §1.4)", {
  f <- make_xo()
  par <- em_ini(f$np, rep(0, f$np))
  testthat::local_mocked_bindings(log_lik_fun = function(...) -50, .package = "evinf")

  prof <- evinf:::em_profile_c(f$y, f$xo, par, c_candidates = c(100, 200, 300))
  expect_length(prof$c_hat, 1L)
  expect_identical(prof$c_hat, 100)
})
