# round10 0.7 (review §6): the NB line search's eta interval used to let a
# trial step push alpha_new negative -- undefined for log_lik_fun()'s NB
# branch, so stats::optimise() silently substituted its own "maximum
# positive value" and warned on every such probe. On near-Poisson data
# (alpha_nb already close to 0), this produced a stream of warnings even
# though the *accepted* step (always inside the original eta.int) was fine.

test_that("evinf_eta_interval_positive_alpha() shrinks only the side that risks alpha <= 0", {
  # d > 0: alpha increases with eta, so only very negative eta is dangerous.
  iv <- evinf:::evinf_eta_interval_positive_alpha(c(-1, 1), alpha_old = 0.02, d = 0.5)
  eta_boundary <- -0.02 / 0.5
  expect_equal(iv[2], 1)  # upper end untouched
  expect_gt(iv[1], eta_boundary)  # shrunk strictly inside the boundary
  expect_lt(iv[1], 0)

  # d < 0: alpha decreases with eta, so only very positive eta is dangerous.
  iv2 <- evinf:::evinf_eta_interval_positive_alpha(c(-1, 1), alpha_old = 0.02, d = -0.5)
  eta_boundary2 <- -0.02 / -0.5
  expect_equal(iv2[1], -1)  # lower end untouched
  expect_lt(iv2[2], eta_boundary2)
  expect_gt(iv2[2], 0)

  # d == 0 (e.g. a Poisson count state): interval passes through unchanged.
  expect_identical(evinf:::evinf_eta_interval_positive_alpha(c(-1, 1), 0.02, 0), c(-1, 1))

  # Every eta the shrunk interval could return keeps alpha_new strictly positive.
  alpha_new <- function(eta, alpha_old, d) alpha_old + eta * d
  expect_gt(alpha_new(iv[1], 0.02, 0.5), 0)
  expect_gt(alpha_new(iv2[2], 0.02, -0.5), 0)
})

test_that("a near-Poisson nbinom fit produces no dispersion-boundary warnings (round10 0.7)", {
  # genevzinb2's own extreme-value rows keep the C_EV grid well-defined;
  # only the non-extreme rows are replaced with near-zero-dispersion Poisson
  # draws, pushing the fitted alpha_nb close to the Poisson limit.
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  extreme <- d$y >= 170
  set.seed(7)
  mu_low <- exp(0.3 + 0.1 * d$x1[!extreme])
  d$y[!extreme] <- rpois(sum(!extreme), mu_low)

  expect_no_warning(m <- suppressMessages(evzinb(
    y ~ x1 + x2 + x3, data = d, bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200)
  )))
  expect_true(m$coef$Alpha.NB > 0)
  expect_lt(m$coef$Alpha.NB, 0.05)  # confirms this data really is near-Poisson
})

test_that("default-family fits are unchanged by the eta-interval constraint (round10 0.7)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  # Pinned the same way as round10 0.5's tripwire; alpha_nb here is nowhere
  # near the boundary, so evinf_eta_interval_positive_alpha() is a no-op.
  data(genevzinb2, package = "evinf", envir = environment())
  m2 <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = FALSE, verbose = FALSE
  )))
  expect_equal(m2$log.lik, -254.03389449157694, tolerance = 1e-10)
})
