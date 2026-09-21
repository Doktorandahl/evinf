# round8 A.1 (audit §5.4): log_lik_profile_fun() profiles the whole C_EV
# candidate grid in one C++ call instead of em_profile_c()'s old vapply() over
# log_lik_fun(). This file checks the new function reproduces the old one
# exactly, and that em_profile_c()'s selected c_hat is unaffected on the
# identity fixtures (A.5).

profile_ext_par <- function(m) {
  # evinb() strips Beta.multinom.ZC / x.multinom.zc from the public object
  # (R/evinb.R): the zero-inflation block is fixed, internally, at a design
  # matrix equal to x.nb and coefficients c(-100, 0, ..., 0), which is the
  # convention em_fit() itself uses. Reconstruct it so log_lik_fun() /
  # log_lik_profile_fun() (which always take a zc block) see the same thing
  # the fit did.
  is_evinb <- is.null(m$coef$Beta.multinom.ZC)
  x_zc <- if (is_evinb) m$data$x.nb else m$data$x.multinom.zc
  par <- m$coef
  if (is_evinb) par$Beta.multinom.ZC <- c(-100, rep(0, ncol(m$data$x.nb)))

  x_obj <- list(
    X.multinom.ZC = x_zc,
    X.multinom.PL = m$data$x.multinom.pl,
    X.NB          = m$data$x.nb,
    X.PL          = m$data$x.pl,
    offset.nb     = m$offset_nb
  )
  list(
    ext = evinf:::em_extend_design(x_obj, length(m$data$y)),
    y   = m$data$y,
    par = par
  )
}

# vapply(candidates, log_lik_fun, ...) vs. one log_lik_profile_fun() call, at
# the fitted model's own parameters and its actual c_profile grid.
expect_profile_matches_loop <- function(m, tol = 1e-10) {
  pe <- profile_ext_par(m)
  cands <- m$c_profile$c

  old <- vapply(cands, function(cc) {
    evinf:::log_lik_fun(
      pe$par$Beta.multinom.ZC, pe$par$Beta.multinom.PL, pe$par$Beta.NB,
      pe$par$Alpha.NB, pe$par$Beta.PL, cc,
      pe$ext$zc, pe$ext$pl_mult, pe$ext$nb, pe$ext$pl, pe$y, pe$ext$offset,
      pe$ext$offset_zc, pe$ext$offset_pl_mult, pe$ext$weights
    )
  }, numeric(1))

  new <- as.numeric(evinf:::log_lik_profile_fun(
    pe$par$Beta.multinom.ZC, pe$par$Beta.multinom.PL, pe$par$Beta.NB,
    pe$par$Alpha.NB, pe$par$Beta.PL, cands,
    pe$ext$zc, pe$ext$pl_mult, pe$ext$nb, pe$ext$pl, pe$y, pe$ext$offset,
      pe$ext$offset_zc, pe$ext$offset_pl_mult, pe$ext$weights
  ))

  testthat::expect_equal(new, old, tolerance = tol)
}

test_that("log_lik_profile_fun() matches vapply(candidates, log_lik_fun, ...) on genevzinb2", {
  data(genevzinb2, package = "evinf", envir = environment())
  ctrl <- evinf::evinf_control(c.lim = c(50, 1000), init.C = 200)
  m <- suppressMessages(suppressWarnings(
    evinf::evzinb(y ~ x1 + x2 + x3, data = genevzinb2, control = ctrl,
                 verbose = FALSE, bootstrap = FALSE)
  ))
  expect_profile_matches_loop(m)
})

test_that("log_lik_profile_fun() matches vapply(candidates, log_lik_fun, ...) on an evinb fit", {
  data(genevzinb2, package = "evinf", envir = environment())
  ctrl <- evinf::evinf_control(c.lim = c(50, 1000), init.C = 200)
  m <- suppressMessages(suppressWarnings(
    evinf::evinb(y ~ x1 + x2 + x3, data = genevzinb2, control = ctrl,
                verbose = FALSE, bootstrap = FALSE)
  ))
  expect_profile_matches_loop(m)
})

test_that("log_lik_profile_fun() matches vapply(candidates, log_lik_fun, ...) on an offset() model", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$loff <- log(2)
  ctrl <- evinf::evinf_control(c.lim = c(50, 1000), init.C = 200)
  m <- suppressMessages(suppressWarnings(
    evinf::evzinb(y ~ x1 + x2 + x3 + offset(loff), data = d, control = ctrl,
                 verbose = FALSE, bootstrap = FALSE)
  ))
  expect_profile_matches_loop(m)
})

test_that("log_lik_profile_fun() matches vapply(candidates, log_lik_fun, ...) on hks", {
  skip_on_cran()
  data(hks, package = "evinf", envir = environment())
  f_hks <- osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
    lntpop + brv_AllLag_log + osvAllLagDum + incomp
  f_hks_pareto <- ~ troopLag_log + epdur_log + brv_AllLag_log + lntpop
  m <- suppressMessages(suppressWarnings(
    evinf::evzinb(f_hks, formula_pareto = f_hks_pareto, data = hks,
                 bootstrap = FALSE, verbose = FALSE)
  ))
  expect_profile_matches_loop(m, tol = 1e-8)
})

test_that("log_lik_profile_fun() handles a single-element candidate vector", {
  data(genevzinb2, package = "evinf", envir = environment())
  ctrl <- evinf::evinf_control(c.lim = c(50, 1000), init.C = 200)
  m <- suppressMessages(suppressWarnings(
    evinf::evzinb(y ~ x1 + x2 + x3, data = genevzinb2, control = ctrl,
                 verbose = FALSE, bootstrap = FALSE)
  ))
  pe <- profile_ext_par(m)
  cc <- m$c_profile$c[1]

  old <- evinf:::log_lik_fun(
    pe$par$Beta.multinom.ZC, pe$par$Beta.multinom.PL, pe$par$Beta.NB,
    pe$par$Alpha.NB, pe$par$Beta.PL, cc,
    pe$ext$zc, pe$ext$pl_mult, pe$ext$nb, pe$ext$pl, pe$y, pe$ext$offset,
      pe$ext$offset_zc, pe$ext$offset_pl_mult, pe$ext$weights
  )
  new <- as.numeric(evinf:::log_lik_profile_fun(
    pe$par$Beta.multinom.ZC, pe$par$Beta.multinom.PL, pe$par$Beta.NB,
    pe$par$Alpha.NB, pe$par$Beta.PL, cc,
    pe$ext$zc, pe$ext$pl_mult, pe$ext$nb, pe$ext$pl, pe$y, pe$ext$offset,
      pe$ext$offset_zc, pe$ext$offset_pl_mult, pe$ext$weights
  ))
  expect_length(new, 1L)
  expect_equal(new, old, tolerance = 1e-10)
})

test_that("em_profile_c()'s selected c_hat is unchanged on all 6 identity fixtures (round8 A.5)", {
  # Not just close -- the fixtures' c_ev came from a discrete candidate grid,
  # so log_lik_profile_fun() must pick exactly the same candidate as the old
  # vapply(candidates, log_lik_fun, ...), not merely a numerically close one.
  ctrl_id <- function() evinf::evinf_control(c.lim = c(50, 1000), init.C = 200)
  fit_z_id <- function(formula, data, ...) {
    suppressMessages(suppressWarnings(
      evinf::evzinb(formula, data = data, control = ctrl_id(), verbose = FALSE, ...)
    ))
  }
  fit_i_id <- function(formula, data, ...) {
    suppressMessages(suppressWarnings(
      evinf::evinb(formula, data = data, control = ctrl_id(), verbose = FALSE, ...)
    ))
  }

  data(genevzinb2, package = "evinf", envir = environment())
  gf <- local({
    set.seed(42)
    d <- genevzinb2
    d$g <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
    d
  })
  f3 <- y ~ x1 + log(abs(x2) + 1) + g

  check_c_hat <- function(m, baseline_file) {
    base <- readRDS(test_path("fixtures", baseline_file))
    expect_identical(m$coef$C, base$c_ev, label = baseline_file)
  }

  check_c_hat(fit_z_id(y ~ x1 + x2 + x3, genevzinb2, bootstrap = FALSE),
              "em_baseline_evzinb.rds")
  check_c_hat(fit_i_id(y ~ x1 + x2 + x3, genevzinb2, bootstrap = FALSE),
              "em_baseline_evinb.rds")
  check_c_hat(fit_z_id(f3, gf, bootstrap = FALSE), "em_baseline_evzinb_factor.rds")
  check_c_hat(fit_i_id(f3, gf, bootstrap = FALSE), "em_baseline_evinb_factor.rds")

  skip_on_cran()
  data(hks, package = "evinf", envir = environment())
  f_hks <- osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
    lntpop + brv_AllLag_log + osvAllLagDum + incomp
  f_hks_pareto <- ~ troopLag_log + epdur_log + brv_AllLag_log + lntpop
  check_c_hat(fit_z_id(f_hks, hks, formula_pareto = f_hks_pareto, bootstrap = FALSE),
              "em_baseline_evzinb_hks.rds")
  check_c_hat(fit_i_id(f_hks, hks, formula_pareto = f_hks_pareto, bootstrap = FALSE),
              "em_baseline_evinb_hks.rds")
})
