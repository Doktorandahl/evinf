# audit0.10 §1.8: pdf.pl.type = "exact" was accepted and stored but had no
# effect -- the M-step always used the continuous-Pareto derivatives.

test_that("the exact Pareto gradient/Hessian match central finite differences", {
  ell <- function(beta, C, x, y) evinf:::ell_pl_i_fun(beta, C, x, y)
  eps <- 1e-5

  grid <- list(
    list(beta = c(0.2, -0.3, 0.1), x = c(1, 0.5, -0.7), C = 50, y = 200),
    list(beta = c(-0.5, 0.1, 0.4), x = c(1, -1.2, 0.3), C = 100, y = 500),
    list(beta = c(0.05, 0.2, -0.1), x = c(1, 2, 1), C = 20, y = 21)
  )

  for (g in grid) {
    beta <- g$beta; x <- g$x; C <- g$C; y <- g$y
    ell_at <- function(b) ell(b, C, x, y)

    grad_fd <- vapply(seq_along(beta), function(j) {
      bp <- beta; bp[j] <- bp[j] + eps
      bm <- beta; bm[j] <- bm[j] - eps
      (ell_at(bp) - ell_at(bm)) / (2 * eps)
    }, numeric(1))
    grad_exact <- as.numeric(evinf:::delldbeta_pl_i_fun_exact(beta, C, x, y))
    expect_equal(grad_exact, grad_fd, tolerance = 1e-5)

    hess_fd <- matrix(0, length(beta), length(beta))
    for (j in seq_along(beta)) {
      bp <- beta; bp[j] <- bp[j] + eps
      bm <- beta; bm[j] <- bm[j] - eps
      gp <- vapply(seq_along(beta), function(k) {
        bpp <- bp; bpp[k] <- bpp[k] + eps
        bpm <- bp; bpm[k] <- bpm[k] - eps
        (ell_at(bpp) - ell_at(bpm)) / (2 * eps)
      }, numeric(1))
      gm <- vapply(seq_along(beta), function(k) {
        bmp <- bm; bmp[k] <- bmp[k] + eps
        bmm <- bm; bmm[k] <- bmm[k] - eps
        (ell_at(bmp) - ell_at(bmm)) / (2 * eps)
      }, numeric(1))
      hess_fd[, j] <- (gp - gm) / (2 * eps)
    }
    hess_exact <- evinf:::d2elldbeta2_pl_i_fun_exact(beta, C, x, y)
    expect_equal(as.matrix(hess_exact), hess_fd, tolerance = 1e-5)
  }
})

test_that("pdf.pl.type = 'exact' has no effect on the default fit (audit0.10 §1.8)", {
  m_approx <- fit_evzinb_fast(bootstrap = FALSE)
  m_exact <- fit_evzinb_fast(bootstrap = FALSE,
                             control = evinf::evinf_control(
                               c.lim = c(50, 1000), init.C = 200, pdf.pl.type = "exact"))
  # The default control uses "approx"; passing "exact" must change the fit
  # (this is the regression the audit found -- "exact" was a no-op).
  expect_false(isTRUE(all.equal(m_approx$coef$Beta.PL, m_exact$coef$Beta.PL)))
})

test_that("pdf.pl.type = 'exact' reaches a log-likelihood close to 'approx' on genevzinb2", {
  m_approx <- fit_evzinb_fast(bootstrap = FALSE)
  m_exact <- fit_evzinb_fast(bootstrap = FALSE,
                             control = evinf::evinf_control(
                               c.lim = c(50, 1000), init.C = 200, pdf.pl.type = "exact"))
  expect_gte(m_exact$log.lik, m_approx$log.lik - 1e-6)
})

test_that("pdf.pl.type = 'exact' reaches a log-likelihood close to 'approx' on hks", {
  skip_on_cran()
  data(hks, package = "evinf", envir = environment())
  f_hks <- osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
    lntpop + brv_AllLag_log + osvAllLagDum + incomp
  f_hks_pareto <- ~ troopLag_log + epdur_log + brv_AllLag_log + lntpop

  m_approx <- suppressWarnings(suppressMessages(evinf::evzinb(
    f_hks, formula_pareto = f_hks_pareto, data = hks, bootstrap = FALSE, verbose = FALSE
  )))
  m_exact <- suppressWarnings(suppressMessages(evinf::evzinb(
    f_hks, formula_pareto = f_hks_pareto, data = hks, bootstrap = FALSE, verbose = FALSE,
    control = evinf::evinf_control(pdf.pl.type = "exact")
  )))
  expect_gte(m_exact$log.lik, m_approx$log.lik - 1e-6)
})
