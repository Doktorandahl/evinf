# audit0.10 §1.11: the M-step likelihood is now computed in log space
# throughout (Pareto pmf, NB pmf, mixture combination, softmax state
# probabilities) so it no longer underflows to -Inf for a heavy-tailed /
# large-count fit.

test_that("the stable and naive Pareto log-pmf agree on a moderate grid", {
  naive_ell_pl <- function(a, C, y) log((C / y)^a - (C / (y + 1))^a)

  grid <- expand.grid(a = c(0.05, 0.1, 0.5, 1, 2, 5, 10),
                       y = c(1, 2, 5, 20, 100, 1000, 1e5))
  grid$C <- 1
  naive <- with(grid, naive_ell_pl(a, C, y))
  stable <- vapply(seq_len(nrow(grid)), function(i)
    evinf:::ell_pl_i_fun(log(grid$a[i]), grid$C[i], 1, grid$y[i]), numeric(1))

  expect_equal(stable, naive, tolerance = 1e-10)
})

test_that("the stable Pareto log-pmf stays finite where the naive one underflows to -Inf", {
  # The audit's own y = 1e7 / alpha = 1e-3 / C = 100 example turns out not to
  # underflow on IEEE-754 double arithmetic (the cancellation is real but
  # stays well above 1e-16 relative to the terms). This combination is chosen
  # to actually force naive cancellation to exactly 0 (verified below).
  C <- 100; y <- 1e15; a <- 1e-5
  naive <- log((C / y)^a - (C / (y + 1))^a)
  expect_false(is.finite(naive))

  stable <- evinf:::ell_pl_i_fun(log(a), C, 1, y)
  expect_true(is.finite(stable))
})

test_that("the closed-form and looped NB log-pmf agree", {
  naive_ell_nb <- function(beta, alpha, x, y) {
    mu <- exp(sum(x * beta))
    ell <- 0
    if (y > 0) {
      for (j in 0:(y - 1)) ell <- ell + log(j + 1 / alpha)
      for (j in 1:y) ell <- ell - log(j)
    }
    ell - (1 / alpha) * log(1 + alpha * mu) - y * log(1 + alpha * mu) +
      y * log(alpha) + y * log(mu)
  }

  x <- c(1, 0.3); beta <- c(0.1, -0.2)
  for (alpha in c(0.05, 0.5, 2)) {
    for (y in c(0L, 1L, 5L, 50L, 500L)) {
      expect_equal(
        evinf:::ell_nb_i_fun(beta, alpha, x, y),
        naive_ell_nb(beta, alpha, x, y),
        tolerance = 1e-10,
        info = sprintf("alpha=%s y=%s", alpha, y)
      )
    }
  }
})

test_that("the closed-form NB log-pmf stays accurate at large y (round8 0.9)", {
  # review §7 / round8 0.9: for a large y with a small 1/alpha, the closed
  # form's three lgamma() terms partially cancel. Checked against a
  # Kahan-summed version of the naive loop (plain summation of ~1e6 terms is
  # itself not a trustworthy reference at this scale -- it drifts by ~1e-8 on
  # its own); switched the count-likelihood closed form from
  # lgamma(y+r) - lgamma(r) - lgamma(y+1) to the algebraically identical
  # -log(y+r) - lbeta(y+1, r), which was at least as accurate in every case
  # checked (max relative error ~1.5e-10 for lgamma, ~1e-10 for lbeta at
  # y = 1e6, alpha = 5).
  skip_on_cran()

  kahan_sum <- function(x) {
    s <- 0; c <- 0
    for (v in x) {
      y <- v - c; t <- s + y; c <- (t - s) - y; s <- t
    }
    s
  }
  naive_count_term_kahan <- function(alpha, y) {
    r <- 1 / alpha
    if (y == 0) return(0)
    kahan_sum(log(seq(0, y - 1) + r)) - kahan_sum(log(seq_len(y)))
  }

  x <- c(1, 0.3); beta <- c(0.1, -0.2)
  for (alpha in c(5, 50)) {
    for (y in c(1e5, 1e6)) {
      mu <- exp(sum(x * beta))
      count_term_ref <- naive_count_term_kahan(alpha, y)
      full_ref <- count_term_ref - (1 / alpha) * log(1 + alpha * mu) -
        y * log(1 + alpha * mu) + y * log(alpha) + y * log(mu)
      expect_equal(
        evinf:::ell_nb_i_fun(beta, alpha, x, y),
        full_ref,
        tolerance = 1e-9,
        info = sprintf("alpha=%s y=%s", alpha, y)
      )
    }
  }
})

test_that("a large-count model fits without -Inf (audit0.10 §1.11)", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  # Scale the extreme-value tail up by several orders of magnitude, so the
  # Pareto pmf is evaluated at counts large enough to matter numerically.
  d$y[d$y >= 200] <- d$y[d$y >= 200] * 1e6

  m <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + x2 + x3, data = d,
    control = evinf_control(c.lim = c(50 * 1e6, 1000 * 1e6)),
    bootstrap = FALSE, verbose = FALSE
  )))
  expect_true(is.finite(m$log.lik))
})

test_that("the closed-form and looped NB score/Hessian (alpha component) agree (round8 A.2)", {
  # audit §5.4 / round8 A.2: delldtheta_nb_i_fun()'s and
  # d2elldtheta2_nb_i_fun()'s O(y) loops for the dispersion derivatives are
  # replaced by digamma/trigamma closed forms (for y > 30 in the Hessian --
  # see src/evinf.cpp for why the closed form cancels catastrophically below
  # that, at a small alpha). Grid matches the prompt: y in
  # {0, 1, 5, 100, 1e4, 1e5} x alpha in {1e-3, 0.01, 0.5, 1, 5, 50}.
  old_grad_term <- function(y, alpha) if (y == 0) 0 else sum(1 / (seq(0, y - 1) + 1 / alpha))
  old_hess_term <- function(y, alpha) {
    if (y == 0) return(0)
    sum((seq(0, y - 1) / (1 + alpha * seq(0, y - 1)))^2)
  }

  x <- c(1, 0.3, -0.2); beta <- c(0.1, -0.2, 0.05)
  mu <- exp(sum(x * beta))
  n <- length(x)

  for (alpha in c(1e-3, 0.01, 0.5, 1, 5, 50)) {
    for (y in c(0, 1, 5, 100, 1e4, 1e5)) {
      g <- evinf:::delldtheta_nb_i_fun(beta, alpha, x, y)
      h <- evinf:::d2elldtheta2_nb_i_fun(beta, alpha, x, y)

      expected_g <- if (y == 0) {
        log(1 + alpha * mu) / alpha^2 - mu / (alpha * (1 + alpha * mu))
      } else {
        (log(1 + alpha * mu) - old_grad_term(y, alpha)) / alpha^2 +
          (y - mu) / (alpha * (1 + alpha * mu))
      }
      expected_h <- -old_hess_term(y, alpha) -
        2 / alpha^3 * log(1 + alpha * mu) + (2 / alpha^2) * mu / (1 + alpha * mu) +
        (y + 1 / alpha) * mu^2 / (1 + alpha * mu)^2

      info <- sprintf("y=%s alpha=%s", y, alpha)
      expect_equal(g[length(g)], expected_g, tolerance = 1e-9, info = info)
      expect_equal(h[n + 1, n + 1], expected_h, tolerance = 1e-9, info = info)
    }
  }
})
