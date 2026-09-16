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
