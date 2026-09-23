# round10 H.1 (audit §5.8, plan_extensions_5.1_5.2_5.4.md §B.1): evinf_pmf()/
# evinf_cdf() are the package's one definition of the predictive
# distribution, family-aware across all four (count x zero) combinations.
# mixture_p()/mixture_quantile()/residuals(type = "quantile")/
# classify_states() all build on the same shared pieces (R/evinf_pmf.R).

fit_family <- function(count, zero, ...) {
  data(genevzinb2, package = "evinf", envir = environment())
  suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2, data = genevzinb2, bootstrap = FALSE, verbose = FALSE,
    family = evinf_family(count = count, zero = zero),
    control = evinf_control(c.lim = c(50, 1000), init.C = 200), ...
  )))
}

test_that("evinf_pmf() rows sum to ~1 over a wide support, for every family", {
  data(genevzinb2, package = "evinf", envir = environment())
  wide <- 0:5000
  for (count in c("nbinom", "poisson")) {
    for (zero in c("mixture", "hurdle")) {
      m <- fit_family(count, zero)
      pmf <- evinf:::evinf_pmf(m, support = wide)
      expect_equal(dim(pmf), c(nrow(genevzinb2), length(wide)))
      rs <- rowSums(pmf)
      expect_true(all(rs > 0.98 & rs <= 1 + 1e-8),
                 info = paste(count, zero, "rowsum range", min(rs), max(rs)))
    }
  }
  mi <- suppressMessages(suppressWarnings(evinb(
    y ~ x1 + x2, data = genevzinb2, bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200)
  )))
  pmf <- evinf:::evinf_pmf(mi, support = wide)
  rs <- rowSums(pmf)
  expect_true(all(rs > 0.98 & rs <= 1 + 1e-8))
})

test_that("evinf_cdf() agrees with cumsum(evinf_pmf()), for every family", {
  wide <- 0:3000
  for (count in c("nbinom", "poisson")) {
    for (zero in c("mixture", "hurdle")) {
      m <- fit_family(count, zero)
      pmf <- evinf:::evinf_pmf(m, support = wide)
      cdf <- evinf:::evinf_cdf(m, support = wide)
      expect_equal(cdf, t(apply(pmf, 1, cumsum)), tolerance = 1e-8,
                  info = paste(count, zero))
    }
  }
})

test_that("evinf_pmf(y = ) matches evinf_pmf(support = )'s corresponding column", {
  data(genevzinb2, package = "evinf", envir = environment())
  m <- fit_family("nbinom", "mixture")
  y <- genevzinb2$y
  support <- 0:max(y)
  pmf_support <- evinf:::evinf_pmf(m, support = support)
  pmf_y <- evinf:::evinf_pmf(m, y = y)
  expected <- vapply(seq_along(y), function(i) pmf_support[i, y[i] + 1L], numeric(1))
  expect_equal(pmf_y, expected, tolerance = 1e-10)
})

test_that("evinf_pmf()/evinf_cdf() require exactly one of y / support", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_error(evinf:::evinf_pmf(m), "exactly one of")
  expect_error(evinf:::evinf_pmf(m, y = 1:5, support = 0:10), "exactly one of")
  expect_error(evinf:::evinf_cdf(m), "exactly one of")
})

test_that("classify_states(newdata = ) posterior matches evinf_responsibilities() directly (round10 H.1 refactor)", {
  data(genevzinb2, package = "evinf", envir = environment())
  m <- fit_evzinb_fast(bootstrap = FALSE)
  nd <- genevzinb2[1:10, ]
  cs <- classify_states(m, newdata = nd)

  dp <- evinf:::evinf_dist_params(m, newdata = nd)
  resp <- evinf:::evinf_responsibilities(nd$y, dp$nb_mu, dp$nb_alpha, dp$pl_alpha,
                                         dp$C, dp$probabilities, family = dp$family)
  expect_equal(cs$posterior_zero, unname(resp[, "zero"]), tolerance = 1e-10)
  expect_equal(cs$posterior_count, unname(resp[, "count"]), tolerance = 1e-10)
  expect_equal(cs$posterior_evi, unname(resp[, "evi"]), tolerance = 1e-10)
})

test_that("residuals(type = 'quantile') is unchanged by the evinf_cdf() refactor (round10 H.1)", {
  # Hand-computed the old way (direct mixture_p() calls on object$fitted/
  # object$props, unclamped) and compared against the refactored
  # residuals.evzinb(), on a fit where alpha_pl never collapses (so the new
  # clamp this refactor added is a no-op and the two must match exactly).
  m <- fit_evzinb_fast(bootstrap = FALSE)
  y <- m$data$y
  probs <- m$props
  mu <- m$fitted$mu.nb
  alph <- m$fitted$alpha.pl
  expect_true(min(alph) > (m$control$alpha_pl_floor %||% 0.01))  # clamp is a no-op here

  Fy_old <- mixture_p(y, alph, m$coef$C, mu, m$coef$Alpha.NB, probs,
                      family_count = "nbinom", family_zero = "mixture")
  Fy1_old <- ifelse(y <= 0, 0,
                    mixture_p(y - 1, alph, m$coef$C, mu, m$coef$Alpha.NB, probs,
                             family_count = "nbinom", family_zero = "mixture"))
  Fy_old <- pmin(pmax(Fy_old, 0), 1)
  Fy1_old <- pmin(pmax(Fy1_old, 0), Fy_old)

  Fy_new <- evinf:::evinf_cdf(m, y = y)
  Fy1_new <- ifelse(y <= 0, 0, evinf:::evinf_cdf(m, y = y - 1))

  expect_equal(Fy_new, Fy_old, tolerance = 1e-10)
  expect_equal(Fy1_new, Fy1_old, tolerance = 1e-10)

  set.seed(1)
  r1 <- residuals(m, type = "quantile", seed = 5)
  set.seed(1)
  r2 <- residuals(m, type = "quantile", seed = 5)
  expect_equal(r1, r2)
  expect_length(r1, length(y))
})
