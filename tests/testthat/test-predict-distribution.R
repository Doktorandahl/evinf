# round10 H.2-H.7 (audit sec5.8, plan_extensions_5.1_5.2_5.4.md sec B.1):
# predict(type = "distribution"/"quantile" (vector)/"exceedance"/"draws"),
# all built on evinf_pmf()/evinf_cdf() (R/evinf_pmf.R) via the dedicated
# functions in R/predict_distribution.R.

fit_family <- function(count, zero, ...) {
  data(genevzinb2, package = "evinf", envir = environment())
  suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2, data = genevzinb2, bootstrap = TRUE, n_bootstraps = 4,
    boot_seed = 11, multicore = FALSE, verbose = FALSE,
    family = evinf_family(count = count, zero = zero),
    control = evinf_control(c.lim = c(50, 1000), init.C = 200), ...
  )))
}

test_that("predict(type = 'quantile', quantile = <vector>) is monotone and matches the scalar path column-by-column, for every family", {
  for (count in c("nbinom", "poisson")) {
    for (zero in c("mixture", "hurdle")) {
      m <- fit_family(count, zero)
      qs <- c(.5, .9, .99)
      qv <- suppressWarnings(predict(m, type = "quantile", quantile = qs))
      expect_s3_class(qv, "tbl_df")
      expect_equal(names(qv), c("q50", "q90", "q99"))
      expect_true(all(qv$q50 <= qv$q90 & qv$q90 <= qv$q99, na.rm = TRUE))
      for (j in seq_along(qs)) {
        q1 <- suppressWarnings(predict(m, type = "quantile", quantile = qs[j]))
        expect_equal(qv[[j]], q1, info = paste(count, zero, qs[j]))
      }
    }
  }
  mi <- suppressMessages(suppressWarnings(evinb(
    y ~ x1 + x2, data = genevzinb2, bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(c.lim = c(50, 1000), init.C = 200)
  )))
  qv <- suppressWarnings(predict(mi, type = "quantile", quantile = c(.5, .9)))
  expect_equal(names(qv), c("q50", "q90"))
})

test_that("a length-1 quantile is byte-identical to the pre-H.3 scalar path (round10 H.3 backward-compat requirement)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  q_new <- quantiles_from_evzinb(m, 0.9, newdata = NULL, round = TRUE)
  expect_true(is.numeric(q_new) && !is.matrix(q_new))
  q_cont <- quantiles_from_evzinb(m, 0.9, newdata = NULL, round = FALSE)
  expect_true(is.numeric(q_cont) && !is.matrix(q_cont))
  p1 <- predict(m, type = "quantile", quantile = 0.9)
  expect_equal(p1, q_new)
  expect_false(inherits(p1, "tbl_df"))
})

test_that("predict(type = 'quantile', quantile = <vector>, confint = TRUE) returns a well-formed interval", {
  # With only a handful of bootstrap replicates the empirical percentile
  # interval is itself noisy enough that the point estimate need not fall
  # inside it (same property the package's existing CI columns already
  # have -- test-predict.R only checks the columns exist, not bracketing) --
  # so this checks lo <= hi, not "encloses the point estimate".
  m <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 6)
  qv <- suppressWarnings(predict(m, type = "quantile", quantile = c(.5, .9), confint = TRUE))
  expect_equal(names(qv), c("q50", "q90", "q50_lo", "q50_hi", "q90_lo", "q90_hi"))
  expect_true(all(qv$q50_lo <= qv$q50_hi))
  expect_true(all(qv$q90_lo <= qv$q90_hi))
})

test_that("predict(type = 'distribution') rows sum to ~1 (format = 'matrix') and matches evinf_pmf() (format = 'long'), for every family", {
  support <- 0:5000
  for (count in c("nbinom", "poisson")) {
    for (zero in c("mixture", "hurdle")) {
      m <- fit_family(count, zero)
      dm <- predict(m, type = "distribution", support = support, format = "matrix")
      expect_equal(dim(dm), c(nrow(genevzinb2), length(support)))
      rs <- rowSums(dm)
      expect_true(all(rs > 0.98 & rs <= 1 + 1e-8))

      dl <- predict(m, type = "distribution", support = support, format = "long")
      expect_equal(nrow(dl), nrow(genevzinb2) * length(support))
      expect_equal(names(dl), c(".row", "y", "prob"))
      expect_equal(dl$prob[dl$.row == 1], unname(dm[1, ]), tolerance = 1e-10)
    }
  }
})

test_that("predict(type = 'distribution') default support is capped by max_support", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_message(
    d <- predict(m, type = "distribution", max_support = 50),
    "truncating"
  )
  expect_equal(max(d$y), 50)
})

test_that("predict(type = 'distribution', confint = TRUE) errors (round10 H.2: no confint support)", {
  m <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 3)
  expect_error(predict(m, type = "distribution", confint = TRUE), "not available")
})

test_that("exceedance(threshold = 0) is exactly 1, and exceedance matches 1 - cdf from distribution, for every family (round10 H.4/H.7)", {
  for (count in c("nbinom", "poisson")) {
    for (zero in c("mixture", "hurdle")) {
      m <- fit_family(count, zero)
      thresholds <- c(0, 1, 25, 1000)
      e <- suppressWarnings(predict(m, type = "exceedance", threshold = thresholds))
      expect_equal(names(e), paste0("p_ge_", thresholds))
      expect_true(all(e$p_ge_0 == 1))

      for (k in thresholds[thresholds > 0]) {
        Fkm1 <- evinf:::evinf_cdf(m, y = rep(k - 1, nrow(genevzinb2)))
        expect_equal(e[[paste0("p_ge_", k)]], 1 - Fkm1, tolerance = 1e-10,
                    info = paste(count, zero, "k =", k))
      }
    }
  }
})

test_that("predict(type = 'exceedance', confint = TRUE) returns a well-formed interval, and threshold = 0's interval is degenerate at 1", {
  m <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 6)
  e <- suppressWarnings(predict(m, type = "exceedance", threshold = c(0, 25), confint = TRUE))
  expect_true(all(e$p_ge_0_lo == 1 & e$p_ge_0_hi == 1))
  expect_true(all(e$p_ge_25_lo <= e$p_ge_25_hi))
})

test_that("predict(type = 'draws') returns non-negative integer-valued draws, reproducible when seeded, and does not perturb the caller's RNG state (round10 H.5)", {
  for (count in c("nbinom", "poisson")) {
    for (zero in c("mixture", "hurdle")) {
      m <- fit_family(count, zero)
      set.seed(999)
      pre <- .Random.seed
      d1 <- predict(m, type = "draws", n_draws = 20, seed = 42)
      expect_true(identical(.Random.seed, pre))
      expect_true(all(d1$y >= 0 & d1$y == round(d1$y)))
      expect_equal(nrow(d1), nrow(genevzinb2) * 20)

      d2 <- predict(m, type = "draws", n_draws = 20, seed = 42)
      expect_equal(d1, d2)
    }
  }
})

test_that("predict(type = 'draws', parameter_uncertainty = TRUE) draws from usable bootstrap replicates and is reproducibly seeded", {
  m <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 8)
  d1 <- predict(m, type = "draws", n_draws = 15, parameter_uncertainty = TRUE, seed = 7)
  d2 <- predict(m, type = "draws", n_draws = 15, parameter_uncertainty = TRUE, seed = 7)
  expect_equal(d1, d2)
  expect_true(all(d1$y >= 0))
})

test_that("seeded predictive draws match predict(type = 'distribution')'s pmf (binned chi-square goodness-of-fit, round10 H.7)", {
  # Not a mean comparison (review flagged the predictive mean only exists for
  # alpha_pl > 1, with finite MC error only for alpha_pl > 2) -- a binned
  # chi-square test against the full predictive distribution instead.
  m <- fit_evzinb_fast(bootstrap = FALSE)
  nd <- genevzinb2[1, ]
  n_draws <- 20000
  draws <- predict(m, type = "draws", newdata = nd, n_draws = n_draws, seed = 3)$y

  support <- 0:200
  pmf <- as.numeric(predict(m, type = "distribution", newdata = nd, support = support,
                            format = "matrix"))
  pmf <- c(pmf, 1 - sum(pmf))  # tail bin: y > max(support)

  # One bin per support value (0, 1, ..., 200) plus a tail bin (> 200):
  # length(support) + 1 bins, matching pmf's length.
  breaks <- c(support, max(support) + 1, Inf)
  observed <- table(cut(draws, breaks = breaks, right = FALSE))
  expect_length(observed, length(pmf))
  expected <- pmf * n_draws

  # Merge sparse bins (expected count < 5) into the tail bin, the usual
  # chi-square rule of thumb, so the test statistic isn't dominated by noise
  # in bins the model gives almost no mass.
  keep <- expected >= 5
  keep[length(keep)] <- TRUE  # always keep the tail bin
  obs_c <- c(observed[keep][-sum(keep)], sum(observed[!keep]) + observed[keep][sum(keep)])
  exp_c <- c(expected[keep][-sum(keep)], sum(expected[!keep]) + expected[keep][sum(keep)])
  exp_c <- exp_c * sum(obs_c) / sum(exp_c)  # renormalize for any floating-point slack

  chisq <- sum((obs_c - exp_c)^2 / exp_c)
  df <- length(obs_c) - 1
  pval <- stats::pchisq(chisq, df, lower.tail = FALSE)
  expect_true(pval > 0.001, info = paste("chisq =", chisq, "df =", df, "p =", pval))
})

test_that("keep = passes identifier columns through for quantile (vector) / exceedance / draws / distribution (round10 H.6)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  nd <- genevzinb2[1:5, ]
  nd$loc_id <- letters[1:5]

  qv <- suppressWarnings(predict(m, type = "quantile", quantile = c(.5, .9), newdata = nd, keep = "loc_id"))
  expect_equal(qv$loc_id, nd$loc_id)

  ev <- suppressWarnings(predict(m, type = "exceedance", threshold = c(0, 25), newdata = nd, keep = "loc_id"))
  expect_equal(ev$loc_id, nd$loc_id)

  dr <- predict(m, type = "draws", newdata = nd, n_draws = 3, seed = 1, keep = "loc_id")
  expect_equal(dr$loc_id, nd$loc_id[dr$.row])

  di <- predict(m, type = "distribution", newdata = nd, support = 0:5, keep = "loc_id")
  expect_equal(di$loc_id, nd$loc_id[di$.row])

  expect_error(predict(m, type = "exceedance", threshold = 0, newdata = nd, keep = "nope"),
              "not found")
})

test_that("evinf_pmf()/evinf_cdf() used by the new H.2/H.4 types are safe on bootstrap replicates via confint (regression for the H.1 is_zinb fix)", {
  # predict(..., confint = TRUE) iterates object$bootstraps internally; this
  # is a behavioral (not just unit-level) check that the round10 H.1 follow-up
  # fix (evinf_dist_params()'s explicit is_zinb =) is actually wired through
  # H.3/H.4's confint paths and not just exercised directly in test-evinf-pmf.R.
  m <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 5)
  expect_silent(suppressWarnings(
    predict(m, type = "quantile", quantile = c(.5, .9), confint = TRUE)
  ))
  expect_silent(suppressWarnings(
    predict(m, type = "exceedance", threshold = c(0, 25), confint = TRUE)
  ))
})
