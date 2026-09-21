# audit 4.10 - posterior-state tools

test_that("classify_states() returns the expected shape (evzinb)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  cs <- classify_states(m)
  expect_named(cs, c("prior_zero", "prior_count", "prior_evi",
                     "posterior_zero", "posterior_count", "posterior_evi",
                     "map_prior", "map_posterior", "y"))
  expect_equal(nrow(cs), nobs(m))
  expect_equal(levels(cs$map_prior), c("zero", "count", "evi"))
  expect_true(all(abs(rowSums(cs[, c("prior_zero", "prior_count", "prior_evi")]) - 1) < 1e-6))
})

test_that("classify_states() has no zero state for evinb", {
  mi <- fit_evinb_fast(bootstrap = FALSE)
  cs <- classify_states(mi)
  expect_true(all(cs$prior_zero == 0))
})

test_that("classify_states(newdata=) computes posteriors only with the response", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  data(genevzinb2, package = "evinf", envir = environment())
  nd <- genevzinb2[1:15, ]

  cs1 <- classify_states(m, newdata = nd)
  expect_true(any(is.finite(cs1$posterior_count)))

  cs2 <- classify_states(m, newdata = nd[, c("x1", "x2", "x3")])
  expect_true(all(is.na(cs2$posterior_count)))
  expect_true(all(is.na(cs2$map_posterior)))
})

test_that("threshold rule leaves low-confidence rows unassigned", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  cs <- classify_states(m, rule = "threshold", threshold = 0.99)
  expect_true(anyNA(cs$map_prior))
})

test_that("there is a single Pareto pmf, dpareto_disc(), not a duplicate copy (audit0.10 §1.12, D.4)", {
  expect_false(exists("evinf_pareto_pmf", where = asNamespace("evinf"), inherits = FALSE))

  y <- c(0L, 1L, 40L, 41L, 100L)
  C <- 37.4
  alpha <- 1.8
  prior <- matrix(c(0.2, 0.5, 0.3), nrow = length(y), ncol = 3, byrow = TRUE)
  resp <- evinf:::evinf_responsibilities(y, mu_nb = rep(5, length(y)),
                                         alpha_nb = 1, pl_alpha = alpha,
                                         C = C, prior = prior)
  d_evi <- evinf:::dpareto_disc(y, C, alpha)
  d_zero <- as.numeric(y == 0)
  d_count <- stats::dnbinom(y, mu = 5, size = 1)
  num_evi <- prior[, 3] * d_evi
  denom <- prior[, 1] * d_zero + prior[, 2] * d_count + num_evi
  expect_equal(unname(resp[, "evi"]), unname(num_evi / denom))
})

test_that("state_table() cross-tabulates prior vs posterior", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  st <- state_table(m)
  expect_s3_class(st, "evinf_state_table")
  expect_true(all(dim(st$counts) >= c(1, 1)))
  expect_output(print(st), "Row percentages")
})
