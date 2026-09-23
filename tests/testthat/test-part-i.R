# round10 Part I (audit §5.9): ecosystem -- marginaleffects extensions (I.1),
# augment() (I.2), emmeans (I.3), print.evzinbcomp() (I.4). I.5 (skip texreg)
# is a NEWS-only decision, no test.

# --- I.1: marginaleffects ---------------------------------------------------

test_that("get_predict(type = 'states') returns a long group/estimate frame matching predict(type = 'states')", {
  skip_if_not_installed("marginaleffects")
  m <- fit_evzinb_fast(bootstrap = FALSE)
  data(genevzinb2, package = "evinf", envir = environment())
  st <- predict(m, type = "states")

  gp <- marginaleffects::get_predict(m, newdata = genevzinb2, type = "states")
  expect_setequal(unique(gp$group), c("zero", "count", "evi"))
  expect_equal(nrow(gp), nrow(genevzinb2) * 3L)
  for (g in c("zero", "count", "evi")) {
    sub <- gp[gp$group == g, ]
    expect_equal(sub$estimate, st[[paste0("pr_", g)]], tolerance = 1e-10)
  }

  mi <- fit_evinb_fast(bootstrap = FALSE)
  gpi <- marginaleffects::get_predict(mi, newdata = genevzinb2, type = "states")
  expect_setequal(unique(gpi$group), c("count", "evi"))
})

test_that("get_predict(type = 'quantile') uses the continuous surrogate and avg_predictions() runs", {
  skip_if_not_installed("marginaleffects")
  m <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 5)
  data(genevzinb2, package = "evinf", envir = environment())

  gp <- marginaleffects::get_predict(m, newdata = genevzinb2, type = "quantile")
  cont <- quantiles_from_evzinb(m, 0.5, newdata = genevzinb2, round = FALSE)
  expect_equal(gp$estimate, as.numeric(cont), tolerance = 1e-10)

  gp90 <- marginaleffects::get_predict(m, newdata = genevzinb2, type = "quantile", quantile = 0.9)
  cont90 <- quantiles_from_evzinb(m, 0.9, newdata = genevzinb2, round = FALSE)
  expect_equal(gp90$estimate, as.numeric(cont90), tolerance = 1e-10)

  # round11 A3: avg_predictions()/predictions() route `type` through
  # marginaleffects' own sanitize_type(), which rejects any type not in its
  # internal dictionary for a model's class unless that class is completely
  # absent from the dictionary -- true from 0.22.0 on (bisected against the
  # CRAN archive: 0.21.0 errors with "Must be element of set {'response',
  # 'class','link'}, but is 'quantile'", 0.22.0 does not; NEWS.md does not
  # call the fix out explicitly). get_predict() itself (tested above) never
  # goes through sanitize_type() and always worked, at any version.
  skip_if_not_installed("marginaleffects", "0.22.0")
  ap <- suppressWarnings(marginaleffects::avg_predictions(m, type = "quantile"))
  expect_true(is.finite(ap$estimate))
})

test_that("get_predict(type = 'exceedance') matches predict() (single and several thresholds)", {
  skip_if_not_installed("marginaleffects")
  m <- fit_evzinb_fast(bootstrap = FALSE)
  data(genevzinb2, package = "evinf", envir = environment())

  gp1 <- marginaleffects::get_predict(m, newdata = genevzinb2, type = "exceedance", threshold = 25)
  p1 <- predict(m, newdata = genevzinb2, type = "exceedance", threshold = 25)
  expect_equal(gp1$estimate, p1$p_ge_25, tolerance = 1e-10)

  gp2 <- marginaleffects::get_predict(m, newdata = genevzinb2, type = "exceedance", threshold = c(0, 25))
  expect_setequal(unique(gp2$group), c("p_ge_0", "p_ge_25"))
  expect_true(all(gp2$estimate[gp2$group == "p_ge_0"] == 1))

  expect_error(marginaleffects::get_predict(m, newdata = genevzinb2, type = "exceedance"),
              "threshold")
})

# --- I.2: augment() ----------------------------------------------------------

test_that("augment.evzinb() reproduces classify_states()/predict()/residuals() exactly, on the training data", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  a <- augment(m, seed = 5)

  data(genevzinb2, package = "evinf", envir = environment())
  expect_equal(nrow(a), nrow(genevzinb2))
  expect_equal(a$.fitted, as.numeric(predict(m, type = "harmonic")), tolerance = 1e-10)

  cs <- classify_states(m)
  expect_equal(a$.prob_zero, cs$prior_zero)
  expect_equal(a$.prob_count, cs$prior_count)
  expect_equal(a$.prob_evi, cs$prior_evi)
  expect_equal(a$.state, cs$map_prior)
  expect_equal(a$.post_zero, cs$posterior_zero)
  expect_equal(a$.post_state, cs$map_posterior)

  r <- residuals(m, type = "quantile", seed = 5)
  expect_equal(a$.resid, r)
})

test_that("augment.evinb() drops the zero columns", {
  mi <- fit_evinb_fast(bootstrap = FALSE)
  ai <- augment(mi)
  expect_false(any(c(".prob_zero", ".post_zero") %in% names(ai)))
  expect_true(all(c(".prob_count", ".prob_evi", ".post_count", ".post_evi") %in% names(ai)))
})

test_that("augment(newdata = ) without the response omits .resid/.post_*; with it, includes them", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  data(genevzinb2, package = "evinf", envir = environment())
  nd_no_y <- genevzinb2[1:5, setdiff(names(genevzinb2), "y")]

  a1 <- augment(m, newdata = nd_no_y)
  expect_equal(nrow(a1), 5L)
  expect_false(any(c(".resid", ".post_zero", ".post_state") %in% names(a1)))
  expect_true(all(c(".fitted", ".prob_zero", ".state") %in% names(a1)))

  a2 <- augment(m, newdata = genevzinb2[1:5, ])
  expect_true(all(c(".resid", ".post_zero", ".post_state") %in% names(a2)))
})

test_that("residuals(newdata = ) matches the training-data path when given the same rows (round10 I.2 extension)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  data(genevzinb2, package = "evinf", envir = environment())
  r1 <- residuals(m, type = "quantile", seed = 9)
  r2 <- residuals(m, type = "quantile", seed = 9, newdata = genevzinb2)
  expect_equal(r1, r2)
  r1r <- residuals(m, type = "response")
  r2r <- residuals(m, type = "response", newdata = genevzinb2)
  expect_equal(r1r, r2r)
})

# --- I.3: emmeans ------------------------------------------------------------

test_that("emmeans() on the count component matches predict(type = 'counts') on both the link and response scale", {
  skip_if_not_installed("emmeans")
  m <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 5)
  at <- list(x1 = c(0, 1), x2 = 0, x3 = 0)

  em_resp <- summary(emmeans::emmeans(m, ~x1, at = at, type = "response"))
  nd <- data.frame(x1 = c(0, 1), x2 = 0, x3 = 0)
  expect_equal(em_resp$rate %||% em_resp$response %||% em_resp$count,
              as.numeric(predict(m, newdata = nd, type = "counts")), tolerance = 1e-8)

  em_link <- summary(emmeans::emmeans(m, ~x1, at = at))
  expect_equal(em_link$emmean, log(as.numeric(predict(m, newdata = nd, type = "counts"))),
              tolerance = 1e-8)

  mi <- fit_evinb_fast(bootstrap = TRUE, n_bootstraps = 5)
  em_resp_i <- summary(emmeans::emmeans(mi, ~x1, at = at, type = "response"))
  nm <- intersect(c("rate", "response", "count"), names(em_resp_i))[1]
  expect_equal(em_resp_i[[nm]], as.numeric(predict(mi, newdata = nd, type = "counts")),
              tolerance = 1e-8)
})

test_that("emmeans(component = 'states') errors with a clear message (round10 I.3: state basis not implemented)", {
  skip_if_not_installed("emmeans")
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_error(
    suppressMessages(emmeans::emmeans(m, ~x1, at = list(x1 = 0, x2 = 0, x3 = 0), component = "states")),
    "not implemented|count component only"
  )
})

# --- I.4: print.evzinbcomp() -------------------------------------------------

test_that("print.evzinbcomp() shows a glance() table and the compare_fit() comparison", {
  m <- fit_evzinb_fast(bootstrap = TRUE, n_bootstraps = 5)
  cmp <- suppressMessages(suppressWarnings(compare_models(m, poisson_comparison = FALSE, zip_comparison = FALSE)))
  out <- capture.output(print(cmp))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Fit summary")
  expect_match(txt, "logLik")
  expect_match(txt, "Paired bootstrap fit comparison")
  expect_match(txt, "aic")
  expect_true(length(out) < 35)
})
