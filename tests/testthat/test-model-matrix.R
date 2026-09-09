# Baseline coefficients captured from evinf 0.9.3 (model.frame()[, -1] design):
#   evzinb(y ~ x1 + x2 + x3, genevzinb2, bootstrap = FALSE)
BASE_NB <- c(`(Intercept)` = 2.8161580, x1 = 0.3976191, x2 = 0.7695594, x3 = 0.6497135)
BASE_LOGLIK <- -251.2682021

test_that("model.matrix refactor is numerically identical for plain covariates (audit 1.7)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_equal(m$coef$Beta.NB, BASE_NB, tolerance = 1e-5)
  expect_equal(m$log.lik, BASE_LOGLIK, tolerance = 1e-5)
})

test_that("in-formula transformations and interactions are supported and named (audit 1.7)", {
  m <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + log(abs(x2) + 1) + x1:x3,
    data = { data(genevzinb2, package = "evinf", envir = environment()); genevzinb2 },
    bootstrap = FALSE, verbose = FALSE
  )))
  expect_true(all(c("log(abs(x2) + 1)", "x1:x3") %in% names(m$coef$Beta.NB)))
})

test_that("factor covariates expand to dummies and predict() accepts a level subset (audit 1.7)", {
  d <- genevzinb2_factor()
  m <- suppressWarnings(suppressMessages(evzinb(y ~ x1 + g, data = d,
                                               bootstrap = FALSE, verbose = FALSE)))
  expect_equal(sum(grepl("^g", names(m$coef$Beta.NB))), 2L)

  nd <- d[d$g %in% c("a", "b"), ][1:5, ]
  p <- predict(m, newdata = nd)
  expect_length(p, 5L)
  expect_true(all(is.finite(p)))
})

test_that("predict(model) equals predict(model, newdata = model$data$data) (audit 1.7)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  for (ty in c("harmonic", "explog", "counts", "pareto_alpha", "zi", "evinf", "count_state")) {
    expect_equal(predict(m, type = ty),
                 predict(m, newdata = m$data$data, type = ty),
                 info = ty)
  }
  expect_equal(predict(m, type = "states"),
               predict(m, newdata = m$data$data, type = "states"))
})

test_that("user-supplied init.Beta.NB is actually used (audit 1.8 / R0.6)", {
  d <- { data(genevzinb2, package = "evinf", envir = environment()); genevzinb2 }
  # NB has 3 predictors (+ intercept = 4), EVI has 1: before the fix, a length-4
  # start was discarded because the guard compared against the EVI design size.
  start <- c(1.5, 0.2, 0.3, 0.4)
  m <- suppressMessages(evzinb(
    y ~ x1 + x2 + x3, formula_evi = ~x1, data = d,
    bootstrap = FALSE, verbose = FALSE,
    control = evinf_control(init.Beta.NB = start, c.lim = c(50, 1000))
  ))
  expect_equal(as.numeric(m$ini.val$Beta.NB), start)
})
