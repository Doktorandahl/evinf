# audit 4.1 - standard S3 generics

test_that("coef() flattens with <component>_<term> names plus alpha_nb, c_ev", {
  m <- fit_evzinb_fast(n_bootstraps = 6)
  cf <- coef(m)
  expect_true(all(c("count_(Intercept)", "count_x1", "zero_x1", "evi_x1",
                    "pareto_x1", "alpha_nb", "c_ev") %in% names(cf)))
  expect_equal(unname(coef(m, "count")), unname(m$coef$Beta.NB))

  mi <- fit_evinb_fast(n_bootstraps = 6)
  expect_false(any(grepl("^zero_", names(coef(mi)))))
})

test_that("vcov() is the bootstrap covariance and errors without bootstraps", {
  m <- fit_evzinb_fast(n_bootstraps = 8)
  V <- vcov(m)
  expect_equal(rownames(V), names(coef(m)))
  expect_true(isSymmetric(unname(V)))

  m0 <- fit_evzinb_fast(bootstrap = FALSE)
  expect_error(vcov(m0), "bootstrap")
})

test_that("confint() returns percentile and approx intervals", {
  m <- fit_evzinb_fast(n_bootstraps = 10)
  ci <- confint(m)
  expect_equal(nrow(ci), length(coef(m)))
  expect_equal(ncol(ci), 2L)
  ci2 <- confint(m, parm = c("count_x1", "c_ev"), type = "approx")
  expect_equal(rownames(ci2), c("count_x1", "c_ev"))
})

test_that("logLik()/AIC()/BIC()/nobs() agree with the stored values", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  mi <- fit_evinb_fast(n_bootstraps = 5)
  expect_equal(as.numeric(logLik(m)), m$log.lik)
  expect_equal(attr(logLik(m), "df"), length(m$par.all))
  expect_equal(AIC(m), m$AIC)
  expect_equal(BIC(m), m$BIC)
  expect_equal(AIC(mi), mi$AIC)
  expect_equal(nobs(m), nrow(m$data$x.nb))
})

test_that("formula()/terms()/model.frame() dispatch per component", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_s3_class(formula(m, "zero"), "formula")
  expect_s3_class(terms(m, "pareto"), "terms")
  expect_identical(model.frame(m), m$data$data)
})

test_that("fitted() and residuals() work", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_equal(fitted(m), predict(m, type = "harmonic"))
  expect_equal(residuals(m), m$data$y - fitted(m))
  rq <- residuals(m, type = "quantile", seed = 1)
  expect_length(rq, nobs(m))
  expect_true(all(is.finite(rq)))
  # reproducible
  expect_equal(rq, residuals(m, type = "quantile", seed = 1))
})

test_that("simulate() returns a tidy data frame", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  s <- simulate(m, nsim = 3, seed = 1)
  expect_named(s, c("sim_1", "sim_2", "sim_3"))
  expect_equal(nrow(s), nobs(m))
  expect_equal(s, simulate(m, nsim = 3, seed = 1))
})

test_that("residuals(seed=) and simulate(seed=) restore the caller's RNG state (audit0.10 §1.14, D.6)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)

  set.seed(123)
  before <- .Random.seed
  invisible(residuals(m, type = "quantile", seed = 1))
  expect_identical(.Random.seed, before)

  set.seed(123)
  before2 <- .Random.seed
  invisible(simulate(m, seed = 1))
  expect_identical(.Random.seed, before2)
})

test_that("simulate() attaches a 'seed' attribute following stats::simulate's convention (audit0.10 §1.14, D.6)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)

  s <- simulate(m, seed = 1)
  expect_equal(attr(s, "seed"), structure(1, kind = as.list(RNGkind())))

  set.seed(99)
  before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  s2 <- simulate(m)
  expect_equal(attr(s2, "seed"), before)
})

test_that("update() edits per-component formulas and other arguments", {
  m <- fit_evzinb_fast(n_bootstraps = 3)
  m2 <- suppressWarnings(suppressMessages(
    update(m, formula_pareto = . ~ . - x3, evaluate = TRUE)
  ))
  expect_equal(length(attr(stats::terms(formula(m2, "pareto")), "term.labels")), 2L)
  m3 <- suppressWarnings(suppressMessages(update(m, bootstrap = FALSE)))
  expect_null(m3$bootstraps)
})

test_that("update() does not embed the data frame in the stored call (N6)", {
  data(genevzinb2, package = "evinf", envir = environment())
  ctl <- evinf_control(c.lim = c(50, 1000), init.C = 200)
  m <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = FALSE, verbose = FALSE,
    control = ctl)))

  expect_false(is.data.frame(m$call$data))

  u <- suppressWarnings(suppressMessages(
    update(m, formula_pareto. = . ~ . - x3)))
  expect_false(is.data.frame(u$call$data))
  expect_true(is.symbol(u$call$data))          # `genevzinb2` or `.evinf_data`
  expect_lt(as.numeric(object.size(u$call)), 1e4)

  # the refit reproduces a direct call with the reduced Pareto formula
  d <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + x2 + x3,
    formula_pareto = stats::update.formula(y ~ x1 + x2 + x3, . ~ . - x3),
    data = genevzinb2, bootstrap = FALSE, verbose = FALSE, control = m$control)))
  expect_equal(unname(coef(u, "all")), unname(coef(d, "all")), tolerance = 1e-8)
})

test_that("update(data = ) re-resolves a data-driven c.lim, but keeps a pinned one (audit0.10 §1.9)", {
  data(genevzinb2, package = "evinf", envir = environment())
  m <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = FALSE, verbose = FALSE
  )))
  expect_true(isTRUE(m$c_lim_default))

  new_data <- genevzinb2
  new_data$y <- new_data$y * 5  # shifts the data-driven range substantially

  expect_message(
    u <- suppressWarnings(update(m, data = new_data)),
    "data-driven candidate range"
  )
  expect_false(identical(u$control$c.lim, m$control$c.lim))
  expect_true(isTRUE(u$c_lim_default))

  # a user-pinned c.lim is kept as-is across a data change
  ctl <- evinf_control(c.lim = c(50, 1000), init.C = 200)
  m2 <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = FALSE, verbose = FALSE, control = ctl
  )))
  expect_false(isTRUE(m2$c_lim_default))
  u2 <- suppressWarnings(suppressMessages(update(m2, data = new_data)))
  expect_identical(u2$control$c.lim, c(50, 1000))
  expect_identical(u2$control$init.C, 200)
})

test_that("update()'s formula_nb. change does not propagate to an inherited component formula (audit0.10 §1.9)", {
  data(genevzinb2, package = "evinf", envir = environment())
  m <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + x2, data = genevzinb2, bootstrap = FALSE, verbose = FALSE  # formula_pareto NULL -> inherits from formula_nb
  )))
  expect_equal(sort(all.vars(formula(m, "pareto"))[-1]), c("x1", "x2"))

  u <- suppressWarnings(suppressMessages(update(m, formula_nb. = . ~ . + x3)))
  expect_true("x3" %in% all.vars(formula(u, "count")))
  # the pareto formula keeps its originally-inherited terms; formula_nb.'s
  # change does not propagate to it.
  expect_equal(sort(all.vars(formula(u, "pareto"))[-1]), c("x1", "x2"))
})

test_that("modelsummary works with a grouped shape (skip if absent)", {
  skip_if_not_installed("modelsummary")
  # modelsummary's grouped-shape output path needs broom internally, but
  # broom is only one of modelsummary's own Suggests, not guaranteed to be
  # installed alongside it (observed on CI's minimal package set).
  skip_if_not_installed("broom")
  m <- fit_evzinb_fast(n_bootstraps = 6)
  # As with nnet::multinom, modelsummary needs the y.level grouping in `shape`.
  expect_no_error(
    modelsummary::modelsummary(m, shape = term + y.level ~ model,
                               output = "data.frame")
  )
})
