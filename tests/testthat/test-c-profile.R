# audit 4.5 - threshold diagnostics

test_that("the fitted object carries the C_EV diagnostics", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  prof <- c_profile(m)
  expect_named(prof, c("c", "loglik"))
  expect_gt(nrow(prof), 1L)
  expect_true(all(is.finite(prof$loglik)))

  expect_type(m$c_trace, "double")
  expect_equal(m$n_above_c, sum(m$data$y >= m$coef$C))
  expect_gt(m$n_em_steps, 0L)
})

test_that("glance() reports n_above_c and n_em_steps", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  g <- glance(m)
  expect_true(all(c("n_above_c", "n_em_steps") %in% names(g)))
  expect_equal(g$n_above_c, m$n_above_c)
})

test_that("plot_c_profile() returns a ggplot", {
  skip_if_not_installed("ggplot2")
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_s3_class(plot_c_profile(m), "ggplot")
})

test_that("the fixed default c.lim still works when passed explicitly", {
  data(genevzinb2, package = "evinf", envir = environment())
  m <- suppressMessages(evzinb(y ~ x1 + x2 + x3, data = genevzinb2,
                               bootstrap = FALSE, verbose = FALSE,
                               control = evinf_control(c.lim = c(50, 1000))))
  expect_equal(m$control$c.lim, c(50, 1000))
  expect_false(isTRUE(m$c_lim_default))
})

test_that("genevzinb2 fits with the data-driven c.lim default (with a message)", {
  data(genevzinb2, package = "evinf", envir = environment())
  expect_message(
    m <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, bootstrap = FALSE,
                verbose = FALSE),
    "data-driven candidate range"
  )
  expect_true(isTRUE(m$c_lim_default))
  expect_length(m$control$c.lim, 2L)
  expect_true(isTRUE(m$converge))
})

test_that("bootstrap objects keep only c_trace, not the full profile", {
  m <- fit_evzinb_fast(n_bootstraps = 4)
  ok <- m$bootstraps[[1]]
  expect_null(ok$c_profile)
  expect_false(is.null(ok$c_trace))
})
