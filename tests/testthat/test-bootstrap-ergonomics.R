# audit 4.9 - bootstrap ergonomics

test_that("add_bootstraps() appends replicates and records the seed", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  m2 <- suppressWarnings(suppressMessages(add_bootstraps(m, 4, boot_seed = 7)))
  expect_length(m2$bootstraps, 9L)
  expect_equal(names(m2$bootstraps), paste0("bootstrap_", 1:9))
  expect_length(m2$boot_seeds, 2L)

  expect_error(add_bootstraps(fit_evzinb_fast(bootstrap = FALSE), 3), "no bootstraps")
})

test_that("add_bootstraps() with a seed is reproducible", {
  m <- fit_evzinb_fast(n_bootstraps = 3)
  a <- suppressWarnings(suppressMessages(add_bootstraps(m, 3, boot_seed = 123)))
  b <- suppressWarnings(suppressMessages(add_bootstraps(m, 3, boot_seed = 123)))
  expect_equal(a$bootstraps[[4]]$coef$Beta.NB, b$bootstraps[[4]]$coef$Beta.NB)
})

test_that("two fits with the same boot_seed give identical bootstrap coefficients", {
  data(genevzinb2, package = "evinf", envir = environment())
  ctl <- evinf_control(c.lim = c(50, 1000), init.C = 200)
  a <- suppressMessages(evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 4,
                               boot_seed = 55, verbose = FALSE, control = ctl))
  b <- suppressMessages(evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 4,
                               boot_seed = 55, verbose = FALSE, control = ctl))
  expect_equal(a$bootstraps[[2]]$coef$Beta.PL, b$bootstraps[[2]]$coef$Beta.PL)
})

test_that("failed_bootstraps() returns a tibble", {
  m <- fit_evzinb_fast(n_bootstraps = 5)
  fb <- failed_bootstraps(m)
  expect_s3_class(fb, "tbl_df")
  expect_named(fb, c("id", "message"))

  m0 <- fit_evzinb_fast(bootstrap = FALSE)
  expect_equal(nrow(failed_bootstraps(m0)), 0L)
})

test_that("object records the bootstrap seed", {
  m <- fit_evzinb_fast(n_bootstraps = 3)
  expect_equal(m$boot_seeds, list(123))
})
