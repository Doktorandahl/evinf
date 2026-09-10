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
  m <- fit_evzinb_fast(n_bootstraps = 3)   # fitted with boot_seed = 123
  a <- suppressWarnings(suppressMessages(add_bootstraps(m, 3, boot_seed = 456)))
  b <- suppressWarnings(suppressMessages(add_bootstraps(m, 3, boot_seed = 456)))
  expect_equal(a$bootstraps[[4]]$coef$Beta.NB, b$bootstraps[[4]]$coef$Beta.NB)
})

test_that("add_bootstraps() rejects a seed the model already used (N4)", {
  m <- fit_evzinb_fast(n_bootstraps = 3)   # boot_seed = 123
  expect_error(add_bootstraps(m, 2, boot_seed = 123), "already used")

  m2 <- suppressWarnings(suppressMessages(add_bootstraps(m, 2, boot_seed = 456)))
  expect_error(add_bootstraps(m2, 2, boot_seed = 456), "already used")
  expect_error(add_bootstraps(m2, 2, boot_seed = 123), "already used")

  # a fresh seed appends replicates that differ from the first batch
  m3 <- suppressWarnings(suppressMessages(add_bootstraps(m, 3, boot_seed = 789)))
  first <- as.numeric(m3$bootstraps[[1]]$coef$Beta.NB)
  new   <- as.numeric(m3$bootstraps[[5]]$coef$Beta.NB)
  expect_false(isTRUE(all.equal(first, new)))
})

test_that("bootstrapped models always record a seed, even with boot_seed = NULL (N4)", {
  data(genevzinb2, package = "evinf", envir = environment())
  ctl <- evinf_control(c.lim = c(50, 1000), init.C = 200)
  m <- suppressWarnings(suppressMessages(evzinb(
    y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 3, verbose = FALSE,
    control = ctl)))
  expect_length(m$boot_seeds, 1L)
  expect_false(is.null(m$boot_seeds[[1]]))
  expect_true(is.numeric(m$boot_seeds[[1]]))
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
  expect_named(fb, c("id", "type", "message"))

  m0 <- fit_evzinb_fast(bootstrap = FALSE)
  expect_equal(nrow(failed_bootstraps(m0)), 0L)
})

test_that("object records the bootstrap seed", {
  m <- fit_evzinb_fast(n_bootstraps = 3)
  expect_equal(m$boot_seeds, list(123))
})
