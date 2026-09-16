# audit 4.3 - evinf_control()

test_that("evinf_control() returns a validated classed list", {
  ctrl <- evinf_control(max.no.em.steps = 300, c.lim = c(20, 500))
  expect_s3_class(ctrl, "evinf_control")
  expect_equal(ctrl$max.no.em.steps, 300)
  expect_equal(ctrl$c.lim, c(20, 500))
  expect_null(ctrl$init.C)

  expect_error(evinf_control(max.diff.par = -1), "positive")
  expect_error(evinf_control(c.lim = 5), "c.lim")
  expect_error(evinf_control(c.lim = c(500, 20)), "c.lim")
  expect_error(evinf_control(prune.c.range = 2), "prune.c.range")
  expect_error(evinf_control(prune.c.range = 1), "prune.c.range")
  expect_error(evinf_control(pdf.pl.type = "nope"))
  expect_output(print(ctrl), "evinf_control")
})

test_that("control settings are respected and stored on the fit", {
  m <- fit_evzinb_fast(bootstrap = FALSE,
                       control = evinf_control(c.lim = c(50, 1000), init.C = 200,
                                               max.no.em.steps = 40))
  expect_s3_class(m$control, "evinf_control")
  expect_equal(m$control$max.no.em.steps, 40)
  expect_equal(m$control$c.lim, c(50, 1000))
})

# audit0.10 §1.9: response validation and prune.c.range edge cases.

test_that("evzinb()/evinb() reject a non-integer response", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$y <- d$y + 0.5
  expect_error(
    suppressMessages(evzinb(y ~ x1 + x2 + x3, data = d, bootstrap = FALSE, verbose = FALSE)),
    "non-integer"
  )
  expect_error(
    suppressMessages(evinb(y ~ x1 + x2 + x3, data = d, bootstrap = FALSE, verbose = FALSE)),
    "non-integer"
  )
})

test_that("evzinb()/evinb() reject a negative response", {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$y[1] <- -1
  expect_error(
    suppressMessages(evzinb(y ~ x1 + x2 + x3, data = d, bootstrap = FALSE, verbose = FALSE)),
    "negative"
  )
  expect_error(
    suppressMessages(evinb(y ~ x1 + x2 + x3, data = d, bootstrap = FALSE, verbose = FALSE)),
    "negative"
  )
})

test_that("em_c_candidates() prunes a 3-value grid without sample()'s length-1 trap", {
  y <- c(10, 20, 30)
  # sample.size works out to 1 (keep the only interior candidate): with the
  # unfixed sample(seq(2, 2), size = 1, ...), sample() reads the length-1
  # `2` as "sample from 1:2", not as the single candidate position 2.
  out <- evinf:::em_c_candidates(y, c.lim = c(5, 35), prune.c.range = 0)
  expect_identical(out, c(10, 20, 30))

  # fewer than 3 candidates: pruning is skipped entirely, whatever the setting.
  out2 <- evinf:::em_c_candidates(c(10, 20), c.lim = c(5, 35), prune.c.range = 0.9)
  expect_identical(out2, c(10, 20))
})

test_that("pruning is reproducible and leaves the caller's .Random.seed untouched", {
  y <- seq(1, 100)
  set.seed(42)
  invisible(runif(1))  # make sure .Random.seed exists and has some state
  seed_before <- .Random.seed

  out1 <- evinf:::em_c_candidates(y, c.lim = c(1, 100), prune.c.range = 0.5)
  out2 <- evinf:::em_c_candidates(y, c.lim = c(1, 100), prune.c.range = 0.5)

  expect_identical(out1, out2)
  expect_identical(.Random.seed, seed_before)
})

test_that("deprecated individual arguments still work but warn once", {
  data(genevzinb2, package = "evinf", envir = environment())
  expect_warning(
    m <- suppressMessages(evzinb(y ~ x1 + x2 + x3, data = genevzinb2,
                                 bootstrap = FALSE, verbose = FALSE,
                                 max.no.em.steps = 42, c.lim = c(50, 1000))),
    "deprecated"
  )
  expect_equal(m$control$max.no.em.steps, 42)
  expect_equal(m$control$c.lim, c(50, 1000))
})
