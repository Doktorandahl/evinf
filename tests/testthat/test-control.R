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
