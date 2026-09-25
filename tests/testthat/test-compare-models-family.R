# round11 C1/C2: compare_models() under offsets and each family. Weights are
# already covered thoroughly for the default and Poisson families in
# test-compare-models.R / test-compare-models-poisson.R (round10 0.2, round9
# E.3) -- these tests fill the two cells the dev/argument_propagation.md sweep
# found untested: an offset model, and a hurdle-family model.

test_that("compare_models() carries an offset into the NB/ZINB competitor formulas", {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(31)
  d <- genevzinb2
  d$ex <- runif(nrow(d), 0.5, 2)
  ctrl <- evinf_control(c.lim = c(50, 1000), init.C = 200,
                        max.diff.par = 1e-6, max.no.em.steps = 3000)

  m <- suppressMessages(suppressWarnings(evzinb(
    y ~ x1 + x2 + offset(log(ex)), data = d,
    n_bootstraps = 2, multicore = FALSE, verbose = FALSE, control = ctrl
  )))
  comp <- suppressWarnings(suppressMessages(compare_models(m, multicore = FALSE)))

  # the offset reached the competitor formulas (not silently dropped)
  expect_true("offset(log(ex))" %in% attr(stats::terms(formula(comp$nb$full_run)), "term.labels") ||
                grepl("offset\\(log\\(ex\\)\\)", deparse(formula(comp$nb$full_run))))

  hand_nb <- suppressWarnings(MASS::glm.nb(y ~ x1 + x2 + offset(log(ex)), data = d))
  expect_equal(unname(coef(comp$nb$full_run)), unname(coef(hand_nb)), tolerance = 1e-4)
})

test_that("compare_models() runs for a hurdle-family fit; the ZINB/ZIP baseline stays the standard toolkit model (round11 C1 gap)", {
  # round12 B1: the maintainer's resolution of the round11 C1 gap was to add
  # a hurdle competitor *alongside* the existing ZINB one (see
  # test-compare-models-hurdle.R), not to replace it -- so a hurdle-family
  # fit's zinb_comparison stays pscl::zeroinfl(), the standard-toolkit
  # baseline, regardless of whether `object` itself is a mixture or a hurdle
  # model.
  m <- fit_evzinb_fast(family = evinf_family(zero = "hurdle"), n_bootstraps = 2)
  comp <- suppressWarnings(suppressMessages(compare_models(m, multicore = FALSE)))

  expect_s3_class(comp$zinb$full_run, "zeroinfl")
  expect_s3_class(comp$hurdle$full_run, "hurdle")
})
