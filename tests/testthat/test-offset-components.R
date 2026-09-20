# round9 D.1 (audit §5.6): offset() in the zero-inflation and EVI-inflation
# components. See test-offset.R for the pre-existing count-component offset
# tests and the Pareto-component offset error.

make_offset_data2 <- function() {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$loff <- rep(log(2), nrow(d))
  d
}

tight_ctrl2 <- function(...) {
  evinf_control(c.lim = c(50, 1000), init.C = 200,
                max.diff.par = 1e-4, max.no.em.steps = 2000, ...)
}

test_that("offset() is accepted in formula_zi and formula_evi", {
  d <- make_offset_data2()
  m <- suppressMessages(evzinb(
    y ~ x1 + x2 + x3, formula_zi = ~ x1 + offset(loff),
    formula_evi = ~ x1 + offset(loff), data = d,
    bootstrap = FALSE, verbose = FALSE, control = tight_ctrl2()
  ))
  expect_s3_class(m, "evzinb")
  expect_true(is.finite(m$log.lik))
  expect_equal(length(m$offset_zc), nrow(d))
  expect_equal(length(m$offset_pl_mult), nrow(d))
  expect_true(all(m$offset_zc == log(2)))
  expect_true(all(m$offset_pl_mult == log(2)))
})

test_that("offset() is accepted in evinb()'s formula_evi", {
  d <- make_offset_data2()
  m <- suppressMessages(evinb(
    y ~ x1 + x2 + x3, formula_evi = ~ x1 + offset(loff), data = d,
    bootstrap = FALSE, verbose = FALSE, control = tight_ctrl2()
  ))
  expect_s3_class(m, "evinb")
  expect_true(is.finite(m$log.lik))
  expect_true(all(m$offset_pl_mult == log(2)))
})

test_that("an EVI offset shifts predict(type = 'evinf') in the expected direction", {
  # Isolate the offset's marginal effect on one fitted model's coefficients by
  # varying only the offset value in newdata (mirrors the existing NB-offset
  # test's "doubling the offset column doubles the predicted count" style).
  # The EVI logit is the log-odds of the EVI state against the count state,
  # so pushing the offset up must push pr_evi up.
  d <- make_offset_data2()
  m <- suppressMessages(evzinb(
    y ~ x1 + x2 + x3, formula_evi = ~ x1 + offset(loff), data = d,
    bootstrap = FALSE, verbose = FALSE, control = tight_ctrl2()
  ))
  nd <- d[1:5, ]
  base_p <- predict(m, newdata = nd, type = "evinf")

  nd_bigger <- nd
  nd_bigger$loff <- nd$loff + 5
  bigger_p <- predict(m, newdata = nd_bigger, type = "evinf")

  expect_true(all(bigger_p > base_p))
})

test_that("a zero-inflation offset shifts predict(type = 'zi') in the expected direction", {
  d <- make_offset_data2()
  m <- suppressMessages(evzinb(
    y ~ x1 + x2 + x3, formula_zi = ~ x1 + offset(loff), data = d,
    bootstrap = FALSE, verbose = FALSE, control = tight_ctrl2()
  ))
  nd <- d[1:5, ]
  base_p <- predict(m, newdata = nd, type = "zi")

  nd_bigger <- nd
  nd_bigger$loff <- nd$loff + 5
  bigger_p <- predict(m, newdata = nd_bigger, type = "zi")

  expect_true(all(bigger_p > base_p))
})

test_that("predict(newdata = ) requires and uses the EVI offset column", {
  d <- make_offset_data2()
  m <- suppressMessages(evzinb(
    y ~ x1 + x2 + x3, formula_evi = ~ x1 + offset(loff), data = d,
    bootstrap = FALSE, verbose = FALSE, control = tight_ctrl2()
  ))
  nd <- d[1:5, ]
  expect_length(predict(m, newdata = nd, type = "evinf"), 5L)
  expect_error(
    predict(m, newdata = nd[, c("y", "x1", "x2", "x3")], type = "evinf"),
    "offset"
  )
})

test_that("a defaulted zi/evi formula does not inherit formula_nb's offset (no inheritance)", {
  d <- make_offset_data2()
  m <- suppressMessages(evzinb(
    y ~ x1 + x2 + x3 + offset(loff), data = d,
    bootstrap = FALSE, verbose = FALSE, control = tight_ctrl2()
  ))
  expect_null(attr(stats::terms(m$formulas$formula_zi), "offset"))
  expect_null(attr(stats::terms(m$formulas$formula_evi), "offset"))
  expect_true(all(m$offset_zc == 0))
  expect_true(all(m$offset_pl_mult == 0))

  mi <- suppressMessages(evinb(
    y ~ x1 + x2 + x3 + offset(loff), data = d,
    bootstrap = FALSE, verbose = FALSE, control = tight_ctrl2()
  ))
  expect_null(attr(stats::terms(mi$formulas$formula_evi), "offset"))
  expect_true(all(mi$offset_pl_mult == 0))
})

test_that("bootstrap replicates index the EVI offset by boot_id, not by position", {
  # A fit with the offset column permuted to different rows, refit with the
  # same boot_seed (so the *same* boot_id draws are used), must give
  # different bootstrap coefficients if the offset is genuinely carried
  # through and indexed alongside the resampled rows.
  d <- make_offset_data2()
  set.seed(99)
  d$evoff <- stats::rnorm(nrow(d))

  m1 <- suppressMessages(evzinb(
    y ~ x1 + x2 + x3, formula_evi = ~ x1 + offset(evoff), data = d,
    n_bootstraps = 3, boot_seed = 7, multicore = FALSE, verbose = FALSE,
    control = tight_ctrl2()
  ))

  d2 <- d
  d2$evoff <- sample(d$evoff)
  m2 <- suppressMessages(evzinb(
    y ~ x1 + x2 + x3, formula_evi = ~ x1 + offset(evoff), data = d2,
    n_bootstraps = 3, boot_seed = 7, multicore = FALSE, verbose = FALSE,
    control = tight_ctrl2()
  ))

  c1 <- unname(m1$bootstraps[[1]]$coef$Beta.multinom.PL)
  c2 <- unname(m2$bootstraps[[1]]$coef$Beta.multinom.PL)
  expect_false(isTRUE(all.equal(c1, c2)))
})

test_that("default (no-offset) fits are numerically unchanged by D.1 (identity check)", {
  m <- fit_evzinb_fast(bootstrap = FALSE)
  expect_true(all(m$offset_zc == 0))
  expect_true(all(m$offset_pl_mult == 0))
  mi <- fit_evinb_fast(bootstrap = FALSE)
  expect_true(all(mi$offset_pl_mult == 0))
})
