# audit 4.4 - offset in the count (NB) component

make_offset_data <- function() {
  data(genevzinb2, package = "evinf", envir = environment())
  d <- genevzinb2
  d$zoff <- rep(0, nrow(d))
  d$loff <- rep(log(2), nrow(d))
  d
}

tight_ctrl <- function(...) {
  evinf_control(c.lim = c(50, 1000), init.C = 200,
                max.diff.par = 1e-4, max.no.em.steps = 2000, ...)
}

test_that("offset(rep(0, n)) reproduces the no-offset fit (audit 4.4)", {
  d <- make_offset_data()
  m0 <- suppressMessages(evzinb(y ~ x1 + x2 + x3, data = d, bootstrap = FALSE,
                                verbose = FALSE, control = tight_ctrl()))
  mz <- suppressMessages(evzinb(y ~ x1 + x2 + x3 + offset(zoff), data = d,
                                bootstrap = FALSE, verbose = FALSE,
                                control = tight_ctrl()))
  expect_equal(unname(mz$coef$Beta.NB), unname(m0$coef$Beta.NB), tolerance = 1e-6)
  expect_equal(mz$log.lik, m0$log.lik, tolerance = 1e-6)
})

test_that("a constant log(2) offset shifts the NB intercept by -log(2) (audit 4.4)", {
  d <- make_offset_data()
  m0 <- suppressMessages(evzinb(y ~ x1 + x2 + x3, data = d, bootstrap = FALSE,
                                verbose = FALSE, control = tight_ctrl()))
  ml <- suppressMessages(evzinb(y ~ x1 + x2 + x3 + offset(loff), data = d,
                                bootstrap = FALSE, verbose = FALSE,
                                control = tight_ctrl()))
  expect_equal(ml$coef$Beta.NB[["(Intercept)"]],
               m0$coef$Beta.NB[["(Intercept)"]] - log(2), tolerance = 5e-3)
  expect_equal(unname(ml$coef$Beta.NB[-1]), unname(m0$coef$Beta.NB[-1]),
               tolerance = 1e-3)
  expect_true(isTRUE(ml$has_offset))
})

test_that("predict(newdata=) needs the offset column (audit 4.4)", {
  d <- make_offset_data()
  ml <- suppressMessages(evzinb(y ~ x1 + x2 + x3 + offset(loff), data = d,
                                bootstrap = FALSE, verbose = FALSE,
                                control = tight_ctrl()))
  nd <- d[1:5, ]
  expect_length(predict(ml, newdata = nd, type = "counts"), 5L)

  # doubling the offset column doubles the predicted NB count
  nd2 <- nd; nd2$loff <- nd$loff + log(2)
  expect_equal(unname(predict(ml, newdata = nd2, type = "counts")),
               unname(2 * predict(ml, newdata = nd, type = "counts")),
               tolerance = 1e-6)

  expect_error(
    predict(ml, newdata = nd[, c("y", "x1", "x2", "x3")], type = "counts"),
    "offset"
  )
})

test_that("offset() in a non-count component is an error (audit 4.4)", {
  d <- make_offset_data()
  expect_error(
    suppressMessages(evzinb(y ~ x1, formula_zi = ~ x1 + offset(loff), data = d,
                            bootstrap = FALSE)),
    "count component"
  )
  expect_error(
    suppressMessages(evinb(y ~ x1, formula_pareto = ~ x1 + offset(loff), data = d,
                           bootstrap = FALSE)),
    "count component"
  )
})

test_that("bootstrapping works with an offset (audit 4.4)", {
  d <- make_offset_data()
  m <- suppressMessages(evzinb(y ~ x1 + x2 + x3 + offset(zoff), data = d,
                               n_bootstraps = 3, boot_seed = 1, verbose = FALSE,
                               control = evinf_control(c.lim = c(50, 1000),
                                                       init.C = 200)))
  expect_length(m$bootstraps, 3L)
})
