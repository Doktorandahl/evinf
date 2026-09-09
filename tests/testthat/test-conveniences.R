# audit 4.13 - gof_map_evinf(), bare `block`, and hks documentation.

test_that("gof_map_evinf() has a row for every glance() column", {
  gm <- gof_map_evinf()
  expect_s3_class(gm, "tbl_df")
  expect_identical(names(gm), c("raw", "clean", "fmt"))

  m <- fit_evzinb_fast(bootstrap = FALSE)
  glance_cols <- names(generics::glance(m))
  # every gof_map row refers to a real glance column
  expect_true(all(gm$raw %in% glance_cols))
  # and the key statistics are present
  expect_true(all(c("nobs", "npar", "aic", "bic", "logLik", "converged",
                    "n_above_c", "n_bootstraps", "n_failed_bootstraps") %in% gm$raw))
})

test_that("gof_map_evinf(extra = ) appends rows", {
  extra <- data.frame(raw = "n_em_steps", clean = "EM steps", fmt = 0)
  gm <- gof_map_evinf(extra = extra)
  expect_true("n_em_steps" %in% gm$raw)
  expect_equal(nrow(gm), nrow(gof_map_evinf()) + 1L)
})

test_that("the bundled gm_evzinb data object equals gof_map_evinf()", {
  data(gm_evzinb, package = "evinf", envir = environment())
  expect_equal(as.data.frame(gm_evzinb), as.data.frame(gof_map_evinf()))
})

test_that("block accepts a bare column name and a string equivalently", {
  d <- genevzinb2_factor()
  d$cluster <- rep(seq_len(nrow(d) / 5), each = 5)

  m_bare <- suppressMessages(evinf::evzinb(
    y ~ x1 + x2 + x3, data = d, control = .fast_control(),
    bootstrap = TRUE, n_bootstraps = 5, multicore = FALSE,
    boot_seed = 99, verbose = FALSE, block = cluster
  ))
  m_str <- suppressMessages(evinf::evzinb(
    y ~ x1 + x2 + x3, data = d, control = .fast_control(),
    bootstrap = TRUE, n_bootstraps = 5, multicore = FALSE,
    boot_seed = 99, verbose = FALSE, block = "cluster"
  ))
  expect_identical(m_bare$block, "cluster")
  expect_identical(m_str$block, "cluster")
  expect_equal(coef(m_bare, "all"), coef(m_str, "all"))
})

test_that("block accepts a string held in a variable", {
  d <- genevzinb2_factor()
  d$cluster <- rep(seq_len(nrow(d) / 5), each = 5)
  b <- "cluster"
  m <- suppressMessages(evinf::evzinb(
    y ~ x1 + x2 + x3, data = d, control = .fast_control(),
    bootstrap = FALSE, verbose = FALSE, block = b
  ))
  expect_identical(m$block, "cluster")
})

test_that("a bare block name is resolved to the data column, not a global (N5)", {
  d <- genevzinb2_factor()
  d$country <- rep(seq_len(nrow(d) / 5), each = 5)
  country <- c("shadowing", "global", "vector")  # would break eval_tidy()
  m <- suppressMessages(evinf::evzinb(
    y ~ x1 + x2 + x3, data = d, control = .fast_control(),
    bootstrap = FALSE, verbose = FALSE, block = country
  ))
  expect_identical(m$block, "country")
})

test_that("hks help lists exactly the columns that ship", {
  data(hks, package = "evinf", envir = environment())
  expect_setequal(
    names(hks),
    c("osvAll", "troopLag", "policeLag", "militaryobserversLag", "brv_AllLag_log",
      "osvAllLagDum", "incomp", "epduration", "lntpop", "troopLag_log",
      "epdur_log", "policeLag_log", "militaryobserversLag_log")
  )
})
