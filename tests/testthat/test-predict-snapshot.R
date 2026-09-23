# audit0.10 F.5: predict.evzinb() / predict.evinb() were merged into a shared
# evinf_predict_engine(). This test regenerates every type x pred x confint x
# return_bootstraps combination captured by
# fixtures/make_predict_snapshot.R (frozen, pre-refactor, post-D.5 so it
# reflects the conf_level = 0.95 default) from the fixed fixture models
# below, and checks the current code reproduces it byte-identically --
# column names and order included -- protecting against any behavior change
# from the refactor, now and in the future.

skip_if_not(
  file.exists(test_path("fixtures", "predict_snapshot_pre_refactor.rds")),
  "predict snapshot fixture missing"
)

mz <- readRDS(test_path("fixtures", "predict_snapshot_evzinb.rds"))
mi <- readRDS(test_path("fixtures", "predict_snapshot_evinb.rds"))
snap <- readRDS(test_path("fixtures", "predict_snapshot_pre_refactor.rds"))

# round11 A4: "explog" is excluded from this byte-identical check. The
# fixture froze the *pre-A4* explog values, which were silently clamped by
# the global alpha_pl_floor (0.01) -- exactly the bug A4 fixes by making
# clamp_alpha_pl = FALSE (the new default) use alpha_pl as-is. That is an
# intentional behavior change, not a refactor regression, so the fixture is
# expected to disagree on this type; see test-alpha-pl-clamp.R for explog's
# own (post-A4) behavior contract.
vector_types <- c("harmonic", "counts", "pareto_alpha", "evinf", "count_state")

capture_one <- function(m, ty, pr, confint, return_bootstraps) {
  q_arg <- if (ty %in% c("quantile", "all")) 0.75 else NULL
  res <- tryCatch(
    suppressWarnings(predict(m, type = ty, pred = pr, quantile = q_arg,
                              confint = confint, return_bootstraps = return_bootstraps,
                              multicore = FALSE)),
    error = function(e) paste("ERROR:", conditionMessage(e))
  )
  # round11 A4: the fixture predates the "clamp_alpha_pl" result attribute
  # and, for type = "all", the (now intentionally different) explog column --
  # strip both before comparing so this test stays about the F.5 refactor,
  # not A4's behavior change (see the vector_types comment above).
  strip_clamp_attr <- function(x) {
    if (!is.null(attr(x, "clamp_alpha_pl"))) attr(x, "clamp_alpha_pl") <- NULL
    x
  }
  res <- strip_clamp_attr(res)
  if (is.list(res) && !is.data.frame(res) && !is.null(res$ci)) {
    res$ci <- strip_clamp_attr(res$ci)
  }
  if (identical(ty, "all")) {
    if (is.data.frame(res)) {
      res$explog <- NULL
    } else if (is.list(res) && !is.null(res$ci)) {
      res$ci$explog <- NULL
    }
  }
  res
}

check_class <- function(m, evzinb) {
  snapshot <- if (evzinb) snap$evzinb else snap$evinb
  types <- c(vector_types, if (evzinb) "zi", "states", "all", "quantile")
  preds <- c("original", "bootstrap_median", "bootstrap_mean")

  for (ty in types) {
    for (pr in preds) {
      key_base <- paste(ty, pr, sep = "__")

      key <- paste0(key_base, "__confintFALSE")
      expected <- snapshot[[key]]
      if (identical(ty, "all") && is.data.frame(expected)) {
        expected$explog <- NULL
      }
      expect_equal(capture_one(m, ty, pr, FALSE, FALSE), expected, info = key)

      if (ty %in% c("states", "all")) next

      key <- paste0(key_base, "__confintTRUE__rbFALSE")
      expect_equal(capture_one(m, ty, pr, TRUE, FALSE), snapshot[[key]], info = key)

      key <- paste0(key_base, "__confintTRUE__rbTRUE")
      expect_equal(capture_one(m, ty, pr, TRUE, TRUE), snapshot[[key]], info = key)
    }
  }
}

test_that("predict.evzinb() matches the pre-refactor snapshot for every combination (audit0.10 F.5)", {
  check_class(mz, evzinb = TRUE)
})

test_that("predict.evinb() matches the pre-refactor snapshot for every combination (audit0.10 F.5)", {
  check_class(mi, evzinb = FALSE)
})
