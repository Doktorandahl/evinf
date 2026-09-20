# Regenerates tests/testthat/fixtures/predict_snapshot_*.rds
#
# These freeze predict.evzinb()/predict.evinb() output for every type x pred
# x confint x return_bootstraps combination, generated from the code just
# before the evinf_predict_engine() merge (audit0.10 F.5, post-D.5 so it
# captures the conf_level = 0.95 default). test-predict-snapshot.R
# regenerates the same combinations from whatever code is current and checks
# they still match byte-for-byte (column names/order included), as an
# ongoing regression test, not just a one-time refactor check.
#
# Run once, from the package root, with the pre-refactor code checked out:
#   Rscript tests/testthat/fixtures/make_predict_snapshot.R
# then commit the .rds files. Do NOT regenerate them after the refactor,
# EXCEPT for a deliberate behavior change in predict() itself (as opposed to
# a refactor that should leave output untouched) -- e.g. round9 0.1 (review
# §2) intentionally changed harmonic/explog predictions wherever the fitted
# Pareto alpha collapses toward 0, from Inf / an astronomical finite number
# to a value clamped at alpha_pl_floor. In that case regenerate only
# predict_snapshot_pre_refactor.rds against the *unchanged* model fixtures
# (predict_snapshot_evzinb.rds / predict_snapshot_evinb.rds must not be
# refit), and check that every other combination still reproduces
# byte-for-byte before committing.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-evinf.R")

mz <- fit_evzinb_fast(bootstrap = TRUE)
mi <- fit_evinb_fast(bootstrap = TRUE)
saveRDS(mz, "tests/testthat/fixtures/predict_snapshot_evzinb.rds")
saveRDS(mi, "tests/testthat/fixtures/predict_snapshot_evinb.rds")

vector_types <- c("harmonic", "explog", "counts", "pareto_alpha", "evinf", "count_state")
zi_type <- "zi"

capture_all <- function(m, evzinb) {
  out <- list()
  types <- c(vector_types, if (evzinb) zi_type, "states", "all", "quantile")
  preds <- c("original", "bootstrap_median", "bootstrap_mean")

  for (ty in types) {
    for (pr in preds) {
      key_base <- paste(ty, pr, sep = "__")
      q_arg <- if (ty %in% c("quantile", "all")) 0.75 else NULL

      # confint = FALSE
      key <- paste0(key_base, "__confintFALSE")
      out[[key]] <- tryCatch(
        suppressWarnings(predict(m, type = ty, pred = pr, quantile = q_arg,
                                  confint = FALSE, multicore = FALSE)),
        error = function(e) paste("ERROR:", conditionMessage(e))
      )

      if (ty %in% c("states", "all")) next  # confint unsupported for these

      # confint = TRUE, return_bootstraps = FALSE
      key <- paste0(key_base, "__confintTRUE__rbFALSE")
      out[[key]] <- tryCatch(
        suppressWarnings(predict(m, type = ty, pred = pr, quantile = q_arg,
                                  confint = TRUE, return_bootstraps = FALSE,
                                  multicore = FALSE)),
        error = function(e) paste("ERROR:", conditionMessage(e))
      )

      # confint = TRUE, return_bootstraps = TRUE
      key <- paste0(key_base, "__confintTRUE__rbTRUE")
      out[[key]] <- tryCatch(
        suppressWarnings(predict(m, type = ty, pred = pr, quantile = q_arg,
                                  confint = TRUE, return_bootstraps = TRUE,
                                  multicore = FALSE)),
        error = function(e) paste("ERROR:", conditionMessage(e))
      )
    }
  }
  out
}

snap_z <- capture_all(mz, evzinb = TRUE)
snap_i <- capture_all(mi, evzinb = FALSE)

n_err_z <- sum(vapply(snap_z, function(x) is.character(x) && length(x) == 1 && startsWith(x, "ERROR:"), logical(1)))
n_err_i <- sum(vapply(snap_i, function(x) is.character(x) && length(x) == 1 && startsWith(x, "ERROR:"), logical(1)))
cat("evzinb combos:", length(snap_z), " errors:", n_err_z, "\n")
cat("evinb combos:", length(snap_i), " errors:", n_err_i, "\n")

saveRDS(list(evzinb = snap_z, evinb = snap_i),
        "tests/testthat/fixtures/predict_snapshot_pre_refactor.rds")
message("done")
