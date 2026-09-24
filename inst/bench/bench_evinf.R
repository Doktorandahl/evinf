# round8 A.4 (audit §5.4): before/after benchmark for work package A.
# round10 J.3 (audit §5.11, decision D8): switched from the bundled
# genevzinb2 / hks datasets to evinf_bench_data(n, seed), a synthetic
# generator with known true parameters -- lets this script (and anyone doing
# their own scale testing) ask for whatever n it needs without the installed
# package carrying an extra .rda per size. Also compares each timing against
# inst/bench/reference_timings.rds and prints (never fails the build on) a
# warning when a case is more than 2x slower than its stored reference.
#
# Tracked under inst/bench/ (round9 0.5, review §3) so it's reproducible from
# the repository and runnable in CI, rather than living only in the
# gitignored dev/ directory. The headline speedups it reports are
# BLAS-dependent (observed ~8.4x on Apple Accelerate vs. ~1.9x on
# Linux/OpenBLAS for the same fit) -- quote a range across platforms, not a
# single number, when citing these figures; reference_timings.rds was
# recorded on one contributor's machine as a starting baseline (see the
# header of that file) and should be refreshed from actual CI runs once a
# maintainer has a few of those to average.
#
# Run from the package root, with the current source tree loaded:
#
#   Rscript inst/bench/bench_evinf.R
#
# (also run automatically in CI -- see .github/workflows/benchmark.yaml).
#
# To get "before" numbers for work package A specifically, check out the
# commit just before round8 A.1 (069616e, "Recompute fitted y.hat.pl_* from
# the final props...", the last Part 0 commit) and re-run; the candidate-grid
# sizes reported for each case are what A.1 scales with (that comparison
# predates evinf_bench_data() -- use genevzinb2 / hks as that commit did).

suppressMessages(devtools::load_all(".", quiet = TRUE))

time_it <- function(expr) {
  t <- system.time(force(expr))
  unname(t["elapsed"])
}

ref_path <- file.path("inst", "bench", "reference_timings.rds")
reference <- if (file.exists(ref_path)) readRDS(ref_path) else NULL
results <- list()

# Records `elapsed` under `label`, prints it, and -- when a reference value
# exists for `label` -- compares against it: a diagnostic, never a build
# failure (no stop()/quit(status = 1) anywhere in this script).
report_timing <- function(label, elapsed, extra = "") {
  results[[label]] <<- elapsed
  ref <- reference[[label]]
  line <- sprintf("%s: %.4fs%s", label, elapsed, extra)
  if (!is.null(ref) && is.finite(ref) && ref > 0) {
    ratio <- elapsed / ref
    line <- sprintf("%s  [reference %.4fs, %.2fx]", line, ref, ratio)
    if (ratio > 2) {
      line <- paste0(line, "  *** SLOWER THAN 2x REFERENCE ***")
    }
  }
  cat(line, "\n")
}

cat("=== A.4 / J.3 benchmark ===\n\n")

## 1: small full fit (genevzinb2-scale) ---------------------------------------
d_small <- evinf_bench_data(100, seed = 1)
t_small <- time_it(
  m_small <- suppressMessages(suppressWarnings(
    evzinb(y ~ x1 + x2 + x3, data = d_small, bootstrap = FALSE, verbose = FALSE)
  ))
)
report_timing("small_fit", t_small, sprintf(" (n=100, c grid size %d)", nrow(m_small$c_profile)))

## 2: medium full fit (hks-scale) ----------------------------------------------
d_medium <- evinf_bench_data(3700, seed = 2)
t_medium <- time_it(
  m_medium <- suppressMessages(suppressWarnings(
    evzinb(y ~ x1 + x2 + x3, data = d_medium, bootstrap = FALSE, verbose = FALSE)
  ))
)
report_timing("medium_fit", t_medium, sprintf(" (n=3700, c grid size %d)", nrow(m_medium$c_profile)))

## 3: one em_profile_c() call on the medium case -------------------------------
x_obj_medium <- list(
  X.multinom.ZC = m_medium$data$x.multinom.zc,
  X.multinom.PL = m_medium$data$x.multinom.pl,
  X.NB          = m_medium$data$x.nb,
  X.PL          = m_medium$data$x.pl,
  offset.nb     = m_medium$offset_nb
)
cands_medium <- m_medium$c_profile$c
t_profile_medium <- time_it(
  evinf:::em_profile_c(m_medium$data$y, x_obj_medium, m_medium$coef, cands_medium)
)
report_timing("medium_em_profile_c", t_profile_medium,
             sprintf(" (grid size %d)", length(cands_medium)))

## 4: one em_step() / update_bfgs_fun() call on the medium case ----------------
ext_medium <- evinf:::em_extend_design(x_obj_medium, length(m_medium$data$y))
t_step_medium <- time_it(
  evinf:::em_step(m_medium$data$y, ext_medium, m_medium$coef, m_medium$control, fixed_zc = FALSE)
)
report_timing("medium_em_step", t_step_medium)

## 5: large-count case (20,000 rows) -------------------------------------------
d_large <- evinf_bench_data(20000, seed = 3)
t_large_fit <- time_it(
  m_large <- suppressMessages(suppressWarnings(
    evzinb(y ~ x1 + x2 + x3, data = d_large, bootstrap = FALSE, verbose = FALSE)
  ))
)
report_timing("large_fit", t_large_fit, sprintf(" (n=20000, c grid size %d)", nrow(m_large$c_profile)))

x_obj_large <- list(
  X.multinom.ZC = m_large$data$x.multinom.zc,
  X.multinom.PL = m_large$data$x.multinom.pl,
  X.NB          = m_large$data$x.nb,
  X.PL          = m_large$data$x.pl,
  offset.nb     = m_large$offset_nb
)
cands_large <- m_large$c_profile$c
t_profile_large <- time_it(
  evinf:::em_profile_c(m_large$data$y, x_obj_large, m_large$coef, cands_large)
)
report_timing("large_em_profile_c", t_profile_large, sprintf(" (grid size %d)", length(cands_large)))

cat("\n=== done ===\n")

if (is.null(reference)) {
  cat("\n(no reference_timings.rds found; nothing was compared. See that file's\n",
      "header for how to record one.)\n", sep = "")
}
