# round8 A.4 (audit §5.4): before/after benchmark for work package A.
#
# Tracked under inst/bench/ (round9 0.5, review §3) so it's reproducible from
# the repository and runnable in CI, rather than living only in the
# gitignored dev/ directory. The headline speedups it reports are
# BLAS-dependent (observed ~8.4x on Apple Accelerate vs. ~1.9x on
# Linux/OpenBLAS for the same hks fit) -- quote a range across platforms, not
# a single number, when citing these figures.
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
# sizes reported for each case are what A.1 scales with.

suppressMessages(devtools::load_all(".", quiet = TRUE))

time_it <- function(expr) {
  t <- system.time(force(expr))
  unname(t["elapsed"])
}

data(genevzinb2, package = "evinf")
data(hks, package = "evinf")

f_hks <- osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
  lntpop + brv_AllLag_log + osvAllLagDum + incomp
f_hks_pareto <- ~ troopLag_log + epdur_log + brv_AllLag_log + lntpop

ctrl_gen <- evinf_control(c.lim = c(50, 1000), init.C = 200)

cat("=== A.4 benchmark ===\n\n")

## 1: full fit on genevzinb2 -------------------------------------------------
t_gen <- time_it(
  m_gen <- suppressMessages(suppressWarnings(
    evzinb(y ~ x1 + x2 + x3, data = genevzinb2, control = ctrl_gen,
          verbose = FALSE, bootstrap = FALSE)
  ))
)
cat(sprintf("genevzinb2 full fit: %.3fs (c grid size %d)\n",
            t_gen, nrow(m_gen$c_profile)))

## 2: full fit on hks ---------------------------------------------------------
t_hks <- time_it(
  m_hks <- suppressMessages(suppressWarnings(
    evzinb(f_hks, formula_pareto = f_hks_pareto, data = hks,
          bootstrap = FALSE, verbose = FALSE)
  ))
)
cat(sprintf("hks full fit: %.3fs (c grid size %d)\n",
            t_hks, nrow(m_hks$c_profile)))

## 3: one em_profile_c() call on hks ------------------------------------------
x_obj_hks <- list(
  X.multinom.ZC = m_hks$data$x.multinom.zc,
  X.multinom.PL = m_hks$data$x.multinom.pl,
  X.NB          = m_hks$data$x.nb,
  X.PL          = m_hks$data$x.pl,
  offset.nb     = m_hks$offset_nb
)
cands_hks <- m_hks$c_profile$c
t_profile_hks <- time_it(
  evinf:::em_profile_c(m_hks$data$y, x_obj_hks, m_hks$coef, cands_hks)
)
cat(sprintf("hks em_profile_c() (grid size %d): %.4fs\n",
            length(cands_hks), t_profile_hks))

## 4: one em_step() / update_bfgs_fun() call on hks --------------------------
ext_hks <- evinf:::em_extend_design(x_obj_hks, length(m_hks$data$y))
t_step_hks <- time_it(
  evinf:::em_step(m_hks$data$y, ext_hks, m_hks$coef, m_hks$control, fixed_zc = FALSE)
)
cat(sprintf("hks em_step(): %.4fs\n", t_step_hks))

## 5: large-count case (~20,000 rows, simulated from the hks fit) ------------
set.seed(1)
big_data <- hks[sample(nrow(hks), 20000, replace = TRUE), ]
big_data$osvAll <- simulate(m_hks, newdata = big_data, seed = 1)$sim_1

t_big_fit <- time_it(
  m_big <- suppressMessages(suppressWarnings(
    evzinb(f_hks, formula_pareto = f_hks_pareto, data = big_data,
          bootstrap = FALSE, verbose = FALSE)
  ))
)
cat(sprintf("large-count (n=20000) full fit: %.3fs (c grid size %d)\n",
            t_big_fit, nrow(m_big$c_profile)))

x_obj_big <- list(
  X.multinom.ZC = m_big$data$x.multinom.zc,
  X.multinom.PL = m_big$data$x.multinom.pl,
  X.NB          = m_big$data$x.nb,
  X.PL          = m_big$data$x.pl,
  offset.nb     = m_big$offset_nb
)
cands_big <- m_big$c_profile$c
t_profile_big <- time_it(
  evinf:::em_profile_c(m_big$data$y, x_obj_big, m_big$coef, cands_big)
)
cat(sprintf("large-count (n=20000) em_profile_c() (grid size %d): %.4fs\n",
            length(cands_big), t_profile_big))

cat("\n=== done ===\n")
