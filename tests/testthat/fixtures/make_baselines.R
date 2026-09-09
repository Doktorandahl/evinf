# Regenerates tests/testthat/fixtures/em_baseline_*.rds
#
# These capture the estimation output of the CURRENT code so that the
# R/a019.R -> R/em_*.R refactor (round 3, Part B) can be shown to be
# numerically identical (see tests/testthat/test-em-identity.R).
#
# Run once, from the package root, with the pre-refactor code checked out:
#   Rscript tests/testthat/fixtures/make_baselines.R
# then commit the .rds files. Do NOT regenerate them after the refactor.

suppressMessages(devtools::load_all(".", quiet = TRUE))

fixtures_dir <- "tests/testthat/fixtures"

# The quantities the refactor must preserve, pulled off a fitted evzinb/evinb.
capture_fit <- function(m) {
  list(
    beta_nb        = m$coef$Beta.NB,
    beta_zc        = m$coef$Beta.multinom.ZC,   # NULL for evinb
    beta_pl_mult   = m$coef$Beta.multinom.PL,
    beta_pl        = m$coef$Beta.PL,
    alpha_nb       = m$coef$Alpha.NB,
    c_ev           = m$coef$C,
    log_lik        = m$log.lik,
    aic            = m$AIC,
    bic            = m$BIC,
    props          = m$props,
    resp           = m$resp,
    fitted         = m$fitted,
    converge       = m$converge,
    n_above_c      = m$n_above_c,
    loglik_recomputed = isTRUE(m$loglik_recomputed),
    c_profile      = m$c_profile
  )
}

ctrl <- evinf_control(c.lim = c(50, 1000), init.C = 200)

fit_z <- function(formula, data, ...) {
  suppressMessages(suppressWarnings(
    evzinb(formula, data = data, control = ctrl, verbose = FALSE, ...)
  ))
}
fit_i <- function(formula, data, ...) {
  suppressMessages(suppressWarnings(
    evinb(formula, data = data, control = ctrl, verbose = FALSE, ...)
  ))
}

data(genevzinb2, package = "evinf")
set.seed(42)
gf <- genevzinb2
gf$g <- factor(sample(c("a", "b", "c"), nrow(gf), replace = TRUE))

## 1-2: plain fits, bootstrap = FALSE
saveRDS(capture_fit(fit_z(y ~ x1 + x2 + x3, genevzinb2, bootstrap = FALSE)),
        file.path(fixtures_dir, "em_baseline_evzinb.rds"))
saveRDS(capture_fit(fit_i(y ~ x1 + x2 + x3, genevzinb2, bootstrap = FALSE)),
        file.path(fixtures_dir, "em_baseline_evinb.rds"))

## 3: factor + in-formula transform
f3 <- y ~ x1 + log(abs(x2) + 1) + g
saveRDS(capture_fit(fit_z(f3, gf, bootstrap = FALSE)),
        file.path(fixtures_dir, "em_baseline_evzinb_factor.rds"))
saveRDS(capture_fit(fit_i(f3, gf, bootstrap = FALSE)),
        file.path(fixtures_dir, "em_baseline_evinb_factor.rds"))

## 4: bootstrap coefficient matrices (n = 3, fixed seed)
## NB: regenerated in round 5 Part A - the parallel backend moved from
## foreach/%dorng% to future/furrr, so the per-bootstrap resamples (and hence
## these coefficient matrices) changed. The sampling distribution is unchanged.
mzb <- fit_z(y ~ x1 + x2 + x3, genevzinb2, bootstrap = TRUE,
             n_bootstraps = 3, boot_seed = 123, multicore = FALSE)
mib <- fit_i(y ~ x1 + x2 + x3, genevzinb2, bootstrap = TRUE,
             n_bootstraps = 3, boot_seed = 123, multicore = FALSE)
saveRDS(list(
  evzinb = suppressWarnings(coefficient_extractor(mzb, "all")),
  evinb  = suppressWarnings(coefficient_extractor(mib, "all"))
), file.path(fixtures_dir, "em_baseline_boot.rds"))

## 5: hks, ISQ replication formulas (bigger; test is skip_on_cran)
data(hks, package = "evinf")
f_hks <- osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
  lntpop + brv_AllLag_log + osvAllLagDum + incomp
f_hks_pareto <- ~ troopLag_log + epdur_log + brv_AllLag_log + lntpop
saveRDS(capture_fit(fit_z(f_hks, hks, formula_pareto = f_hks_pareto,
                          bootstrap = FALSE)),
        file.path(fixtures_dir, "em_baseline_evzinb_hks.rds"))
saveRDS(capture_fit(fit_i(f_hks, hks, formula_pareto = f_hks_pareto,
                          bootstrap = FALSE)),
        file.path(fixtures_dir, "em_baseline_evinb_hks.rds"))

message("baselines written to ", fixtures_dir)
