# Shared fixtures for the test suite.
#
# Fast models use genevzinb2 (100 rows), few bootstraps, no parallelism, and a
# fixed bootstrap seed so results are reproducible. They pin c.lim = c(50, 1000)
# (the pre-0.10.0 default) so structural tests are unaffected by the data-driven
# c.lim default; tests of that default pass their own control. hks tests are
# marked skip_on_cran() in the test file itself.

# Pin the pre-0.10.0 defaults (fixed c.lim and init.C) for reproducibility.
.fast_control <- function(...) {
  evinf::evinf_control(c.lim = c(50, 1000), init.C = 200, ...)
}

fit_evzinb_fast <- function(formula_nb = y ~ x1 + x2 + x3,
                            ...,
                            control = .fast_control(),
                            bootstrap = TRUE,
                            n_bootstraps = 5) {
  data(genevzinb2, package = "evinf", envir = environment())
  suppressMessages(evinf::evzinb(
    formula_nb,
    data = genevzinb2,
    control = control,
    bootstrap = bootstrap,
    n_bootstraps = n_bootstraps,
    multicore = FALSE,
    boot_seed = 123,
    verbose = FALSE,
    ...
  ))
}

fit_evinb_fast <- function(formula_nb = y ~ x1 + x2 + x3,
                           ...,
                           control = .fast_control(),
                           bootstrap = TRUE,
                           n_bootstraps = 5) {
  data(genevzinb2, package = "evinf", envir = environment())
  suppressMessages(evinf::evinb(
    formula_nb,
    data = genevzinb2,
    control = control,
    bootstrap = bootstrap,
    n_bootstraps = n_bootstraps,
    multicore = FALSE,
    boot_seed = 123,
    verbose = FALSE,
    ...
  ))
}

genevzinb2_factor <- function() {
  data(genevzinb2, package = "evinf", envir = environment())
  set.seed(42)
  genevzinb2$g <- factor(sample(c("a", "b", "c"), nrow(genevzinb2), replace = TRUE))
  genevzinb2
}
