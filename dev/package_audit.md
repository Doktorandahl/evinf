# evinf package audit and development suggestions

Prepared 2026-09-08 alongside the first draft of the JSS vignette (`vignettes/evinf.qmd`).
Based on a read of every file in `R/`, `DESCRIPTION`, `NAMESPACE`, `inst/`, `man/` (spot checks), and `src/` at package version 0.9.3 (last commit 2025-11-13). Nothing has been executed; R is not available on the machine this was run from, so every item below is from reading the code and should be confirmed with a run before fixing.

Items are grouped by severity. Within each group, the ones most relevant to the vignette come first.

---

## 1. Substantive statistical issues (fix before the vignette is finalized)

### 1.1 `lr_test()` reports half the LR statistic

`R/lr_tests.R`:

```r
dplyr::mutate(statistic = .data$loglik_full - .data$loglik_restricted,
              prob = 1 - pchisq(.data$statistic, .data$df))
```

`object$log.lik` is the plain log likelihood (the same file computes `AIC <- 2*n.par - 2*func.val` from it in `a019.R`), so the LR statistic should be `2 * (loglik_full - loglik_restricted)`. As written, the test is conservative by a factor of two, and the same un-doubled statistic feeds the bootstrap branch (`ks.test(..., pchisq, df)`, `chisq_mean`, `chisq_median`, `prop_sig`), where the chi-square reference is then wrong as well.

Note that Table B2 in the ISQ online appendix reproduces exactly this: `ll_restricted = -5590.83, ll_full = -5574.38, statistic = 16.44`. The p-values in that table would all still round to 0.00 after doubling, so the published conclusions are unaffected, but the vignette (Section 2.3, eq. LR test) should describe the correct statistic and the code should match it.

### 1.2 `compare_models()` fits the bootstrapped ZINB with the wrong formula

`R/zinb_comparison.R`: the full-sample ZINB uses `f_zinb` (`y ~ nb_terms | zi_terms`), but `inner_zinb()` calls `pscl::zeroinfl(formulas$formula_zi, ...)`, i.e. the zero-inflation formula *without* the `|` part. `zeroinfl()` then uses the ZI regressors for both the count and zero parts. Whenever `formula_nb != formula_zi`, the bootstrapped ZINB coefficients, OOB predictions, and fit statistics are for a different model than the full-sample fit. Pass `f_zinb` through to `inner_zinb()`.

### 1.3 `summary.evzinb()` / `summary.evinb()`: `p_value = "both"` never reaches the approximate branch

```r
if (p_value %in% c('bootstrapped','both')) { ... } else if (p_value %in% c('approx','both')) { ... }
```

`"both"` is consumed by the first condition. Split into two independent `if` blocks.

### 1.4 `summary()` errors on non-bootstrapped models

`n_failed_bootstraps`, `nb_boot`, etc. are only defined inside `if (!is.null(object$bootstraps))`, but are used unconditionally afterwards (`res$n_failed_bootstraps`, the `standard_error` block). A model fitted with `bootstrap = FALSE` cannot be summarized. Same pattern in `tidy.*` when `standard_error = TRUE` (the default). Either guard these or set `standard_error`/`p_value` defaults conditionally.

### 1.5 `lr_test()` passes the wrong data object for `evinb` models

```r
evinb(..., data = object$data, ...)   # should be object$data$data
```

The evzinb branch correctly uses `object$data$data`. `object$data` is a list, so the evinb branch will fail.

### 1.6 `compare_models()` does not work for `evinb` objects

`iv_zi <- all.vars(object$formulas$formula_zi)[-1]` gives `character(0)` (not `NULL`) when `formula_zi` is `NULL`, so the `is.null(iv_zi)` guard never fires and `f_zinb` becomes `y ~ x1 + x2 |`, which `as.formula()` rejects. The result object is also always named `out$evzinb`. Either support evinb explicitly (`zinb_comparison = FALSE` default for evinb, name the slot `model`) or error early with a clear message.

### 1.7 Design matrices are built from `model.frame()[, -1]`, not `model.matrix()`

`run_evzinb()`/`run_evinb()` (and `prob_from_*`, `counts_from_evzinb`, `fitted_alpha_from_evzinb` for `newdata`) use `as.matrix(mf[, -1])`. Consequences:

* factor/character covariates produce a character matrix and the C++ code will fail or silently misbehave;
* interactions (`x1:x2`, `x1*x2`) and `I()`/`poly()` terms are not expanded;
* in-formula transformations such as `log(x)` *do* work numerically (model.frame evaluates them) but coefficient names come from `all.vars(formula)[-1]`, so the coefficient is labelled `x`, not `log(x)`; and for `y ~ x1 + log(x1)` two coefficients get the same name `x1`;
* `offset()` is impossible (see 3.9 below).

Switching to `model.matrix()` (dropping the intercept column, since the C++ side adds it) and naming coefficients from `colnames(model.matrix(...))` fixes all of this in one place. Until then, the vignette states that transformations must be pre-computed as columns (Section 5, "Prepare covariates before fitting").

### 1.8 Length check for `init.Beta.NB` uses the wrong model frame

Both `run_evzinb()` and `run_evinb()`:

```r
if (is.null(init.Beta.NB) | length(init.Beta.NB) != ncol(mf_evi)) {   # should be ncol(mf_nb)
```

User-supplied NB starting values are silently discarded whenever the NB and EVI formulas differ in length.

---

## 2. Bugs and inconsistencies that affect usability but not published results

2.1 `evinb()`: `names(boots) <- paste('bootstrap_', ...)` is executed *after* `out <- c(full_run, list(bootstraps = boots))`, so the bootstraps in the returned object are unnamed (evzinb does it in the right order).

2.2 `mr_inner()` in `zinb_comparison.R` computes `obj[!(names(obj) %in% c(...))]` but never assigns it, so the intended slimming of NB/ZINB bootstrap objects is a no-op (memory only; harmless otherwise).

2.3 `prob_from_evzinb()`: `if (min(pr_count < -1e20))` is `min()` of a logical vector and is always 0/FALSE; presumably intended `min(pr_count) < -1e-20` (or similar). The negative-probability guard never triggers.

2.4 `print.evinb()` prints "Fitted EVZINB model with formulas:". `print.evzinb()` prints the NB right-hand side in the ZI line (`paste(x$formulas$formula_nb)[3]` on the `ZI:` row).

2.5 `predict.nbboot()` exists but has no `@export`/`S3method(predict, nbboot)` in NAMESPACE and no Rd, so `predict(comp$nb)` dispatches to `predict.glm` on a list and fails. `predict.zinbboot` is exported and documented.

2.6 `summary()` returns an object of class `summary.evzinb` / `summary.evinb` but there is no `print.summary.*` method, so `summary(model)` prints a raw nested list. This is the first thing a new user sees.

2.7 Component vocabulary differs across functions and will confuse users (and complicates the vignette):

| function | count | zero | EV inflation | Pareto |
|---|---|---|---|---|
| `tidy()` | `"count"` | `"zi"` | `"evi"` | `"pareto"` |
| `coefficient_extractor()` | `"nb"` | `"zi"` | `"evinf"` | `"pareto"` |
| `predict()` | `"counts"` / `"count_state"` | `"zi"` | `"evinf"` | `"pareto_alpha"` |
| `summary()` list names | `negative_binomial` | `zero_inflation` | `extreme_value_inflation` | `pareto` |
| `props`/`resp` colnames | `count` | `zero` | – | `pareto` |
| `predict(type="states")` colnames | `pr_count` | `pr_zc` | – | `pr_pareto` |

Pick one set (I'd suggest `"count"`, `"zero"`, `"evi"`, `"pareto"`) and accept the others as aliases for a deprecation cycle.

2.8 `tidy()` default is `component = "zi"`, so a bare `tidy(model)` (and therefore a bare `modelsummary(model)`) shows only the zero-inflation component. `"all"` is the more natural default.

2.9 `tidy.evzinb()`/`tidy.evinb()` accept and validate `confint` and `conf_level` but never use them: no `conf.low`/`conf.high` columns are ever produced (the arguments appear only in the signature and the `match.arg()` call). Either implement (bootstrap percentile intervals are trivial from `*_boot`) or drop the arguments. The vignette currently says `tidy()` "optionally" returns confidence intervals, following the documentation; that sentence should be revisited once this is settled.

2.10 `glance.evinb()` rounds `alpha` and `logLik`; `glance.evzinb()` does not. Neither includes a `converged` flag or `n_bootstraps`/`n_failed_bootstraps`, which would be useful GOF rows.

2.11 `predict.evzinb(multicore = TRUE)` registers a parallel backend and never stops it; with `multicore = FALSE` the `%dopar%` calls emit "executing %dopar% sequentially: no parallel backend registered" warnings unless a backend happens to be registered. Same in `lr_test(bootstrap = TRUE)` (evzinb branch registers, evinb branch does not; neither stops).

2.12 `run_evinb()` has `c.lim = c(70, 300)` as default while `evinb()` has `c(50, 1000)`; harmless because `evinb()` always passes it, but confusing. `run_evinb()` uses `match.arg()` on `pdf.pl.type`; `evinb()`/`evzinb()` do not validate it.

2.13 `a019.R`, end phase: `func.val <- max(log.lik.vec.all)` while warm-up uses `max(log.lik.vec)` (the current c-profile). The reported `log.lik`, AIC and BIC are therefore the best value seen anywhere in the trace, not necessarily the log likelihood at the returned parameters. With the monotone line search these should coincide, but it is worth asserting.

2.14 `summary.*`: `dplyr::left_join(props, props_boot)` and several other joins have no `by =`, producing "Joining with `by = ...`" messages at every call.

2.15 `revzinb_fit()` draws Pareto values one observation at a time in a `foreach` loop (`mistr::rpareto(1, C_est, alphs[j])`); `mistr::rpareto(n, C, alphs)` is vectorised and would be orders of magnitude faster.

---

## 3. Documentation and metadata

3.1 `inst/CITATION` and `inst/REFERENCES.bib`: the ISQ paper is cited as `year = 2023, volume = "x", number = "x", pages = "x"`. It is now *International Studies Quarterly* 68(4), sqae137, 2024, doi:10.1093/isq/sqae137. The `hks` Rd uses key `randahl2023evzinb`.

3.2 `man/hks.Rd`: `policeLag` is described as "thousands of troops" and `militaryobserversLag` as "The number of UN military troops" — should be police and military observers respectively.

3.3 `man/genevzinb.Rd`, `man/genevzinb2.Rd`: "one dependent and three dependent variables" → independent; title says "EVZBINB".

3.4 `prune.c.range` is documented as `length(c.lim) * prune.c.range` in `evzinb()`/`run_evinb()` and as `length(c.lim) * (1 - prune.c.range)` in `run_evzinb()`/`evinb()`. The code keeps `ceiling(length(c.range) * (1 - prune.c.range))` values, so the `(1 - p)` version is correct. Also, `c.lim` is documented as "Integer range defining the possible values of C" — clearer: "numeric vector of length 2; the candidate set for C is the unique observed values of the response within this range".

3.5 The `man/` files are stale relative to `R/`: e.g. `man/evzinb.Rd` and `man/evinb.Rd` lack `prune.c.range`. Re-run `devtools::document()` (RoxygenNote is 7.2.3; current roxygen2 is 7.3.x).

3.6 `NEWS.md` stops at 0.8.0 ("Initial CRAN submission"); DESCRIPTION is 0.9.3; CRAN-SUBMISSION says 0.8.10. `README.md` is empty (8 bytes). The pkgdown `docs/` are from Aug 2023.

3.7 `print.evzinbcomp()` roxygen: "@param x A fitted evinb model", "@return An evinb print function", and the example prints an evinb model rather than a comparison object. `coefficient_extractor.zinbboot` says "A fitted evinb model".

3.8 Several roxygen examples use `multicore = TRUE, ncores = 2` — CRAN policy allows at most 2 cores in examples, so this is fine, but the examples fit 10 bootstraps of a 4-equation model and will be slow on CRAN check machines; consider `\donttest{}` and use `hks` (real data) in at least one example.

3.9 DESCRIPTION:
* `Depends: generics` should be `Imports: generics` (the package only needs the generics for `tidy`/`glance`; re-export them so users don't need `library(generics)`/`library(broom)`).
* `MLmetrics` is imported for four one-line metrics; `stringi` *and* `stringr` are both imported for one call each (`stri_detect_*` in `lr_tests.R`, `str_split` in `zinb_comparison.R`), both of which base R handles. `methods`, `Rdpack` fine.
* Missing `Suggests: modelsummary, ggplot2, knitr, rmarkdown/quarto, testthat`; add `VignetteBuilder` once the vignette is in place.
* Maintainer email is `david.randahl@pcr.uu.se`; the pintervals manuscript uses `david.randahl@fhs.se`.
* `Description` contains "'tidy'" as if it were a package.

3.10 Repository hygiene: an 18 MB `.RData` and a `.Rhistory` sit in the package root; `src/*.o`, `src/evinf.so`, `.DS_Store` are present (gitignored, but `.Rbuildignore` does not exclude `.RData`/`.Rhistory`, which R CMD build will complain about). `src/a019_old.cpp` is the *only* C++ source and should be renamed (`evinf.cpp`); `src/coef_extractor.cpp` is a two-line stub with no functions. No `tests/` directory.

---

## 4. Suggested additional functionality

Roughly in order of value for usability, and with an eye to what a JSS reviewer would expect.

### 4.1 Standard S3 generics
`coef()`, `vcov()` (bootstrap covariance), `confint()`, `logLik()`, `AIC()`, `BIC()`, `nobs()`, `fitted()`, `residuals()` (e.g. Pearson-type on the harmonic prediction, or randomized quantile residuals from the mixture CDF already in `mixture_p()`), `simulate()` (wrapping `revzinb_fit`), `update()`, `formula()`, `terms()`, `model.frame()`. Besides being expected by users, `coef`/`vcov`/`nobs` are what `modelsummary`, `marginaleffects`, `texreg`, `sandwich`, etc. dispatch on, so this is the cheapest way to plug into that ecosystem.

### 4.2 `print.summary.evzinb()` and a better `print.evzinb()`
A formatted summary with one block per component (estimate, boot SE, boot p, stars), followed by alpha_NB, C_EV, state proportions, n obs / n boot / n failed, log-lik, AIC, BIC. `print.evzinb` could also show C_EV and the number of observations above it, which is the single most useful diagnostic.

### 4.3 A `control` object
`evzinb()` currently has ~30 arguments. An `evinf_control()` constructor (like `glm.control()`/`zeroinfl.control()`) holding `max.diff.par`, `max.no.em.steps*`, `c.lim`, `prune.c.range`, `max.upd.par.*`, `no.m.bfgs.steps.*`, `pdf.pl.type`, `eta.int`, and the `init.*` values would make the fitting signature `evzinb(formula_nb, formula_zi, formula_evi, formula_pareto, data, bootstrap, n_bootstraps, block, boot_seed, multicore, ncores, control, verbose)` and easier to document. Alternatively adopt a single two-sided formula interface `y ~ nb | zi | evi | pareto`, mirroring `pscl::zeroinfl`'s `y ~ count | zero`.

### 4.4 `model.matrix()`-based design (see 1.7)
Enables factors, interactions, splines, in-formula transforms, and `offset()`. An **exposure offset** (e.g. `offset(log(population))`) in the NB and/or Pareto component is highly relevant for conflict data and is currently impossible.

### 4.5 Threshold diagnostics
* `plot_c_profile()` / `c_profile()`: log-likelihood profile over the candidate `c` values from the final EM iteration (already computed in `log.lik.vec` inside `zerinfl.nb.pl.regression.fun`; just needs to be stored). This is the natural way to justify `c.lim` and to show whether C_EV is well identified.
* Store `c_trace` and the log-likelihood trace in the object; add `converged`, `n_em_steps`, `n_above_c` to `glance()`.
* A data-driven default for `c.lim` (e.g. from quantiles of the positive values) with a message, instead of the hard-coded `c(50, 1000)`.

### 4.6 Plot methods
`plot(model, type = ...)` for the three figures every user will want: (i) state probabilities over a covariate (Figure C1 in the appendix), (ii) predicted harmonic mean and chosen quantiles over a covariate with bootstrap ribbons (Figure 1 in the paper), (iii) bootstrap coefficient densities per component (Figure B2). A `predict_grid()`/`plot_effects()` helper that takes a variable name and holds the others at means/medians/modes would remove most of the boilerplate in the vignette's Section 4.5. `ggplot2` in Suggests.

### 4.7 Comparison summaries
`compare_models()` returns the raw objects. Add `compare_fit(comp)` returning a tidy tibble of paired bootstrap differences in AIC, BIC, RMSE, RMSLE with median difference and "proportion EVZINB better" (Table B3), and a corresponding density plot (Figure B3). `oob_evaluation()` could take the `evzinbcomp` object and return one column per model directly. Also let `compare_models()` accept evinb objects (1.6).

### 4.8 Marginal effects
`a019.R` contains unused `marginal.effect.*.fun()` code. A user-facing `marginal_effects()`/`avg_slopes()`-style function computing bootstrap-based average marginal effects of a covariate on (a) each state probability, (b) the harmonic-mean prediction, and (c) selected quantiles would deliver on the paper's closing remark about "a unified framework for analysis of the effect of covariates across states". Implementing `coef`/`vcov`/`predict(newdata=)` cleanly may already get most of this for free via `marginaleffects`.

### 4.9 Bootstrap ergonomics
* `add_bootstraps(model, n)` to extend an existing fit rather than refitting everything when 100 turns out to be too few.
* Progress via `progressr` and parallelism via `future`/`furrr` (uniform on Windows/macOS/Linux, respects the user's own `plan()`), replacing the OS-branching `doParallel` code that is duplicated in `evzinb()`, `evinb()`, `compare_models()`, `lr_test()`, `predict()`.
* Store the bootstrap seed and the `boot_id` indices (already done) so that `compare_models()` and `lr_test(bootstrap = TRUE)` remain reproducible; document this.
* A `failed_bootstraps(model)` accessor returning the error messages, which is useful diagnostic information the package currently discards with `try()`.

### 4.10 Posterior-state tools
`classify_states(model, rule = "map")` returning the MAP state per observation plus the posterior probabilities, and a `state_table()` cross-tabulating prior vs posterior classification. The vignette currently reaches into `model$fitted$posterior_pareto` for this.

### 4.11 Covariate-dependent C_EV
Named in the paper's discussion as the natural next extension. Even a two-stage approximation (`c_EV = exp(x' delta)` with the grid search done per stratum of a single categorical covariate) would be a publishable addition and could be flagged as experimental.

### 4.12 Testing and infrastructure
A `testthat` suite using `genevzinb2` (fast) and `hks` (realistic): recover simulation parameters within tolerance, `predict(newdata = data) == fitted`, `tidy()` shapes, `lr_test` statistic equals `2 * diff(logLik)`, `compare_models` ZINB formula, evinb/evzinb equivalence when the ZI intercept is fixed at -Inf, and snapshot tests for `print`. GitHub Actions for R CMD check on the three OSs; rebuild pkgdown; add the vignette (requires `VignetteBuilder: quarto` and `Suggests: quarto`, or an `.Rmd` with `rticles::jss_article`).

### 4.13 Smaller conveniences
* `hks` data has variables `lntroopLag`, `lnepdur` pre-computed for the Pareto component; document why (ties in with 1.7).
* `gm_evzinb` should include rows for `parameter` (C_EV) with a label like "C_EV", `alpha` as "alpha_NB", and be exported as a function `gof_map_evinf()` so it can be extended with `n_failed_bootstraps`.
* Vectorise the `for (i in 1:n)` fitted-value loop at the end of `zerinfl.nb.pl.regression.fun()`.
* Accept `block` as a bare column name (tidy evaluation) as well as a string.
