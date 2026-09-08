# evinf 0.9.4

This release works through sections 1-3 of the internal package audit: it fixes
substantive statistical issues, a number of usability bugs, and brings the
documentation and metadata up to date. A `testthat` suite has been added.

## Breaking changes / deprecations

* **Likelihood-ratio statistic.** `lr_test()` now reports the standard
  `2 * (logLik_full - logLik_restricted)` statistic (and uses it for the
  bootstrap reference distribution). Statistics reported by earlier versions were
  half the correct value; p-values change accordingly but the qualitative
  conclusions in Randahl & Vegelius (2024) are unaffected. Degrees of freedom are
  now the number of design-matrix columns removed across components (so removing
  a four-level factor contributes three degrees of freedom).
* **Component vocabulary.** Model components are now named `"count"`, `"zero"`,
  `"evi"` and `"pareto"` everywhere (`tidy()`, `coefficient_extractor()`,
  `summary()`). The old names `"nb"`, `"zi"` and `"evinf"` still work but emit a
  deprecation warning. `predict(type = )` additionally accepts `"zero"` for
  `"zi"` and `"evi"` for `"evinf"`.
* **`tidy()` default** is now `component = "all"` (was `"zi"`), so a bare
  `tidy(model)` / `modelsummary(model)` shows every component. The `y.level`
  values from `tidy(component = "all")` are `zero`, `evi`, `count`, `pareto`.
* **Column renames** (old names kept as duplicates for this release only):
  `predict(type = "states")` columns are `pr_zero` / `pr_count` / `pr_evi`;
  `model$fitted$prob_pareto` / `posterior_pareto` become `prob_evi` /
  `posterior_evi`; the `props` / `resp` matrices are labelled `zero` / `count` /
  `evi`.
* **`compare_models()`** returns the fitted model in `$model`; `$evzinb` remains
  as an alias. It now also accepts `evinb` objects (with `zinb_comparison`
  forced to `FALSE`).
* `generics` moved from `Depends` to `Imports`; `tidy()` and `glance()` are
  re-exported, so `library(broom)` / `library(generics)` is no longer required.
* Dependencies `MLmetrics`, `stringr` and `stringi` removed.
* Minimum R version raised to 3.5 (data files re-serialised).

## Bug fixes

* `compare_models()` fitted the bootstrapped ZINB with the zero-inflation formula
  for *both* model parts; it now uses the full `count | zero` formula, matching
  the full-sample fit.
* `summary()` and `tidy()` no longer error on models fitted with
  `bootstrap = FALSE`; they return the point estimates with a message.
* `summary(p_value = "both")` now returns both the bootstrap and the approximate
  p-values (the approximate branch was previously unreachable).
* `lr_test()` works on `evinb` objects (it passed the wrong data object before).
* Design matrices are built with `model.matrix()` rather than
  `model.frame()[, -1]`, so factor and character covariates, interactions,
  `poly()` / splines and in-formula transformations such as `log(x)` now work and
  are labelled correctly. `na.omit()` is applied once to the union of all
  component variables, keeping the components row-consistent.
* User-supplied `init.Beta.NB` is no longer silently discarded when the NB and
  EVI formulas differ in length.
* `evinb()` bootstraps in the returned object are now named.
* The negative-probability guard in the predictor actually fires now.
* `predict.nbboot()` / `predict.zinbboot()` are registered as S3 methods.
* `summary()` gains a `print` method; `print.evzinb()` / `print.evinb()` show the
  correct component formulas plus `C_EV` and the count of observations at or
  above `C_EV`.
* `glance()` methods return unrounded values and gain `converged`,
  `n_bootstraps` and `n_failed_bootstraps` columns; `gm_evzinb` updated to match.
* Parallel backends are always unregistered on exit, and no
  "executing sequentially" / "Joining with `by`" messages are emitted.
* `run_evinb()`'s default `c.lim` matches `evinb()`; `pdf.pl.type` is validated.
* The fitting routines check that the reported log-likelihood is the
  log-likelihood at the returned parameters and recompute `logLik` / `AIC` /
  `BIC` (with a warning) if not.
* `tidy(confint = )` is implemented (percentile or approximate intervals).
* `revzinb_fit()` / `revinb_fit()` draw Pareto values in a single vectorised call.
* `evinb` objects no longer count the (fixed) zero-inflation coefficients as free
  parameters: `glance()$npar`, `AIC`, `BIC` and the approximate-test degrees of
  freedom are corrected accordingly.
* `predict()` on `nbboot` / `zinbboot` objects with `pred = "original"` used a
  stray loop index and errored; it now predicts from the full-sample model.

## Documentation

* `inst/CITATION` and `inst/REFERENCES.bib` updated to
  *International Studies Quarterly* 68(4): sqae137 (2024),
  doi:10.1093/isq/sqae137.
* `hks`: `policeLag` and `militaryobserversLag` are now correctly described as UN
  police and military observers; all personnel counts are documented as being in
  thousands.
* `genevzinb` / `genevzinb2`: "independent" (not "dependent") covariates;
  "EVZINB" spelling.
* `prune.c.range` documented consistently across all four fitting functions;
  `c.lim` reworded.
* `man/` regenerated with roxygen2 7.3.3; `src/a019_old.cpp` renamed to
  `src/evinf.cpp`; empty `src/coef_extractor.cpp` removed.
* All model-fitting examples wrapped in `\donttest{}` and reduced to
  `n_bootstraps = 5`, `multicore = FALSE`.
* `README` filled in.


# evinf 0.9.3

* Bug fixes.


# evinf 0.9.2

* `predict()` can return the raw bootstrap draws (`return_bootstraps = TRUE`).


# evinf 0.9.1

* Added `prune.c.range` for thinning the candidate set of C when it contains many
  unique values.
* Informative error when every candidate value of C yields an infinite
  log-likelihood.
* Code formatting.


# evinf 0.9.0

* Package renamed from `evzinb` to `evinf`.
* `revzinb_fit()` / `revinb_fit()` for drawing from a fitted model.
* Default parallel-processing behaviour of `predict()` changed.
* Various minor fixes to non-standard evaluation and naming.


# evinf 0.8.0

* Initial CRAN submission.
