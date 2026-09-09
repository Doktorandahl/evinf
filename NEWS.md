# evinf 0.10.0

Implements section 4 of the internal package audit (new functionality) and the
review follow-ups in `dev/review_round1.md`.

## New features

* `evinf_control()` bundles the EM tuning settings; `evzinb()` / `evinb()` gain
  a `control` argument. The individual tuning arguments still work but are
  deprecated in favour of `control`.
* The candidate range for C_EV is now chosen from the data by default
  (`c.lim = NULL`); it is printed with a message. `c_profile()` /
  `plot_c_profile()` show the log-likelihood profile over that range;
  `glance()` gains `n_above_c` and `n_em_steps`; the fitted object carries
  `$c_profile`, `$c_trace`, `$loglik_trace`.
* `offset()` in the count-component formula is supported
  (`mu_NB = exp(x'b + offset)`), e.g. `y ~ x + offset(log(exposure))`.
* Standard S3 methods: `coef()`, `vcov()` (bootstrap covariance), `confint()`,
  `logLik()` / `AIC()` / `BIC()`, `nobs()`, `formula()`, `terms()`,
  `model.frame()`, `fitted()`, `residuals()` (response and randomized quantile),
  `simulate()`, `update()`. `revzinb_fit()` / `revinb_fit()` are superseded by
  `simulate()`.
* `add_bootstraps()` extends an existing fit with more replicates;
  `failed_bootstraps()` returns the error messages of the replicates that
  failed. Bootstrap seeds are recorded in `$boot_seeds`.
* Progress reporting for the bootstrap loops (and for `lr_test(bootstrap =
  TRUE)` and `compare_models()`) through `progressr` (shown when `verbose =
  TRUE` or a `progressr` handler is set), replacing the old `cat()` messages.
* `classify_states()` returns the prior and posterior state classification per
  observation; `state_table()` cross-tabulates them.
* `predict_grid()` builds a one-variable prediction grid (other covariates held
  at their mean/median or modal level, or pinned via `at`) and evaluates the
  model over it in long form, optionally with bootstrap confidence intervals.
* `plot()` methods for fitted models: `type = "states"` (prior state
  probabilities over a covariate), `"prediction"` (harmonic-mean prediction with
  a bootstrap ribbon and quantile lines), `"coefficients"` (bootstrap
  coefficient densities), `"ppc"` (observed-vs-expected binned frequencies) and
  `"ppc_quantiles"` (observed vs. simulated tail quantiles with a 5-95\% band).
* `compare_fit()` summarises the paired bootstrap differences (`compared -
  evinf`) in AIC, BIC and out-of-bag RMSE / RMSLE for an `evzinbcomp` object,
  with the proportion of bootstraps favouring the extreme-value model;
  `plot()`, `tidy()` and `glance()` methods for `evzinbcomp` and an
  `oob_evaluation()` method that tabulates the out-of-bag error per model.
* `marginal_effects()` computes average marginal effects (central difference for
  numeric covariates, level-vs-reference contrasts for factors) on the harmonic
  mean, the state probabilities or a predicted quantile, with bootstrap
  confidence intervals.
* `marginaleffects` compatibility: `evzinb` / `evinb` models register with
  `marginaleffects` on load, so `marginaleffects::avg_slopes()` and friends work
  with bootstrap delta-method standard errors.
* `gof_map_evinf()` returns the `modelsummary` goodness-of-fit map (with an
  `extra` argument to append rows); the bundled `gm_evzinb` data object is now
  generated from it and gains an `obs_above_c_ev` row.
* `evzinb()` / `evinb()` accept `block` as a bare column name
  (`block = id`) as well as a string.
* `print(summary(model))` reports the number of observations at or above C_EV
  and notes when the log-likelihood was recomputed.

## Bug fixes

* `marginal_effects(type = "quantile")` no longer returns zeros. The mixture
  quantile is rounded to an integer for `predict()`, but a central finite
  difference of that step function is zero almost everywhere; the marginal-
  effects path now differentiates the continuous (unrounded) quantile.
  `quantiles_from_*()` gain a `round` argument (`TRUE` by default;
  `predict(type = "quantile")` is unchanged).
* `marginal_effects(method = "difference")` is now implemented for numeric
  covariates (average change from a `delta`-unit increase, default 1 unit); it
  was previously accepted and ignored. `type = "quantile"` defaults to
  `method = "difference"`.
* `marginal_effects(type = "quantile")` caps the rows it averages over at
  `n_max` (new argument, default 500; automatic above `nrow * n_bootstraps =
  2e5`, or whenever set explicitly; `n_max = Inf` disables it), with a message.
  The per-observation `mistr` quantile machinery made this prohibitively slow
  on data the size of `hks`.
* `evzinb()` / `evinb()` resolve a bare `block` name to the matching column of
  `data` before trying to evaluate it, so an unrelated object of the same name
  in the calling environment no longer shadows the column.
* `add_bootstraps()` errors if `boot_seed` was already used for this model (the
  `%dorng%` stream is fully determined by the seed, so it would silently
  duplicate existing draws and shrink the bootstrap SEs). `evzinb()` / `evinb()`
  now always record the seed actually used (drawing one when `boot_seed = NULL`),
  so `object$boot_seeds` is complete.
* `mixture_p()` (the mixture CDF behind `residuals(type = "quantile")`) now
  uses the discretised Pareto CDF that matches the pmf the likelihood uses
  (`1 - (C/(y+1))^a` for `y >= C`), rather than the continuous `1 - (C/y)^a`.
  Randomized quantile residuals shift by a small amount in the extreme-value
  tail (on `genevzinb2`, 12 of 100 residuals move, max change 0.015).
* `evinb()` no longer runs the EVZINB EM step (which updates the zero-inflation
  block) during the warm-up phase -- a copy-paste artefact from the pre-0.10.0
  code. The zero-inflation intercept is now held fixed throughout, as the model
  intends. Estimates are unchanged to machine precision on the tested data (the
  fixed intercept starts far enough from zero that the spurious update was
  numerically inert).
* Block-bootstrap: the block variable is now included in the single `na.omit()`
  step, so a model whose block column has missing values where the model
  variables do not no longer mis-indexes the resamples (previously could error
  or silently sample the wrong rows).
* `lr_test()` refits the restricted models with the full model's control
  settings (candidate range for C, tolerances, `pdf.pl.type`, ...) and warm-
  started from the full-model estimates. Previously it used the package
  defaults, so a model fitted with a non-default `c.lim` was compared against a
  restricted fit that searched a different candidate set for C_EV (the LR
  statistic could even come out negative).
* The internal check that the reported log-likelihood equals the log-likelihood
  at the returned parameters no longer warns from inside every bootstrap; it
  recomputes silently, records `object$loglik_recomputed`, and (only for the
  full-sample fit, only when `verbose = TRUE`) emits a single message.
* A user-registered parallel backend (`doParallel`, `doFuture`, ...) is left
  untouched when a function is called with `multicore = FALSE`.
* A non-finite value produced by an in-formula transformation (e.g. `log(x)`
  with `x <= 0`) now raises an informative error naming the column instead of
  misbehaving in the C++ code.
* User-supplied `init.Beta.NB` is used (regression-test hardened).
* `evinb` parameter count and `predict.*boot(pred = "original")` indexing fixes
  (carried over from the 0.9.4 audit follow-ups).
* `predict()` on the `zinb` slot of a `compare_models()` result with
  `type = "counts"` (or `type = "all"`, or `type = "counts", confint = TRUE`)
  no longer errors with "$ operator is invalid for atomic vectors".

## Breaking changes / deprecations

* The minimum R version is now 4.1.0. The `NAMESPACE` uses delayed S3 method
  registration for `insight` / `marginaleffects`, which needs R >= 3.6.0; 4.1
  is the common floor for packages using that pattern.
* The approximate bootstrap runtime estimate is printed only when
  `verbose = TRUE` (it was always printed before 0.9.4; this note was missing
  from the 0.9.4 changelog).
* `predict(type = "all")` and every `confint = TRUE` output now lead with the
  canonical state-probability columns (`pr_zero`, `pr_count`, `pr_evi`), keeping
  the deprecated `pr_zc` / `pr_pareto` as trailing duplicates — matching what
  `predict(type = "states")` already did.
* Default fits change slightly: with `c.lim = NULL` the candidate set for C_EV
  is data-driven and `init.C` defaults to the median of that set, so estimates
  from a default 0.9.x call are not bit-identical to a default 0.10.0 call. Pin
  `evinf_control(c.lim = ..., init.C = ...)` to reproduce an older fit.
* `gm_evzinb` gains a row (`obs_above_c_ev`) and its `parameter` row is now
  formatted with 0 decimals; regenerate any cached copy with `gof_map_evinf()`.
* The `hks` column `brv_AllLag` is renamed `brv_AllLag_log` (it is on the log
  scale, like the other `_log` columns). There is no alias — code referring to
  the old name will get an immediate "object not found" / "column not found"
  error.
* `predict()` on the `nb` / `zinb` slot of a `compare_models()` result no longer
  accepts `type = "evinf"` (a plain NB / ZINB has no extreme-value state);
  `match.arg()` now rejects it.
* The three internal `marginal.effect.*` helpers (never exported, superseded by
  `marginal_effects()`) were removed.

## Documentation

* New `pkgdown` reference index grouping the exported functions (Fitting,
  Summaries & tables, Prediction & effects, Diagnostics & plots, Model
  comparison, Simulation, Data).
* `?hks` now documents the columns as they actually ship (the previous help
  page listed several columns that are not in the data and omitted the `_log`
  transforms) and explains that `log1p(troopLag)` etc. can be used directly in
  the formula since 0.9.4.
* GitHub Actions workflows for `R CMD check`, test coverage and the pkgdown
  site; coverage and check badges in the README.
* The per-observation fitted-value loop in the estimation routines is
  vectorised (no change to results).
* Internal: the twelve near-identical bootstrap-refit loops in
  `compare_models()` are replaced by a single `boot_refit_family()` helper
  (output unchanged).
* Internal: the 1,800-line `R/a019.R` is split into six documented,
  de-duplicated files (`R/em_candidates.R`, `R/em_profile_c.R`, `R/em_step.R`,
  `R/em_fixed_c.R`, `R/em_fitted.R`, `R/em_fit.R`), one EM/ECME driver
  `em_fit()` for both models instead of the two ~90%-duplicated
  `*.regression.fun()` pairs. Estimates are unchanged (a numerical-identity
  test suite pins the pre-split output). The three unused
  `prediction.*.fun()` legacy helpers were removed.
* `evinf_control(prune.c.range = ...)` now also thins the \eqn{C_{EV}}
  candidate set for `evinb()` (it was silently ignored there before), and the
  "more than 100 candidates" warning now fires for `evinb()` too.
* `?hks`: the `_log` columns are documented as a rule (`log1p()` for the
  personnel counts, `log()` for the episode duration), and the help notes that
  the bundled data carry no conflict identifier.
* `?marginal_effects` / `?tidy.evzinb` / `?gof_map_evinf` gain guidance on
  harmonic-mean effect intervals and on the `modelsummary` `shape` argument.


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
