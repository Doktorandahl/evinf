# Changelog

## evinf 0.10.0

Implements section 4 of the internal package audit (new functionality)
and the review follow-ups in `dev/review_round1.md`.

### New features

- Bootstrap replicates whose extreme-value tail collapses (smallest
  fitted Pareto shape below `evinf_control(alpha_floor = )`, default
  `0.001`) are flagged rather than silently skewing summaries. They
  carry `$degenerate` / `$degenerate_reason`;
  [`failed_bootstraps()`](../reference/failed_bootstraps.md) lists them
  (new `type` column, `"error"` / `"degenerate"`);
  [`glance()`](https://generics.r-lib.org/reference/glance.html) gains
  `n_degenerate_bootstraps` and `n_bootstraps` now counts only usable
  replicates; [`print()`](https://rdrr.io/r/base/print.html) /
  [`summary()`](https://rdrr.io/r/base/summary.html) report the count.
  Every bootstrap-summary function
  ([`coef()`](https://rdrr.io/r/stats/coef.html) bootstrap extractor,
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`marginal_effects()`](../reference/marginal_effects.md)) excludes
  degenerate replicates by default and takes
  `exclude_degenerate = FALSE` to keep them.
- [`evinf_control()`](../reference/evinf_control.md) bundles the EM
  tuning settings; [`evzinb()`](../reference/evzinb.md) /
  [`evinb()`](../reference/evinb.md) gain a `control` argument. The
  individual tuning arguments still work but are deprecated in favour of
  `control`.
- The candidate range for C_EV is now chosen from the data by default
  (`c.lim = NULL`); it is printed with a message.
  [`c_profile()`](../reference/c_profile.md) /
  [`plot_c_profile()`](../reference/plot_c_profile.md) show the
  log-likelihood profile over that range;
  [`glance()`](https://generics.r-lib.org/reference/glance.html) gains
  `n_above_c` and `n_em_steps`; the fitted object carries `$c_profile`,
  `$c_trace`, `$loglik_trace`.
- [`offset()`](https://rdrr.io/r/stats/offset.html) in the
  count-component formula is supported (`mu_NB = exp(x'b + offset)`),
  e.g. `y ~ x + offset(log(exposure))`.
- Standard S3 methods: [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) (bootstrap covariance),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) /
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) /
  [`BIC()`](https://rdrr.io/r/stats/AIC.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`formula()`](https://rdrr.io/r/stats/formula.html),
  [`terms()`](https://rdrr.io/r/stats/terms.html),
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) (response and
  randomized quantile),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`update()`](https://rdrr.io/r/stats/update.html).
  [`revzinb_fit()`](../reference/revzinb_fit.md) /
  [`revinb_fit()`](../reference/revinb_fit.md) are superseded by
  [`simulate()`](https://rdrr.io/r/stats/simulate.html).
- [`add_bootstraps()`](../reference/add_bootstraps.md) extends an
  existing fit with more replicates;
  [`failed_bootstraps()`](../reference/failed_bootstraps.md) returns the
  error messages of the replicates that failed. Bootstrap seeds are
  recorded in `$boot_seeds`.
- Progress reporting for the bootstrap loops (and for
  `lr_test(bootstrap = TRUE)` and
  [`compare_models()`](../reference/compare_models.md)) through
  `progressr` (shown when `verbose = TRUE` or a `progressr` handler is
  set), replacing the old [`cat()`](https://rdrr.io/r/base/cat.html)
  messages.
- [`classify_states()`](../reference/classify_states.md) returns the
  prior and posterior state classification per observation;
  [`state_table()`](../reference/state_table.md) cross-tabulates them.
- [`predict_grid()`](../reference/predict_grid.md) builds a one-variable
  prediction grid (other covariates held at their mean/median or modal
  level, or pinned via `at`) and evaluates the model over it in long
  form, optionally with bootstrap confidence intervals.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
  fitted models: `type = "states"` (prior state probabilities over a
  covariate), `"prediction"` (harmonic-mean prediction with a bootstrap
  ribbon and quantile lines), `"coefficients"` (bootstrap coefficient
  densities), `"ppc"` (observed-vs-expected binned frequencies) and
  `"ppc_quantiles"` (observed vs. simulated tail quantiles with a 5-95%
  band).
- [`compare_fit()`](../reference/compare_fit.md) summarises the paired
  bootstrap differences (`compared - evinf`) in AIC, BIC and out-of-bag
  RMSE / RMSLE for an `evzinbcomp` object, with the proportion of
  bootstraps favouring the extreme-value model;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
  [`glance()`](https://generics.r-lib.org/reference/glance.html) methods
  for `evzinbcomp` and an
  [`oob_evaluation()`](../reference/oob_evaluation.md) method that
  tabulates the out-of-bag error per model.
- [`marginal_effects()`](../reference/marginal_effects.md) computes
  average marginal effects (central difference for numeric covariates,
  level-vs-reference contrasts for factors) on the harmonic mean, the
  state probabilities or a predicted quantile, with bootstrap confidence
  intervals. For `type = "quantile"`, whose bootstrap distribution is
  heavy-tailed, the reported `std.error` is a robust scale from the
  percentile interval rather than the raw standard deviation.
- `marginaleffects` compatibility: `evzinb` / `evinb` models register
  with `marginaleffects` on load, so
  [`marginaleffects::avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html)
  and friends work with bootstrap delta-method standard errors.
- [`gof_map_evinf()`](../reference/gof_map_evinf.md) returns the
  `modelsummary` goodness-of-fit map (with an `extra` argument to append
  rows); the bundled `gm_evzinb` data object is now generated from it
  and gains an `obs_above_c_ev` row.
- [`evzinb()`](../reference/evzinb.md) /
  [`evinb()`](../reference/evinb.md) accept `block` as a bare column
  name (`block = id`) as well as a string.
- `print(summary(model))` reports the number of observations at or above
  C_EV and notes when the log-likelihood was recomputed.

### Bug fixes

- `predict(type = "quantile")` and `marginal_effects(type = "quantile")`
  parametrised the negative-binomial part of the mixture by `Alpha.NB`
  where every other part of the package (the likelihood,
  `residuals(type = "quantile")`,
  [`simulate()`](https://rdrr.io/r/stats/simulate.html)) uses
  `size = 1 / Alpha.NB`. Predicted quantiles now use the correct
  dispersion; the effect is small for `evzinb` (`Alpha.NB` near 1) and
  can be large for an over-dispersed `evinb` fit.
- `marginal_effects(type = "quantile")` no longer returns zeros. The
  mixture quantile is rounded to an integer for
  [`predict()`](https://rdrr.io/r/stats/predict.html), but a central
  finite difference of that step function is zero almost everywhere; the
  marginal- effects path now differentiates the continuous (unrounded)
  quantile. `quantiles_from_*()` gain a `round` argument (`TRUE` by
  default; `predict(type = "quantile")` is unchanged).
- `marginal_effects(method = "difference")` is now implemented for
  numeric covariates (average change from a `delta`-unit increase,
  default 1 unit); it was previously accepted and ignored.
  `type = "quantile"` defaults to `method = "difference"`.
- `marginal_effects(type = "quantile")` caps the rows it averages over
  at `n_max` (new argument, default 500; automatic above
  `nrow * n_bootstraps = 2e5`, or whenever set explicitly; `n_max = Inf`
  disables it), with a message. The per-observation mixture-quantile
  solve made this prohibitively slow on data the size of `hks`.
- [`update()`](https://rdrr.io/r/stats/update.html) no longer embeds the
  whole data frame in the refitted model’s stored call.
  [`evzinb()`](../reference/evzinb.md) /
  [`evinb()`](../reference/evinb.md) keep `data` as the expression the
  user passed; [`update()`](https://rdrr.io/r/stats/update.html)
  re-evaluates it, falling back to the embedded copy (bound to a symbol)
  only when the expression can no longer be resolved. `model$call` and
  `str(model)` stay small on large data.
- [`evzinb()`](../reference/evzinb.md) /
  [`evinb()`](../reference/evinb.md) resolve a bare `block` name to the
  matching column of `data` before trying to evaluate it, so an
  unrelated object of the same name in the calling environment no longer
  shadows the column.
- [`add_bootstraps()`](../reference/add_bootstraps.md) errors if
  `boot_seed` was already used for this model (the `%dorng%` stream is
  fully determined by the seed, so it would silently duplicate existing
  draws and shrink the bootstrap SEs).
  [`evzinb()`](../reference/evzinb.md) /
  [`evinb()`](../reference/evinb.md) now always record the seed actually
  used (drawing one when `boot_seed = NULL`), so `object$boot_seeds` is
  complete.
- `mixture_p()` (the mixture CDF behind `residuals(type = "quantile")`)
  now uses the discretised Pareto CDF that matches the pmf the
  likelihood uses (`1 - (C/(y+1))^a` for `y >= C`), rather than the
  continuous `1 - (C/y)^a`. Randomized quantile residuals shift by a
  small amount in the extreme-value tail (on `genevzinb2`, 12 of 100
  residuals move, max change 0.015).
- [`evinb()`](../reference/evinb.md) no longer runs the EVZINB EM step
  (which updates the zero-inflation block) during the warm-up phase – a
  copy-paste artefact from the pre-0.10.0 code. The zero-inflation
  intercept is now held fixed throughout, as the model intends.
  Estimates are unchanged to machine precision on the tested data (the
  fixed intercept starts far enough from zero that the spurious update
  was numerically inert).
- Block-bootstrap: the block variable is now included in the single
  [`na.omit()`](https://rdrr.io/r/stats/na.fail.html) step, so a model
  whose block column has missing values where the model variables do not
  no longer mis-indexes the resamples (previously could error or
  silently sample the wrong rows).
- [`lr_test()`](../reference/lr_test.md) refits the restricted models
  with the full model’s control settings (candidate range for C,
  tolerances, `pdf.pl.type`, …) and warm- started from the full-model
  estimates. Previously it used the package defaults, so a model fitted
  with a non-default `c.lim` was compared against a restricted fit that
  searched a different candidate set for C_EV (the LR statistic could
  even come out negative).
- The internal check that the reported log-likelihood equals the
  log-likelihood at the returned parameters no longer warns from inside
  every bootstrap; it recomputes silently, records
  `object$loglik_recomputed`, and (only for the full-sample fit, only
  when `verbose = TRUE`) emits a single message.
- A user-registered parallel backend (`doParallel`, `doFuture`, …) is
  left untouched when a function is called with `multicore = FALSE`.
- A non-finite value produced by an in-formula transformation
  (e.g. `log(x)` with `x <= 0`) now raises an informative error naming
  the column instead of misbehaving in the C++ code.
- User-supplied `init.Beta.NB` is used (regression-test hardened).
- `evinb` parameter count and `predict.*boot(pred = "original")`
  indexing fixes (carried over from the 0.9.4 audit follow-ups).
- [`predict()`](https://rdrr.io/r/stats/predict.html) on the `zinb` slot
  of a [`compare_models()`](../reference/compare_models.md) result with
  `type = "counts"` (or `type = "all"`, or
  `type = "counts", confint = TRUE`) no longer errors with “\$ operator
  is invalid for atomic vectors”.

### Breaking changes / deprecations

- The parallel backend is now **`future` / `furrr`** instead of
  `foreach` / `doParallel` / `doRNG` (which are no longer imported). Set
  a plan once for your session —
  `future::plan(future::multisession, workers = 8)` — and every
  bootstrap loop ([`evzinb()`](../reference/evzinb.md),
  [`evinb()`](../reference/evinb.md),
  [`add_bootstraps()`](../reference/add_bootstraps.md),
  `lr_test(bootstrap = TRUE)`,
  [`compare_models()`](../reference/compare_models.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`marginal_effects()`](../reference/marginal_effects.md)) uses it. The
  `multicore` argument now defaults to `NULL` (respect the current
  plan); `multicore = TRUE` still works as a shortcut that sets a
  temporary `multisession` plan for the one call and restores the
  previous plan on exit, and `multicore = FALSE` still forces sequential
  execution. See the new “Parallel processing” section in
  [`?evzinb`](../reference/evzinb.md).
- Because the RNG is now derived with `furrr`’s per-element L’Ecuyer
  streams, the bootstrap resamples drawn for a given `boot_seed` differ
  from earlier versions (the resampling distribution is unchanged, and
  results no longer depend on the number of workers). Pin results you
  need to reproduce.
- `future` and `furrr` are new hard dependencies; `foreach`,
  `doParallel`, `doRNG` and **`mistr`** are dropped.
- The extreme-value component of the mixture is now a single,
  package-wide **discretised Pareto** (`R/dist_pareto.R`): pmf for
  integer . Previously `predict(type = "quantile")` /
  `quantiles_from_*()` / `marginal_effects(type = "quantile")` built the
  mixture with `mistr` and inverted it there. The new path bisects the
  same mixture CDF used by `residuals(type = "quantile")` and the
  likelihood, so it also **corrects the negative-binomial dispersion**
  those quantiles used (they parametrised the NB by `Alpha.NB` rather
  than the canonical `size = 1 / Alpha.NB`). For an `evzinb` fit, where
  `Alpha.NB` is near 1, predicted quantiles move by at most a count or
  two; for an `evinb` fit with a strongly over-dispersed count component
  the old quantiles could be off by a large factor.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) /
  [`revzinb_fit()`](../reference/revzinb_fit.md) /
  [`revinb_fit()`](../reference/revinb_fit.md) draw the extreme-value
  part with `floor(C * U^(-1/alpha))` (was `round(...)`), so simulated
  extreme-value counts are ~0.5 lower on average and never below `C`.
- The minimum R version is now 4.1.0. The `NAMESPACE` uses delayed S3
  method registration for `insight` / `marginaleffects`, which needs R
  \>= 3.6.0; 4.1 is the common floor for packages using that pattern.
- The approximate bootstrap runtime estimate is printed only when
  `verbose = TRUE` (it was always printed before 0.9.4; this note was
  missing from the 0.9.4 changelog).
- `predict(type = "all")` and every `confint = TRUE` output now lead
  with the canonical state-probability columns (`pr_zero`, `pr_count`,
  `pr_evi`), keeping the deprecated `pr_zc` / `pr_pareto` as trailing
  duplicates — matching what `predict(type = "states")` already did.
- Default fits change slightly: with `c.lim = NULL` the candidate set
  for C_EV is data-driven and `init.C` defaults to the median of that
  set, so estimates from a default 0.9.x call are not bit-identical to a
  default 0.10.0 call. Pin `evinf_control(c.lim = ..., init.C = ...)` to
  reproduce an older fit.
- `gm_evzinb` gains a row (`obs_above_c_ev`) and its `parameter` row is
  now formatted with 0 decimals; regenerate any cached copy with
  [`gof_map_evinf()`](../reference/gof_map_evinf.md).
- The `hks` column `brv_AllLag` is renamed `brv_AllLag_log` (it is on
  the log scale, like the other `_log` columns). There is no alias —
  code referring to the old name will get an immediate “object not
  found” / “column not found” error.
- [`predict()`](https://rdrr.io/r/stats/predict.html) on the `nb` /
  `zinb` slot of a [`compare_models()`](../reference/compare_models.md)
  result no longer accepts `type = "evinf"` (a plain NB / ZINB has no
  extreme-value state);
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) now rejects it.
- The three internal `marginal.effect.*` helpers (never exported,
  superseded by
  [`marginal_effects()`](../reference/marginal_effects.md)) were
  removed.

### Documentation

- [`?vcov.evzinb`](../reference/vcov.evzinb.md),
  [`?confint.evzinb`](../reference/confint.evzinb.md) and the new
  `?marginaleffects-methods` page note that `c_ev` is estimated on a
  discrete grid, so continuous delta-method standard errors that perturb
  it should be read with care;
  [`marginal_effects()`](../reference/marginal_effects.md) (bootstrap)
  is the recommended route for effect uncertainty.
  [`?confint.evzinb`](../reference/confint.evzinb.md) gains an example
  using `c_ev` / `alpha_nb`.
- New `pkgdown` reference index grouping the exported functions
  (Fitting, Summaries & tables, Prediction & effects, Diagnostics &
  plots, Model comparison, Simulation, Data).
- [`?hks`](../reference/hks.md) now documents the columns as they
  actually ship (the previous help page listed several columns that are
  not in the data and omitted the `_log` transforms) and explains that
  `log1p(troopLag)` etc. can be used directly in the formula since
  0.9.4.
- GitHub Actions workflows for `R CMD check`, test coverage and the
  pkgdown site; coverage and check badges in the README.
- The per-observation fitted-value loop in the estimation routines is
  vectorised (no change to results).
- Internal: the twelve near-identical bootstrap-refit loops in
  [`compare_models()`](../reference/compare_models.md) are replaced by a
  single `boot_refit_family()` helper (output unchanged).
- Internal: the 1,800-line `R/a019.R` is split into six documented,
  de-duplicated files (`R/em_candidates.R`, `R/em_profile_c.R`,
  `R/em_step.R`, `R/em_fixed_c.R`, `R/em_fitted.R`, `R/em_fit.R`), one
  EM/ECME driver [`em_fit()`](../reference/em_fit.md) for both models
  instead of the two ~90%-duplicated `*.regression.fun()` pairs.
  Estimates are unchanged (a numerical-identity test suite pins the
  pre-split output). The three unused `prediction.*.fun()` legacy
  helpers were removed.
- `evinf_control(prune.c.range = ...)` now also thins the candidate set
  for [`evinb()`](../reference/evinb.md) (it was silently ignored there
  before), and the “more than 100 candidates” warning now fires for
  [`evinb()`](../reference/evinb.md) too.
- [`?hks`](../reference/hks.md): the `_log` columns are documented as a
  rule ([`log1p()`](https://rdrr.io/r/base/Log.html) for the personnel
  counts, [`log()`](https://rdrr.io/r/base/Log.html) for the episode
  duration), and the help notes that the bundled data carry no conflict
  identifier.
- [`?marginal_effects`](../reference/marginal_effects.md) /
  [`?tidy.evzinb`](../reference/tidy.evzinb.md) /
  [`?gof_map_evinf`](../reference/gof_map_evinf.md) gain guidance on
  harmonic-mean effect intervals and on the `modelsummary` `shape`
  argument.

## evinf 0.9.4

This release works through sections 1-3 of the internal package audit:
it fixes substantive statistical issues, a number of usability bugs, and
brings the documentation and metadata up to date. A `testthat` suite has
been added.

### Breaking changes / deprecations

- **Likelihood-ratio statistic.** [`lr_test()`](../reference/lr_test.md)
  now reports the standard `2 * (logLik_full - logLik_restricted)`
  statistic (and uses it for the bootstrap reference distribution).
  Statistics reported by earlier versions were half the correct value;
  p-values change accordingly but the qualitative conclusions in Randahl
  & Vegelius (2024) are unaffected. Degrees of freedom are now the
  number of design-matrix columns removed across components (so removing
  a four-level factor contributes three degrees of freedom).
- **Component vocabulary.** Model components are now named `"count"`,
  `"zero"`, `"evi"` and `"pareto"` everywhere
  ([`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  [`coefficient_extractor()`](../reference/coefficient_extractor.md),
  [`summary()`](https://rdrr.io/r/base/summary.html)). The old names
  `"nb"`, `"zi"` and `"evinf"` still work but emit a deprecation
  warning. `predict(type = )` additionally accepts `"zero"` for `"zi"`
  and `"evi"` for `"evinf"`.
- **[`tidy()`](https://generics.r-lib.org/reference/tidy.html) default**
  is now `component = "all"` (was `"zi"`), so a bare `tidy(model)` /
  `modelsummary(model)` shows every component. The `y.level` values from
  `tidy(component = "all")` are `zero`, `evi`, `count`, `pareto`.
- **Column renames** (old names kept as duplicates for this release
  only): `predict(type = "states")` columns are `pr_zero` / `pr_count` /
  `pr_evi`; `model$fitted$prob_pareto` / `posterior_pareto` become
  `prob_evi` / `posterior_evi`; the `props` / `resp` matrices are
  labelled `zero` / `count` / `evi`.
- **[`compare_models()`](../reference/compare_models.md)** returns the
  fitted model in `$model`; `$evzinb` remains as an alias. It now also
  accepts `evinb` objects (with `zinb_comparison` forced to `FALSE`).
- `generics` moved from `Depends` to `Imports`;
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) and
  [`glance()`](https://generics.r-lib.org/reference/glance.html) are
  re-exported, so [`library(broom)`](https://broom.tidymodels.org/) /
  [`library(generics)`](https://generics.r-lib.org) is no longer
  required.
- Dependencies `MLmetrics`, `stringr` and `stringi` removed.
- Minimum R version raised to 3.5 (data files re-serialised).

### Bug fixes

- [`compare_models()`](../reference/compare_models.md) fitted the
  bootstrapped ZINB with the zero-inflation formula for *both* model
  parts; it now uses the full `count | zero` formula, matching the
  full-sample fit.
- [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) no longer
  error on models fitted with `bootstrap = FALSE`; they return the point
  estimates with a message.
- `summary(p_value = "both")` now returns both the bootstrap and the
  approximate p-values (the approximate branch was previously
  unreachable).
- [`lr_test()`](../reference/lr_test.md) works on `evinb` objects (it
  passed the wrong data object before).
- Design matrices are built with
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) rather
  than `model.frame()[, -1]`, so factor and character covariates,
  interactions, [`poly()`](https://rdrr.io/r/stats/poly.html) / splines
  and in-formula transformations such as `log(x)` now work and are
  labelled correctly.
  [`na.omit()`](https://rdrr.io/r/stats/na.fail.html) is applied once to
  the union of all component variables, keeping the components
  row-consistent.
- User-supplied `init.Beta.NB` is no longer silently discarded when the
  NB and EVI formulas differ in length.
- [`evinb()`](../reference/evinb.md) bootstraps in the returned object
  are now named.
- The negative-probability guard in the predictor actually fires now.
- [`predict.nbboot()`](../reference/predict.nbboot.md) /
  [`predict.zinbboot()`](../reference/predict.zinbboot.md) are
  registered as S3 methods.
- [`summary()`](https://rdrr.io/r/base/summary.html) gains a `print`
  method; [`print.evzinb()`](../reference/print.evzinb.md) /
  [`print.evinb()`](../reference/print.evinb.md) show the correct
  component formulas plus `C_EV` and the count of observations at or
  above `C_EV`.
- [`glance()`](https://generics.r-lib.org/reference/glance.html) methods
  return unrounded values and gain `converged`, `n_bootstraps` and
  `n_failed_bootstraps` columns; `gm_evzinb` updated to match.
- Parallel backends are always unregistered on exit, and no “executing
  sequentially” / “Joining with `by`” messages are emitted.
- `run_evinb()`’s default `c.lim` matches
  [`evinb()`](../reference/evinb.md); `pdf.pl.type` is validated.
- The fitting routines check that the reported log-likelihood is the
  log-likelihood at the returned parameters and recompute `logLik` /
  `AIC` / `BIC` (with a warning) if not.
- `tidy(confint = )` is implemented (percentile or approximate
  intervals).
- [`revzinb_fit()`](../reference/revzinb_fit.md) /
  [`revinb_fit()`](../reference/revinb_fit.md) draw Pareto values in a
  single vectorised call.
- `evinb` objects no longer count the (fixed) zero-inflation
  coefficients as free parameters: `glance()$npar`, `AIC`, `BIC` and the
  approximate-test degrees of freedom are corrected accordingly.
- [`predict()`](https://rdrr.io/r/stats/predict.html) on `nbboot` /
  `zinbboot` objects with `pred = "original"` used a stray loop index
  and errored; it now predicts from the full-sample model.

### Documentation

- `inst/CITATION` and `inst/REFERENCES.bib` updated to *International
  Studies Quarterly* 68(4): sqae137 (2024), <doi:10.1093/isq/sqae137>.
- `hks`: `policeLag` and `militaryobserversLag` are now correctly
  described as UN police and military observers; all personnel counts
  are documented as being in thousands.
- `genevzinb` / `genevzinb2`: “independent” (not “dependent”)
  covariates; “EVZINB” spelling.
- `prune.c.range` documented consistently across all four fitting
  functions; `c.lim` reworded.
- `man/` regenerated with roxygen2 7.3.3; `src/a019_old.cpp` renamed to
  `src/evinf.cpp`; empty `src/coef_extractor.cpp` removed.
- All model-fitting examples wrapped in `\donttest{}` and reduced to
  `n_bootstraps = 5`, `multicore = FALSE`.
- `README` filled in.

## evinf 0.9.3

- Bug fixes.

## evinf 0.9.2

- [`predict()`](https://rdrr.io/r/stats/predict.html) can return the raw
  bootstrap draws (`return_bootstraps = TRUE`).

## evinf 0.9.1

- Added `prune.c.range` for thinning the candidate set of C when it
  contains many unique values.
- Informative error when every candidate value of C yields an infinite
  log-likelihood.
- Code formatting.

## evinf 0.9.0

- Package renamed from `evzinb` to `evinf`.
- [`revzinb_fit()`](../reference/revzinb_fit.md) /
  [`revinb_fit()`](../reference/revinb_fit.md) for drawing from a fitted
  model.
- Default parallel-processing behaviour of
  [`predict()`](https://rdrr.io/r/stats/predict.html) changed.
- Various minor fixes to non-standard evaluation and naming.

## evinf 0.8.0

- Initial CRAN submission.
