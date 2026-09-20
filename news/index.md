# Changelog

## evinf 0.10.0

Implements section 4 of the internal package audit (new functionality)
and the review follow-ups in `dev/review_round1.md`.

### New features

- Bootstrap replicates that would silently skew the bootstrap summaries
  are flagged (`$degenerate` / `$degenerate_reason`): a replicate whose
  EM did not converge; one with a non-finite fitted coefficient, or one
  larger in absolute value than `evinf_control(coef_limit = )` (default
  `50`, on the linear-predictor scale); or one whose smallest fitted
  Pareto shape is non-finite or below `evinf_control(alpha_floor = )`
  (default `0.001`).
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
- A bootstrap replicate whose `C_EV` estimate lands on an endpoint of
  the candidate grid is treated as degenerate (dropping those would bias
  the bootstrap distribution of `C_EV` inward); instead the count is
  stored as `object$n_c_on_boundary`, reported by
  [`glance()`](https://generics.r-lib.org/reference/glance.html) and the
  [`print()`](https://rdrr.io/r/base/print.html) boundary note, and a
  single [`warning()`](https://rdrr.io/r/base/warning.html) suggests
  widening `c.lim`.
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
  bootstrap differences (`evinf - compared`, negative favours the
  extreme-value model) in AIC, BIC and out-of-bag RMSE / RMSLE for an
  `evzinbcomp` object, with the proportion of bootstraps favouring the
  extreme-value model;
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

- `lr_test(bootstrap = TRUE)` no longer errors (“Tibble columns must
  have compatible sizes”) when any bootstrap replicate failed or was
  flagged degenerate. It now filters to the usable replicates first,
  like every other bootstrap summary, and gains
  `exclude_degenerate = TRUE` (audit0.10 §1.1); the results also report
  `n_bootstraps_used` next to `n_failed_bootstraps`.

- `summary(standard_error = FALSE)` / `tidy(standard_error = FALSE)` no
  longer error (“Column ‘se’/‘std.error’ not found”) when the model has
  bootstraps: `approx_t_value` is now silently forced to `FALSE` along
  with it, and `print.summary.*()` renders correctly with the
  `se`/`approx_t` columns missing (audit0.10 §1.2).

- `compare_models(razorize = TRUE)` refits NB/ZINB bootstraps on the
  wrong rows whenever the razorised data had fewer rows than the full
  data (which it always does): the bootstrap resample indices were used
  as-is against the smaller `data_razor`, mixing in unrelated or
  out-of-range rows. Resamples are now mapped onto `data_razor`’s own
  rows before refitting. Winsorised out-of-bag error is now computed
  against the raw (unwinsorised) outcome, not the winsorised one, so it
  is comparable to the EVZINB/EVINB OOB error.
  `inner_nb()`/`inner_zinb()`/`boot_refit_one()`/[`oob_evaluation()`](../reference/oob_evaluation.md)’s
  internal [`try()`](https://rdrr.io/r/base/try.html)s no longer print
  to the console on a failed replicate (audit0.10 §1.5).

- `marginal_effects(variables = )` errors on an unrecognised covariate
  name (listing the available ones) instead of silently returning a
  0-row tibble. A numeric covariate that only enters the model’s
  formulas wrapped in
  [`factor()`](https://rdrr.io/r/base/factor.html)/[`as.factor()`](https://rdrr.io/r/base/factor.html)/[`cut()`](https://rdrr.io/r/base/cut.html)
  (e.g. `y ~ factor(g_num)`) is now treated as categorical
  (level-vs-reference contrasts), instead of being perturbed as a number
  and producing unseen factor levels; one wrapped only in
  [`poly()`](https://rdrr.io/r/stats/poly.html)/`ns()`/`bs()` keeps the
  numeric derivative/difference path with its perturbation clamped to
  the covariate’s observed range. A covariate used both ways across
  formula components errors clearly (audit0.10 §1.10).

- A singular M-step Hessian no longer aborts the fit (or, on some BLAS
  backends, silently returns a wildly ill-conditioned step): the Newton
  step now uses `arma::solve(..., solve_opts::no_approx)`, which reports
  failure instead of throwing or returning garbage; on failure it
  retries once with a small ridge, and if that also fails the affected
  block keeps its pre-step value for that EM step (audit0.10 §1.3). This
  is also the usual reason a bootstrap replicate failed.

- The EM outer loop (the C_EV profile update) could in principle run
  forever if the profile oscillated between two candidate values.
  [`evinf_control()`](../reference/evinf_control.md) gains `max.c.iter`
  (default 50), capping each phase separately; hitting the cap in the
  convergence phase sets `converge = FALSE` (see the new `$c_converged`
  / [`glance()`](https://generics.r-lib.org/reference/glance.html)
  column to tell this apart from the inner EM not converging) and, for a
  full-sample fit, a [`warning()`](https://rdrr.io/r/base/warning.html)
  names the last two C_EV values visited.
  [`em_profile_c()`](../reference/em_profile_c.md) also now uses
  [`which.max()`](https://rdrr.io/r/base/which.min.html), so an exact
  tie in the profile no longer returns a length \> 1 `c_hat` and every
  candidate being `NaN`/infinite errors clearly instead of silently
  breaking the outer loop (audit0.10 §1.4).

- The likelihood could underflow to `-Inf` for a large count or a small
  Pareto shape: the discretised Pareto log-pmf now uses a
  cancellation-free form (`a*log(C/y) + log(-expm1(a*log(y/(y+1))))` in
  place of `log((C/y)^a - (C/(y+1))^a)`), the NB/Pareto mixture is
  combined with a log-sum-exp instead of summing on the natural scale,
  the NB log-pmf uses [`lgamma()`](https://rdrr.io/r/base/Special.html)
  instead of an O(y) loop, and the multinomial state probabilities (in
  the C++ M-step and in `prob_from_evzinb()` / `prob_from_evinb()`) use
  a numerically stable softmax. Coefficients move by at most ~2e-10 and
  log-likelihood by at most ~7e-12 on the identity fixtures – well under
  the 1e-6 the round’s numerical-identity gate would have required
  flagging (audit0.10 §1.11).

- [`evzinb()`](../reference/evzinb.md) /
  [`evinb()`](../reference/evinb.md) now validate the response after
  [`na.omit()`](https://rdrr.io/r/stats/na.fail.html): a negative,
  non-finite or non-integer value errors, naming the count and the
  problem, instead of being silently misread (a non-integer count is
  effectively ceiling()’d by the NB part of the C++ code but used
  exactly by the Pareto part; a negative one produced a cryptic
  “argument is of length zero” several layers down).
  `evinf_control(prune.c.range = )` now requires `[0, 1)` (was `[0, 1]`;
  `1` passed validation but crashed with “invalid ‘size’ argument”);
  [`em_c_candidates()`](../reference/em_c_candidates.md) clamps its
  internal sample size to `>= 0`, skips pruning entirely below 3
  candidates, and no longer risks
  [`sample()`](https://rdrr.io/r/base/sample.html)’s length-1-population
  trap (`sample(2, ...)` samples from `1:2`, not the single candidate
  `2`) on a 3-value grid. Pruning’s random draw is now reproducible
  without touching the caller’s `.Random.seed`, via a generalised
  `evinf_seeded_sample()` (audit0.10 §1.9).

- [`lr_test()`](../reference/lr_test.md) no longer emits
  [`evzinb()`](../reference/evzinb.md)/[`evinb()`](../reference/evinb.md)’s
  deprecation warning on every restricted refit (once per model, or once
  per model x bootstrap with `bootstrap = TRUE`). The restricted refits
  now pass a modified copy of the full model’s `control` object (its
  `init.*` fields set from the full-model estimate) instead of ~20
  individual tuning arguments, which also means a setting added to
  [`evinf_control()`](../reference/evinf_control.md) after the model was
  fitted (e.g. `max.c.iter`) is carried over automatically instead of
  needing to be added to a hand-picked argument list (audit0.10 §1.7).

- `evinf_control(pdf.pl.type = "exact")` was accepted and stored but had
  no effect: the M-step always used the continuous-Pareto
  gradient/Hessian for the Pareto block, giving bit-identical estimates
  to `"approx"`. `"exact"` now uses the analytic gradient/Hessian of the
  *discretised* Pareto log-pmf (the pmf the likelihood itself always
  uses) in the Newton step, computed in the same cancellation-free form
  as the likelihood (audit0.10 §1.11). The default stays `"approx"`, and
  default estimates are unchanged. See
  [`?evinf_control`](../reference/evinf_control.md) for what each option
  does (audit0.10 §1.8).

- A bootstrapped p-value of exactly `0` (no draw crossed the estimate)
  is now floored at `1 / B`, where `B` is the number of usable bootstrap
  replicates – the true p-value could be anywhere in `[0, 1/B)`, so
  reporting exactly `0` overstates the precision `B` draws can deliver.
  [`print.summary.evzinb()`](../reference/print.summary.evzinb.md) /
  [`print.summary.evinb()`](../reference/print.summary.evzinb.md) show
  this as `"< 1/B"` (e.g. `"<0.01"` for 100 usable bootstraps) instead
  of the default, misleadingly precise `"<2e-16"` (audit0.10 §1.11).

- `update(model, data = new_data)` kept the *old* data’s data-driven
  `control$c.lim` / `control$init.C` (both were written into
  `object$control` at the original fit), silently reusing a candidate
  range for chosen for data the model no longer uses. If the original
  range was itself data-driven (not pinned by the user via
  `evinf_control(c.lim = )`) and `data` is among the arguments being
  changed, they are now reset to `NULL` and re-resolved from the new
  data, with the usual message.
  [`?update.evzinb`](../reference/update.evzinb.md) also documents that
  changing `formula_nb.` does not propagate to a component formula that
  was originally left `NULL` and inherited from it – that component
  keeps its own formula regardless (audit0.10 §1.11).

- `object$props` / `object$resp` (and the `fitted$prob_*` /
  `posterior_*` vectors derived from them) came from the E-step at the
  *start* of the last EM step, one step behind the returned
  coefficients. They are now recomputed at the actual returned
  parameters and final C_EV before [`evzinb()`](../reference/evzinb.md)
  / [`evinb()`](../reference/evinb.md) return, using the same stable
  softmax as the C++ E-step. Point predictions
  (`predict(type = "harmonic")` etc.) are unaffected – only these
  diagnostic quantities change, by up to ~0.008 on the bundled example
  data (audit0.10 §1.11).

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a fitted
  model errored on ggplot2 \< 3.5.0 (`scale_*_continuous()` didn’t gain
  the `transform =` argument until 3.5.0). The four log1p-scale call
  sites in [`plot.evzinb()`](../reference/plot.evzinb.md) /
  [`plot.evinb()`](../reference/plot.evzinb.md) now pick `transform =`
  or the older `trans =` based on the installed ggplot2 version, and the
  `(>= 3.5.0)` floor on `ggplot2` in `Suggests` is removed (round8 0.1).

- `pdf.pl.type = "exact"`’s gradient/Hessian went to `NaN` for a sharply
  peaked Pareto block (large `alpha * |log(C/y)|`, e.g. `C = 100`,
  `y = 1e4`, `alpha ~ 403`) because `pareto_exact_derivs_fun()` formed
  `(C/y)^alpha` and `(C/(y+1))^alpha` separately in its numerator even
  though the denominator already used a cancellation-free form; both
  underflow to exactly 0 there, giving `0/0`. The `NaN` failed safe –
  the M-step just stopped updating the Pareto block silently – but is
  now fixed by factoring the same term out of the numerator (review §3,
  round8 0.2).

- `object$fitted$y.hat.pl_*` (the state-probability-weighted point
  predictions) were computed from the E-step’s one-step-stale prior
  state probabilities, while `object$props` / `object$resp` (and the
  `fitted$prob_*` / `posterior_*` vectors derived from them) already
  used the final, recomputed ones – so `fitted$y.hat.pl_E.inv.y` and
  `predict(type = "harmonic")` disagreed by a small but nonzero amount
  (`cor() = 0.9999999`, not 1). Both are now computed from the same
  final props (review §4, round8 0.3).

- [`oob_evaluation()`](../reference/oob_evaluation.md) on an
  `evzinbcomp` object masked a degenerate evinf bootstrap replicate’s
  out-of-bag error to `NA` in the `evinf` column only, leaving the
  compared `nb` / `zinb` columns unmasked at that row – so a column-wise
  `na.rm = TRUE` summary (as the vignette’s model-evaluation section
  does) compared medians computed over different replicate sets. Every
  column is now masked at the union of `NA` positions, and the row mask
  is exposed as `attr(out, "excluded")` (review §5, round8 0.4).

- [`compare_fit()`](../reference/compare_fit.md)’s `aic` / `bic` rows
  for a `*_razor` or `*_winsor` slot are not comparable to the evinf
  model’s – a razorised fit is estimated on fewer observations, and a
  winsorised fit on a different outcome – even though the `rmse` /
  `rmsle` rows now are (OOB error is always against the raw outcome).
  Those two metrics now come back `NA` for those slots, with a one-line
  footnote from [`print()`](https://rdrr.io/r/base/print.html) when any
  are present (review §6, round8 0.5).

- The warm-up phase’s C_EV profile could hit `max.c.iter` without
  settling, same as the convergence phase, but nothing recorded it –
  `$c_converged` and its
  [`warning()`](https://rdrr.io/r/base/warning.html) only ever covered
  the convergence-phase loop. Recorded now as `$c_warmup_capped`, with
  its own [`warning()`](https://rdrr.io/r/base/warning.html) for a
  full-sample fit; it does not affect `$converge` (review §7, round8
  0.6).

- `inner_nb()` / `inner_zinb()` (the per-bootstrap refits behind
  [`compare_models()`](../reference/compare_models.md)) computed
  out-of-bag error against `data[-boot_id, ]`, which silently returns
  zero rows rather than all of them when `boot_id` is empty – a resample
  where razorising leaves no drawn row behind. Both now guard explicitly
  and fail the same way a caught
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) /
  [`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html) error
  already does, so every existing `inherits(b, "try-error")` check picks
  it up (review §7, round8 0.7).

- CI: added an `ubuntu-latest` / `oldrel-1` job to `R-CMD-check.yaml` –
  the job that would have caught round8 0.1’s ggplot2-version failure
  (round7 prompt E.4, review §2, round8 0.8).

- The count-likelihood closed form
  `lgamma(y + 1/alpha) - lgamma(1/alpha) - lgamma(y + 1)` (round 7)
  partially cancels for a large `y` with a small `1/alpha` – checked
  against a Kahan-summed reference loop (plain summation of the
  pre-round-7 loop is not itself a trustworthy reference at this scale)
  at up to ~1.5e-10 relative error for `y` up to 1e6 and `alpha` down to
  `1e-3`. Switched to the algebraically identical, cancellation-free
  `-log(y + 1/alpha) - lbeta(y + 1, 1/alpha)`, which was at least as
  accurate in every case checked (review §7, round8 0.9). The
  Kahan-summed reference used to check this also showed that the
  *pre-round-7* loop (plain summation) itself drifts by ~1e-8 at large
  counts – so a fit of large counts from before 0.10.0 had a slightly
  wrong log-likelihood all along, not merely a slower one (round9 0.6).

- A fitted Pareto shape (`alpha_pl`) that collapses toward 0 made
  `predict(type = "explog")` return `Inf` and
  `predict(type = "harmonic")` return an astronomically large finite
  number (e.g. ~1e22 x `C` on `hks` with an ordinary specification)
  instead of erroring or warning – both involve `exp(1 / alpha_pl)` or
  `1 / alpha_pl`. A shared helper (`evinf_clamp_alpha_pl()`) now floors
  `alpha_pl` at `evinf_control(alpha_pl_floor = )` (default `0.01`) in
  `harmonic_calc()`, `explog_calc()`, the derived Pareto-tail summaries
  in `$fitted`
  ([`em_fitted_values()`](../reference/em_fitted_values.md)), and the
  quantile-prediction clamp that already existed for
  `quantiles_from_evzinb()` (extended to `quantiles_from_evinb()`, which
  previously lacked it). The helper warns once per call, naming how many
  observations were clamped – except inside
  [`em_fitted_values()`](../reference/em_fitted_values.md), which runs
  on every EM fit including every bootstrap replicate, so it clamps
  silently there rather than warning on routine fits.
  [`glance()`](https://generics.r-lib.org/reference/glance.html) gains
  `min_alpha_pl` (the true, unclamped value) and
  [`print()`](https://rdrr.io/r/base/print.html) adds a note when it
  falls below the floor, so a collapsed extreme-value shape stays
  visible without digging (review §2, round9 0.1).

- The round8 A.2 closed form for the NB dispersion Hessian’s alpha-alpha
  entry switched to its O(y) loop fallback only for `y <= 30`, but the
  cancellation it guards against is governed by the product
  `alpha_nb * y`, not `y` alone – the near-Poisson regime (small
  `alpha_nb`) is exactly where this mattered and exactly where the
  round8 grid (`alpha >= 1e-3`) never checked. A fresh grid down to
  `alpha = 1e-6` found the closed form 18% off a loop reference at
  `y = 31, alpha = 1e-6` (`alpha_nb * y = 3.1e-5`), 0.47% off at
  `y = 50`, and still 0.27% off at `y = 100` – errors a `y > 30` guard
  let straight through. The guard now switches on `alpha_nb * y > 0.1`,
  which keeps a comfortable margin above where the closed form’s
  relative error crosses `1e-9` (observed worst case ~1e-11 once
  `alpha_nb * y` clears that cutoff); a large `y` together with a
  vanishingly small `alpha_nb` still takes the O(y) loop rather than the
  (still available, but not implemented here) small-product series
  expansion – accepted per the round9 brief, since the loop is only
  slow, not wrong. This matters directly for the round9 0.11.0 Poisson
  count family (`alpha_nb -> 0` is its limit) (review §4, round9 0.4).

- `glance()$n_bootstraps` no longer counts every bootstrap replicate —
  it counts the **usable** ones (neither errored nor degenerate). The
  three columns `n_bootstraps`, `n_failed_bootstraps` and
  `n_degenerate_bootstraps` now partition the number of replicates
  requested. The [`gof_map_evinf()`](../reference/gof_map_evinf.md) row
  labels were retitled accordingly (“Usable bootstraps”, “Failed
  bootstraps”, “Degenerate bootstraps”) and given title case throughout;
  the `raw` names (i.e. the
  [`glance()`](https://generics.r-lib.org/reference/glance.html) column
  names) are unchanged.

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

- The likelihood, CDF, quantile prediction, residuals,
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) and the
  posterior state probabilities now all share a single **discretised
  Pareto** (`R/dist_pareto.R`): pmf for integer . The
  `'harmonic'`/`'explog'` point predictions still use the
  continuous-Pareto harmonic/geometric-mean formulas as a closed-form
  approximation (documented in
  [`?predict.evzinb`](../reference/predict.evzinb.md)), so they are the
  one place that is not computed from the discretised distribution.
  Previously `predict(type = "quantile")` / `quantiles_from_*()` /
  `marginal_effects(type = "quantile")` built the mixture with `mistr`
  and inverted it there. The new path bisects the same mixture CDF used
  by `residuals(type = "quantile")` and the likelihood, so it also
  **corrects the negative-binomial dispersion** those quantiles used
  (they parametrised the NB by `Alpha.NB` rather than the canonical
  `size = 1 / Alpha.NB`). For an `evzinb` fit, where `Alpha.NB` is near
  1, predicted quantiles move by at most a count or two; for an `evinb`
  fit with a strongly over-dispersed count component the old quantiles
  could be off by a large factor.
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

- **[`compare_fit()`](../reference/compare_fit.md) sign convention.**
  The paired bootstrap difference is now `evinf - compared` (was
  `compared - evinf`); `prop_evinf_better` is now `mean(diffs < 0)`. The
  print header and the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) x-axis label
  already said “negative favours evinf” — it was the computation that
  disagreed with them (a probe on `genevzinb2` showed a *positive*
  median with `prop_evinf_better = 1` for a model that clearly wins on
  AIC). [`compare_fit()`](../reference/compare_fit.md) /
  [`plot.evzinbcomp()`](../reference/compare_fit.md) also gain
  `exclude_degenerate = TRUE`, and
  [`oob_evaluation()`](../reference/oob_evaluation.md) now honours it
  for a single evinf model.

- The default `conf_level` for
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`marginal_effects()`](../reference/marginal_effects.md) is now `0.95`
  (was `0.9`), matching
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html). Pass
  `conf_level = 0.9` explicitly to keep the old default.

- **Deprecations planned for removal in 0.11.0:** the pre-0.9.4
  `component = "nb"/"zi"/"evinf"` aliases (use
  `"count"`/`"zero"`/`"evi"`); passing the individual EM tuning
  arguments (`max.diff.par`, `c.lim`, `init.C`, …) directly to
  [`evzinb()`](../reference/evzinb.md)/[`evinb()`](../reference/evinb.md)
  instead of through `control = evinf_control(...)`; and the deprecated
  `pr_zc`/`pr_pareto` column names and
  `fitted$prob_pareto`/`fitted$posterior_pareto` duplicate fields (use
  the canonical `pr_zero`/`pr_evi` and `prob_evi`/`posterior_evi`). All
  of these still work in 0.10.0 and most already warn when used.
  [`revzinb_fit()`](../reference/revzinb_fit.md)/[`revinb_fit()`](../reference/revinb_fit.md)
  are superseded by
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) but are not
  deprecated and have no planned removal.

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
- Internal (audit §5.4, round8 A.1):
  [`em_profile_c()`](../reference/em_profile_c.md) used to call
  `log_lik_fun()` once per `C_EV` candidate, redoing the
  state-probability and count-log-likelihood computation (which doesn’t
  depend on `C`) every time. New C++ `log_lik_profile_fun()` computes
  that shared work once and profiles the whole candidate grid in one
  call; [`em_profile_c()`](../reference/em_profile_c.md)’s return shape,
  [`which.max()`](https://rdrr.io/r/base/which.min.html) behaviour and
  all-infinite error are unchanged. Verified to `1e-10` against the old
  per-candidate calls, and the selected `c_hat` is unchanged (not just
  close) on every identity fixture.
- Internal (audit §5.4, round8 A.2): the two remaining O(y) loops in
  `delldtheta_nb_i_fun()` / `d2elldtheta2_nb_i_fun()` (the NB dispersion
  score and Hessian) are closed forms via
  [`digamma()`](https://rdrr.io/r/base/Special.html)/[`trigamma()`](https://rdrr.io/r/base/Special.html).
  The Hessian’s closed form cancels catastrophically for a small
  `alpha * y` (not just a small `y`), so it’s used only above a cutoff
  on that product, with the loop (already trivially cheap in that
  regime) kept below it – see round9 0.4 below, which corrected the
  cutoff variable this bullet originally (and wrongly) described as
  `y > 30`. Verified to `1e-9` against the old loops over `y` in
  `{0, 1, 5, 100, 1e4, 1e5}` x `alpha` in `{1e-3, 0.01, 0.5, 1, 5, 50}`;
  the M-step gradient/Hessian are unchanged on the identity fixtures.
- Internal (audit §5.4, round8 A.3): `log_lik_fun()` and
  `update_bfgs_fun()`’s per-observation loops rebuilt a linear predictor
  (`trans(X.submat(i, ...))`) and, for the four M-step BFGS blocks,
  allocated a fresh gradient/Hessian block via a per-row function call –
  on every row, every BFGS sub-iteration, and `log_lik_fun()` alone is
  called dozens of times per [`em_step()`](../reference/em_step.md)
  (once per evaluation of each of the four
  [`stats::optimise()`](https://rdrr.io/r/stats/optimize.html) line
  searches in `R/em_step.R`). Replaced with precomputed linear
  predictors (one matrix-vector product per block) and gradient/Hessian
  accumulation as `X' w` / `X' diag(w) X`. A full `hks` fit went from
  1.1x faster (A.1 + A.2 alone) to a platform-dependent 1.9x-8.4x faster
  (round9 0.5, review §3: the original single 8.4x/163x headline was
  measured on Apple Accelerate and did not reproduce on Linux/OpenBLAS,
  where the same `hks` fit was only 1.9x faster – quote a range, not one
  number, and see `inst/bench/bench_evinf.R`, now tracked and
  CI-runnable, for how to reproduce either end of it); see `A.4`’s
  benchmark table for the Apple Accelerate figures. Verified against the
  pre-A.3 code on 5 cases (`genevzinb2` and `hks`, both Pareto types,
  plus an `evinb` fit) at machine precision (max abs diff ~2e-12); the
  identity fixtures are unchanged at `1e-8`, with `c_hat` exactly (not
  just closely) unchanged on those cases. That comparison is on a single
  call at identical fixed inputs, not on a full fit’s trajectory: the
  reordering perturbs each Newton step at the ~1e-12 level, and on data
  with a close second optimum that is enough for the *iterated* EM to
  converge somewhere else. On one factor-covariate fit this moved the
  log-likelihood from -253.5553 to -253.6816 and `min(alpha_pl)` from
  1.91 to 0.973 on some platforms (review §1) – “estimates unchanged”
  above is therefore too strong as a blanket claim; read it as
  “unchanged on every case checked so far”, with a dedicated
  fixed-trajectory regression test now guarding that specific case
  (round9 0.3) and the `alpha_pl -> 0` consequence guarded separately
  (round9 0.1).

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
