## R CMD check results

0 errors | 1 warning | 2 notes

* The warning ("checking whether package 'evinf' can be installed") is a
  compiler-flag warning from R's own bundled header
  (`R_ext/Boolean.h: unknown warning group '-Wfixed-enum-extension'`),
  produced by this local machine's clang version not recognising a pragma
  emitted by a newer R install; it does not originate in this package's
  code and is not expected to reproduce on CRAN's build machines.
* The "New maintainer" note is expected (see below); CRAN will contact the
  old address to confirm.
* The "checking HTML version of manual" note ("Skipping checking math
  rendering: package 'V8' unavailable") is a local-machine limitation (V8
  is not installed here) and is not expected to reproduce on CRAN's check
  machines.

Full local run: examples (including `--run-donttest`) in 90s, tests in
43s (down from ~265s before this round's `skip_on_cran()` pass on the
slowest test files -- everything still runs under CI, where
`NOT_CRAN=true`).

## Changes since the last CRAN release (0.8.10)

This is the first CRAN release since 0.8.10 (2024). 0.9.x, 0.10.0 and 0.11.0
were development versions never released to CRAN, so everything below is,
as far as CRAN is concerned, a change since 0.8.10 -- not a changelog for a
single prior release. The full per-version detail is in `NEWS.md`; this is
the arc:

* **Parallel backend:** `foreach`/`doParallel`/`doRNG` replaced by
  `future`/`furrr` (`multicore`/`ncores` arguments, `future::plan()`).
  Bootstrap resamples for a given seed differ from 0.8.10 (the resampling
  distribution itself is unchanged). Chunked scheduling and a runtime
  projection message were added for large bootstrap counts.
* **Model families:** `evinf_family()` generalises the count process
  (negative binomial or Poisson) and the zero process (zero-inflation
  mixture or hurdle); `compare_models()` gained matching Poisson/ZIP
  competitor baselines.
* **Weights and offsets:** `weights =` on `evzinb()`/`evinb()`, carried
  through bootstrap resampling, `lr_test()`'s restricted refits, and every
  `compare_models()` competitor fit (including winsorised/razorised
  variants and their own bootstrap refits). `offset()` is supported in each
  model component's formula and is preserved through formula reduction
  (`lr_test()`) and competitor-formula construction (`compare_models()`).
* **Panel/time-series resampling:** `bootstrap_scheme =` adds
  `"cluster"`/`"moving_block"`/`"stationary"` resampling (`block =`/
  `time =`) alongside the original i.i.d. bootstrap, with validation that
  fires before fitting rather than surfacing as every replicate failing.
* **Diagnostics:** `evinf_control(n_starts = )` for multi-start EM fits
  (guards against the likelihood's local optima), `plot(type = "trace")`,
  and `check_evinf()`, a single diagnostic report over convergence,
  `C_EV` boundary behaviour, the fitted Pareto shape floor, and bootstrap
  health.
* **Forecasting types:** `predict()` gained `type = "distribution"`
  (the full predictive pmf), `"quantile"` (vectorised over several
  probabilities at once), `"exceedance"` (`P(Y >= k)`), and `"draws"`
  (predictive samples, optionally propagating bootstrap parameter
  uncertainty), all built on one family-aware predictive-distribution
  backbone shared with `residuals()`, `simulate()` and `classify_states()`.
  `predict(type = "explog")` gained `clamp_alpha_pl` (round11 A4) after a
  release-blocking review found the type could silently return an
  astronomical or infinite value on a collapsed fitted Pareto shape.
* **Ecosystem integrations:** `marginaleffects` (`get_predict()` for the
  new prediction types), `augment()` (the `broom`/`generics` convention),
  `emmeans` (count-component basis, via delayed S3 registration), and
  `modelsummary` tables via `tidy()`/`glance()`.
* **Numerous EM correctness fixes:** a numerically stable log-space
  likelihood, guarded Hessian solves in the C++ optimiser, capped profiling
  loops, and (round11) a hardened `lr_test()` formula-reduction path and a
  `predict()` argument-validation pass that rejects unknown arguments and
  warns on ones irrelevant to the requested `type`. Default (unweighted)
  fits are bit-identical to 0.8.10 on every model this was checked against.
* **Breaking changes:** `compare_fit()`'s paired-bootstrap sign convention
  flipped to `evinf - compared`; default `conf_level` for `predict()`/
  `marginal_effects()` changed from 0.9 to 0.95 (matching `confint()`/
  `tidy()`); the legacy `object$fitted$y.hat.pl_*` fields were dropped in
  favour of `predict()`/`fitted()`. See `NEWS.md` for the complete list.
* **Dropped dependencies:** `mistr`, `foreach`, `doParallel`, `doRNG`,
  `MLmetrics`, `stringr`, `stringi`.
* **New dependencies:** `future`, `furrr`, `parallelly`, `progressr`.

## Maintainer e-mail change

The maintainer's e-mail address has changed from `david.randahl@pcr.uu.se`
to `david.randahl@fhs.se`. CRAN will contact the old address to confirm
this submission.

## Reverse dependencies

`evinf` has no reverse dependencies on CRAN, verified with
`tools::package_dependencies("evinf", reverse = TRUE, db = available.packages())`
(returned `character(0)`).
