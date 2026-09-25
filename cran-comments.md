## R CMD check results

0 errors | 1 warning | 3 notes

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
* The "unable to verify current time" note (round12) is this local
  machine's clock-verification service being unreachable from its network;
  a known, environment-specific `R CMD check --as-cran` note unrelated to
  the package, and not expected to reproduce on CRAN's own check machines.

Full local run (round12): examples (including `--run-donttest`) in ~105s,
tests in ~45s -- both essentially unchanged from round 11's ~90s/~43s
despite the new `hurdle_comparison` tests, since bootstrap-heavy files stay
`skip_on_cran()`-gated (everything still runs under CI, where
`NOT_CRAN=true`; the full CI suite, including the new `hurdle_comparison`
tests, is exercised there). Also verified on R 4.3 (the package's declared
floor) via a dedicated CI job added this round -- see the R-version note
below.

## R-version note (round12)

`lr_test()` crashed on R < 4.4 (`stats::drop.terms()` errors when a model
component's every term is dropped): fixed by building the emptied formula
by hand in that case instead of relying on `drop.terms()`'s version-
dependent tolerance for an empty result. Round 11 believed this already
fixed based on a reproduction run on R 4.5.2, which is why a CI job pinned
to the package's declared floor (`R (>= 4.1.0)`, tested here as R 4.3) was
added alongside the existing `release`/`oldrel-1` jobs -- `oldrel-1` floats
with each R release and does not itself guard the declared floor.

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
  competitor baselines, and (round12) a `pscl::hurdle()` competitor
  alongside the existing ZINB/ZIP one for a hurdle-family fit.
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
  loops, and a hardened `lr_test()` formula-reduction path -- offsets and
  no-intercept specifications preserved through formula reduction (round10),
  and the reduction no longer crashes on R < 4.4 when a component's every
  term is dropped (round12, with a CI job now pinned to the package's
  declared R floor to catch this class of regression) -- and (round11) a
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
