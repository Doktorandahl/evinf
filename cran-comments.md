## R CMD check results

0 errors | 1 warning | 0 notes

* The warning ("checking whether package 'evinf' can be installed") is a
  compiler-flag warning from R's own bundled header
  (`R_ext/Boolean.h: unknown warning group '-Wfixed-enum-extension'`),
  produced by this local machine's clang version not recognising a pragma
  emitted by a newer R install; it does not originate in this package's
  code and is not expected to reproduce on CRAN's build machines.

## Changes since the last CRAN release (0.8.10)

This is a substantial update. Highlights:

* **Breaking changes:** the parallel backend is now `future`/`furrr` (the
  `multicore`/`ncores` arguments and `future::plan()` replace `foreach`/
  `doParallel`/`doRNG`); bootstrap resamples for a given seed therefore
  differ from earlier versions (the resampling distribution itself is
  unchanged). `compare_fit()`'s paired-bootstrap sign convention flipped to
  `evinf - compared`. Default `conf_level` for `predict()`/
  `marginal_effects()` changed from 0.9 to 0.95, matching `confint()`/
  `tidy()`. See `NEWS.md` for the complete list (EM tuning settings moved
  into `evinf_control()`, the `C_EV` candidate range is data-driven by
  default, quantile prediction/simulation now use one discretised-Pareto
  implementation package-wide instead of `mistr`, and more).
* **Dropped dependencies:** `mistr`, `foreach`, `doParallel`, `doRNG`,
  `MLmetrics`, `stringr`, `stringi`.
* **New dependencies:** `future`, `furrr`, `parallelly`, `progressr`.
* Numerous correctness fixes to the EM fitting routine (numerically stable
  log-space likelihood, guarded Hessian solves in the C++ optimiser,
  capped profiling loops) and several smaller bug fixes throughout
  (documented in `NEWS.md`).

## Maintainer e-mail change

The maintainer's e-mail address has changed from `david.randahl@pcr.uu.se`
to `david.randahl@fhs.se`. CRAN will contact the old address to confirm
this submission.

## Reverse dependencies

`evinf` has no reverse dependencies on CRAN, verified with
`tools::package_dependencies("evinf", reverse = TRUE, db = available.packages())`
(returned `character(0)`).
