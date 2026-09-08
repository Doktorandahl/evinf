# Prompt for Claude Code: evinf cleanup round (audit sections 1–3)

You are working in the R package `evinf` (extreme value and zero-inflated negative binomial regression). The repository root is the working directory. Read `dev/package_audit.md` first; this task covers **sections 1, 2, and 3** of that document only. Section 4 (new functionality) is out of scope — do not add features beyond what is needed to fix the listed issues.

## Ground rules

- Work in plan mode first. Produce a written plan that lists every audit item in §1–§3 with (a) the file(s) and function(s) it touches, (b) the intended fix in one or two sentences, (c) whether the fix changes user-facing behaviour or output, and (d) how it will be tested. Wait for my approval before editing anything.
- I code in R with tidyverse idioms. Keep the existing style (`dplyr::`, `purrr::`, explicit namespacing, roxygen2 docs). Do not reformat files you are not otherwise changing.
- Every bug fix gets a `testthat` test that fails on the current code and passes after the fix. Create `tests/testthat/` (it does not exist). Use `genevzinb2` (100 rows) for fast tests and keep any single test under ~10 seconds; use `n_bootstraps` of 5–10 and `multicore = FALSE` in tests. Use `hks` only in tests explicitly marked `skip_on_cran()`.
- After the code changes: run `devtools::document()` to regenerate `NAMESPACE` and `man/`, then `devtools::check()`. Report the check output (errors, warnings, notes) verbatim in your final message. Do not suppress notes by editing `.Rbuildignore` unless the note is about a file that genuinely should be ignored.
- Do not change the C++ code in `src/` except to rename `a019_old.cpp` (see §3.10). Do not change the estimation algorithm in `R/a019.R` except where an audit item explicitly points there (§2.13).
- Do not touch `vignettes/` in this round.
- Bump the version to 0.9.4 and write a `NEWS.md` entry that lists every user-visible change, grouped as "Bug fixes", "Breaking changes / deprecations", and "Documentation". Bring `NEWS.md` up to date for 0.9.x in general using `git log`.
- Commit in small, logically separate commits (one audit item or one closely related group per commit) with messages that reference the audit item number, e.g. `Fix LR statistic (audit 1.1)`. Do not push.

## Decisions already made (apply these; do not ask again)

**§1.1 LR statistic.** Use `2 * (loglik_full - loglik_restricted)` everywhere in `lr_test()`, including the bootstrap branch (`ks.test`, `chisq_mean`, `chisq_median`, `prop_sig`, and the per-bootstrap `statistic` column in `boot_results`). Add a test asserting `statistic == 2 * (logLik_full - logLik_restricted)` and `prob == pchisq(statistic, df, lower.tail = FALSE)`. Mention in NEWS that previously reported LR statistics were half the correct value.

**§1.2 compare_models ZINB formula.** Build `f_zinb` once and pass it to `inner_zinb()`; `inner_zinb()` must use it for the bootstrap fits. Test: with `formula_zi` a strict subset of `formula_nb`, the bootstrap ZINB objects must have `length(coef(boot)$zero) == length(all.vars(formula_zi))` (intercept included).

**§1.3 summary p_value = "both".** Two independent `if` blocks.

**§1.4 summary/tidy on non-bootstrapped models.** These must work. When `object$bootstraps` is `NULL`: `standard_error`, `p_value`, `approx_t_value`, `bootstrapped_props`, and `confint` are silently treated as off (return estimate-only tables), `n_failed_bootstraps` is `NA_integer_`, and a single `message()` says that bootstrap-based quantities are unavailable. `coef_type`/`coef` values other than `"original"` should `stop()` with an informative error when there are no bootstraps.

**§1.5 / §1.6 evinb support.** `lr_test()` must use `object$data$data` for evinb. `compare_models()` must accept evinb objects: for evinb set `zinb_comparison` to `FALSE` with a message if the user left it at the default (`TRUE`), and `stop()` if the user explicitly asked for it. Rename the first slot of the returned list from `evzinb` to `model`, but keep `evzinb` as an alias (i.e. both `out$model` and `out$evzinb` work — implement via a second list element pointing to the same object, or document the rename and update `print.evzinbcomp`). Update every internal reference and the examples.

**§1.7 model.matrix.** This is the largest change. Implement it fully in this round:
- In `run_evzinb()`, `run_evinb()`, `prob_from_evzinb()`, `prob_from_evinb()`, `counts_from_evzinb()`, `fitted_alpha_from_evzinb()`, `quantiles_from_*()`, and the `bootrun_*` functions, construct design matrices with `stats::model.matrix(formula, data)` and drop the `"(Intercept)"` column (the C++ code prepends the intercept). Coefficient names come from `colnames(model.matrix(...))`.
- Store the `terms` object of each component in the fitted object (`object$terms$nb`, `$zi`, `$evi`, `$pareto`) and use `model.frame(terms, newdata, xlev = ...)` + `model.matrix()` for `newdata` so factor levels and contrasts are consistent between fit and prediction. Store `xlevels` per component.
- `formula_var_remover()` in `lr_tests.R` must still work when a term is `log(x)`, `x1:x2`, or a factor; removing variable `x` should remove every term whose `all.vars()` contains `x`. Degrees of freedom must be counted as the number of *columns* dropped from the design matrix across components (a 4-level factor removed from one component is 3 df), not the number of terms.
- `compare_models()` must pass the same formulas to `MASS::glm.nb` and `pscl::zeroinfl` (they already accept full formula syntax).
- `evzinb()`/`evinb()`: `na.omit` behaviour must remain row-consistent across all four components and the `block` variable (currently handled by selecting `all.vars()` of all formulas; keep that logic but base it on the union of `all.vars()` from the four formulas).
- Tests: (i) fit with `y ~ x1 + x2 + x3` gives numerically identical coefficients before and after the refactor (write the "before" values into the test from the current package version); (ii) `y ~ x1 + log(abs(x2) + 1) + x1:x3` fits and names coefficients `log(abs(x2) + 1)` and `x1:x3`; (iii) a factor covariate with three levels yields two dummy columns and `predict(newdata = )` works when newdata contains a subset of the levels; (iv) `predict(model)` equals `predict(model, newdata = model$data$data)`.
- Offsets are *not* required in this round, but do not write code that would make adding `offset()` later harder (i.e. keep `model.offset()` retrievable from the stored model frame or terms).

**§1.8** Use `ncol(mf_nb)`-equivalent (now: number of NB design columns + 1).

**§2.1, 2.2, 2.3, 2.4, 2.5, 2.10, 2.12, 2.14, 2.15** Fix as described in the audit. For §2.3 the intended guard is `if (min(pr_count) < -1e-8) stop(...) else pr_count[pr_count < 0] <- 0`. For §2.5 add `@export` and an Rd for `predict.nbboot`. For §2.10 make both `glance` methods return unrounded values and add columns `converged` (logical), `n_bootstraps`, and `n_failed_bootstraps` (both `NA` when not bootstrapped); update `gm_evzinb` (in `data/`; regenerate the `.rda` with a script saved under `data-raw/`) with rows for these. For §2.15 vectorise with `mistr::rpareto(n, C_est, alphs)`.

**§2.6 print.summary.** Add `print.summary.evzinb()` and `print.summary.evinb()`. Layout: header line with model type and formulas; one block per component with columns Estimate, Std. Error, approx t, p-value (whichever were computed), using `printCoefmat()` with signif. stars; then `alpha_NB`, `C_EV`, mean state proportions, observations, parameters, log-likelihood, AIC, BIC, number of bootstraps and failed bootstraps. Keep `print.evzinb()` short but add `C_EV` and the number of observations at or above `C_EV`.

**§2.7 component vocabulary.** Canonical values are `"count"`, `"zero"`, `"evi"`, `"pareto"`, `"all"` for the `component` argument in `tidy()`, `coefficient_extractor()`, and (new) `summary()`'s list names (`count`, `zero`, `evi`, `pareto`). Accept the old values (`"nb"`, `"zi"`, `"evinf"`) as aliases via an internal `normalize_component()` helper that emits a `lifecycle`-style `warning()` (plain `warning()`, do not add lifecycle as a dependency) saying the alias is deprecated. `predict(type = )` keeps its current values (`"counts"`, `"zi"`, `"evinf"`, `"count_state"`, `"pareto_alpha"`, `"states"`) but also accepts `"zero"` as an alias for `"zi"` and `"evi"` for `"evinf"`. Column names of `predict(type = "states")` become `pr_zero`, `pr_count`, `pr_evi`; `props`/`resp` colnames become `zero`, `count`, `evi`; `fitted$prob_pareto`/`posterior_pareto` become `prob_evi`/`posterior_evi` — with the old names kept as duplicates for one release and listed under deprecations in NEWS. The `y.level` values from `tidy(component = "all")` become `zero`, `evi`, `count`, `pareto` in that order.

**§2.8** `tidy()` default `component = "all"`.

**§2.9 tidy confint.** Implement it: `confint = "bootstrapped"` gives percentile intervals from the bootstrap distribution at `conf_level`; `"approx"` gives `estimate ± qt(1 - (1 - conf_level)/2, df = nobs - npar) * std.error`. Columns named `conf.low`, `conf.high`. Same for `tidy.zinbboot`/`tidy.nbboot` if they have the argument.

**§2.11 parallel backends.** Every function that registers a backend must unregister it on exit (`on.exit(doParallel::stopImplicitCluster())` / `parallel::stopCluster(cl)` + `foreach::registerDoSEQ()`). When `multicore = FALSE`, use `%do%` (not `%dopar%`) so no "executing sequentially" warnings appear. Factor the register/unregister logic into one internal helper used by `evzinb()`, `evinb()`, `compare_models()`, `lr_test()`, and `predict.*()`.

**§2.13** Add an internal check that the log-likelihood at the returned parameters equals `out$log.lik` within `1e-6`; if not, recompute `log.lik`, `AIC`, `BIC` from the returned parameters and issue a `warning()`.

**§3 documentation.** Fix every item in §3.1–§3.8. For §3.1 the reference is: Randahl, David, and Johan Vegelius. 2024. "Inference with Extremes: Accounting for Extreme Values in Count Regression Models." *International Studies Quarterly* 68(4): sqae137. doi:10.1093/isq/sqae137. Bib key stays `randahl2023evzinb` in `REFERENCES.bib` only if renaming would break `\insertCite` elsewhere; otherwise rename to `randahl2024inference` and update all uses. For §3.9: move `generics` from `Depends` to `Imports` and re-export `tidy` and `glance` (`@export` on `generics::tidy`/`generics::glance` re-exports); remove `MLmetrics`, `stringr`, and `stringi` (replace with base R: `rmsle <- function(pred, true) sqrt(mean((log1p(pred) - log1p(true))^2))` etc.; `regexpr`/`grepl`/`strsplit`); add `Suggests: testthat (>= 3.0.0), modelsummary, ggplot2`; fix the Description text; **do not** change the maintainer email. For §3.10: add `^\.RData$`, `^\.Rhistory$`, `^dev$`, `^data-raw$` to `.Rbuildignore`; rename `src/a019_old.cpp` → `src/evinf.cpp` and delete the empty `src/coef_extractor.cpp`; rebuild `RcppExports` with `Rcpp::compileAttributes()`; do not delete `.RData` from disk (just ignore it).

**Examples.** Wrap every example that fits a model in `\donttest{}` and reduce to `n_bootstraps = 5`. Replace the examples that use `multicore = TRUE, ncores = 2` with `multicore = FALSE`. Fix the wrong examples noted in §3.7.

## What to report back

1. The approved plan with each item marked done / partially done / skipped, with a reason for anything not fully done.
2. `devtools::check()` output.
3. A table of every user-visible behaviour change (this feeds NEWS and the vignette).
4. Any place where you were unsure whether existing behaviour was intentional — list it rather than guessing.
