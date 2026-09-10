# Average marginal effects for an evzinb / evinb model

Average, over the data, of the numerical effect of each covariate on the
chosen predicted quantity, with bootstrap confidence intervals obtained
by recomputing the effect on every bootstrap fit.

## Usage

``` r
marginal_effects(
  object,
  variables = NULL,
  type = c("harmonic", "states", "quantile"),
  quantile = NULL,
  at = list(),
  method = c("derivative", "difference"),
  eps = 1e-04,
  delta = 1,
  conf_level = 0.9,
  newdata = NULL,
  n_max = 500,
  exclude_degenerate = TRUE,
  multicore = NULL,
  ncores = NULL
)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` model with bootstraps.

- variables:

  Covariates to compute effects for (default: all raw covariates).

- type:

  `"harmonic"`, `"states"` or `"quantile"`.

- quantile:

  Quantile for `type = "quantile"`.

- at:

  Named list of values at which to hold covariates before differencing.

- method:

  For numeric covariates, `"derivative"` (central finite difference of
  step `eps`) or `"difference"` (average change from a `delta`-unit
  increase in the covariate). Defaults to `"derivative"`, except for
  `type = "quantile"` where it defaults to `"difference"` (a one-unit
  change is what a quantile effect means in practice). Factor covariates
  always use level-vs-reference differences.

- eps:

  Step size for the central finite difference (`method = "derivative"`).

- delta:

  Covariate increase for `method = "difference"` (default `1` unit; pass
  e.g. `sd(x)` for a one-SD change).

- conf_level:

  Confidence level for the bootstrap intervals.

- newdata:

  Data to average over (default: the estimation data).

- n_max:

  For `type = "quantile"` only: cap the number of rows the effect is
  averaged over (the per-observation quantile machinery is slow). The
  cap is applied automatically when
  `nrow(newdata) * n_bootstraps > 2e5`, or whenever `n_max` is passed
  explicitly; `n_max = Inf` disables it. The subsample is reproducible
  from the model's bootstrap seed.

- exclude_degenerate:

  Drop bootstrap replicates flagged degenerate (default `TRUE`); see the
  `alpha_floor` argument of [`evinf_control`](evinf_control.md).

- multicore:

  Bootstrap parallelisation shortcut. The default (`NULL`) respects
  whatever `future` plan is currently set (see the **Parallel
  processing** section). `TRUE` sets a temporary
  [`multisession`](https://future.futureverse.org/reference/multisession.html)
  plan for the duration of the call; `FALSE` forces sequential
  execution. The previous plan is always restored on exit.

- ncores:

  Number of workers when `multicore = TRUE`. Default (`NULL`) is one
  less than the number of available cores. Ignored when `multicore` is
  `NULL` or `FALSE`.

## Value

A tibble with columns `variable`, `contrast`, `type`, `estimate`,
`std.error`, `conf.low`, `conf.high`.

## Details

The confidence interval is a percentile interval: `conf.low` /
`conf.high` are the empirical quantiles of the effect recomputed on each
bootstrap fit. `std.error` is the standard deviation of those recomputed
effects, except for `type = "quantile"` where the bootstrap distribution
is heavy-tailed and `std.error` is instead a robust scale read off the
interval,
`(conf.high - conf.low) / (2 * qnorm(1 - (1 - conf_level)/2))`. It is
reported for reference and is not used to build the interval.

**Harmonic-mean effects and bootstrap intervals.** The harmonic-mean
prediction is \\C(1 + \alpha)/\alpha\\, which blows up on bootstrap fits
whose Pareto shape \\\alpha\\ is very small. The point estimate for
`type = "harmonic"` is well behaved, but the percentile interval can be
extremely wide because a handful of bootstrap draws are enormous. For
inference about how a covariate shifts the outcome, prefer
`type = "quantile"` (effect on a predicted quantile) or
`type = "states"` (effect on the state probabilities), whose bootstrap
distributions are bounded.

**Quantile effects.** The mixture quantile is a step function of a
count, so its derivative is either 0 or huge; `type = "quantile"`
therefore defaults to `method = "difference"` (the average change in the
predicted quantile from a `delta`-unit increase in the covariate). The
derivative is still available with `method = "derivative"` but its
bootstrap distribution can be very heavy-tailed. The quantile effect is
evaluated per observation by bisection on the mixture CDF, which is
still the slowest prediction type, so it is averaged over at most
`n_max` rows (see that argument).

## Parallel processing

Bootstrap fits (and the per-bootstrap work in
[`add_bootstraps`](add_bootstraps.md), [`lr_test`](lr_test.md),
[`compare_models`](compare_models.md),
[`predict.evzinb`](predict.evzinb.md) and `marginal_effects`) are
dispatched with furrr on top of a future plan. Set the plan once for
your session and leave `multicore` at its default:


    future::plan(future::multisession, workers = 8)
    progressr::handlers(global = TRUE)   # opt in to progress bars
    model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 1000)

Any `future` backend works (`multisession`, `cluster`, `callr`, a HPC
batchtools plan, ...). As a convenience `multicore = TRUE` sets a
temporary `multisession` plan for the single call and restores the
previous plan on exit. Large models may need
`options(future.globals.maxSize = <bytes>)` to raise the default limit
on the data shipped to each worker.

## Reproducibility

The bootstrap **resample indices** (`object$bootstraps[[i]]$boot_id`)
are fully determined by `boot_seed` and are independent of the future
backend, the number of workers, and chunking — a fixed seed always draws
the same resamples.

The **fitted coefficients** are reproducible only to within
floating-point noise. BLAS operations are not bit-reproducible across
processes, so a sequential run and a `multisession` run of the same seed
can differ in the last few digits; and because the EM log-likelihood can
have close local optima, an occasional replicate converges to a
different one under a different backend. For replication material,
record the `future` plan and
[`utils::sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
alongside `boot_seed`, and produce the canonical set of estimates under
a single fixed plan.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 10)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
#> Error in eval(expr, p) : inv(): matrix is singular
#> Error in eval(expr, p) : inv(): matrix is singular
#> Error in eval(expr, p) : inv(): matrix is singular
#> Warning: C_EV reached the boundary of the candidate range in 1 of 10 bootstrap replicates; consider widening c.lim.
marginal_effects(model, variables = "x1")
#> # A tibble: 1 × 7
#>   variable contrast type     estimate std.error conf.low conf.high
#>   <chr>    <chr>    <chr>       <dbl>     <dbl>    <dbl>     <dbl>
#> 1 x1       dydx     harmonic     140.     1645.     101.     2969.
# }
```
