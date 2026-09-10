# Predictions from evinb object

Predictions from evinb object

## Usage

``` r
# S3 method for class 'evinb'
predict(
  object,
  newdata = NULL,
  type = c("harmonic", "explog", "counts", "pareto_alpha", "evinf", "count_state",
    "states", "all", "quantile"),
  pred = c("original", "bootstrap_median", "bootstrap_mean"),
  quantile = NULL,
  confint = FALSE,
  conf_level = 0.9,
  multicore = NULL,
  ncores = NULL,
  return_bootstraps = FALSE,
  exclude_degenerate = TRUE,
  ...
)
```

## Arguments

- object:

  An evinb object for which to produce predicted values

- newdata:

  Optional new data (tibble) to produce predicted values from

- type:

  Character string, 'harmonic' for the harmonic mean and 'explog' for
  exponentiated expected log, 'counts' for predicted count of the
  negative binomial component, 'pareto_alpha' for the predicted pareto
  alpha value, 'states' for the predicted component states (prior),
  'count_state' for predicted probability of the count state, 'evinf'
  for predicted probability of the pareto state, 'all' for all predicted
  values, and 'quantile' for quantile prediction.

- pred:

  Type of prediction to be used, defaults to the original prediction
  from the fitted model, with alternatives being the bootstrapped median
  or mean. Note that bootstrap mean may yield infinite values,
  especially when doing quantile prediction

- quantile:

  Quantile for which to produce quantile prediction

- confint:

  Should confidence intervals be made for the predictions? Note: only
  available for vector type predictions and not 'states' and 'all'.

- conf_level:

  What confidence level should be used for confidence intervals

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

- return_bootstraps:

  Should the bootstrapped predictions be returned as well? Useful for
  further custom analyses of the bootstrapped predictions.

- exclude_degenerate:

  Drop bootstrap replicates flagged degenerate (default TRUE); see the
  alpha_floor argument of evinf_control().

- ...:

  Other arguments passed to predict function

## Value

A vector of predicted values for type 'harmonic', 'explog', 'counts',
'pareto_alpha','evinf', 'count_state', and 'quantile' or a tibble of
predicted values for type 'states' and 'all' or if confint=T

## Parallel processing

Bootstrap fits (and the per-bootstrap work in
[`add_bootstraps`](add_bootstraps.md), [`lr_test`](lr_test.md),
[`compare_models`](compare_models.md),
[`predict.evzinb`](predict.evzinb.md) and
[`marginal_effects`](marginal_effects.md)) are dispatched with furrr on
top of a future plan. Set the plan once for your session and leave
`multicore` at its default:


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
model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
predict(model)
#>   [1]   53.549686   60.127971  410.351359   12.431620   39.194583   34.139420
#>   [7]  148.622606   99.329017   36.948813   40.651211    5.503769    1.450365
#>  [13]  261.660557   12.502587   30.214181   13.977638  481.714310   25.192966
#>  [19]  244.014161    4.705262   11.275945    7.750345   14.216249   13.849668
#>  [25]    1.705376   75.022677   12.948875   16.934129   10.507474   31.967200
#>  [31]   18.273388   23.133565   44.422700   39.316315    8.986171  110.014001
#>  [37]    5.114418   23.530260   12.002418   72.783160 1829.771456    4.552935
#>  [43]   12.543468   22.627168   36.658583  115.545186   57.298135   37.853021
#>  [49]  105.045485    3.806636   24.764057   40.402978   17.284007   56.161572
#>  [55]   89.050057   11.247267    2.132761  206.644661   93.762307   23.778490
#>  [61]   18.323217   52.785191    4.373497  115.842081   72.542800   28.570804
#>  [67]   12.004366   30.413402  235.022572   70.932515    2.937371   23.088212
#>  [73]   45.103709  149.571137    3.827968   19.295807   23.837272    6.856760
#>  [79]   14.309719   53.467792   33.234651   15.067921    8.666817   18.419649
#>  [85]   12.801270   30.326189    4.296430  123.654408    3.678261   23.667426
#>  [91]   13.833416   64.543036 1031.629755   77.681179   40.425301   21.961707
#>  [97]  139.887921   11.336808   16.955050   34.099904
predict(model, type='all', quantile = 0.9) # all available predicted values
#> # A tibble: 100 × 8
#>    harmonic explog   q90 pr_count pr_evi pr_pareto  count pareto_alpha
#>       <dbl>  <dbl> <dbl>    <dbl>  <dbl>     <dbl>  <dbl>        <dbl>
#>  1     53.5   58.8   229    0.888 0.112     0.112   22.2         1.56 
#>  2     60.1   60.1   184    0.955 0.0447    0.0447  54.3       682.   
#>  3    410.   410.    971    0.883 0.117     0.117  440.        160.   
#>  4     12.4   13.3     8    0.959 0.0408    0.0408   1.60        2.21 
#>  5     39.2   39.2   184    0.937 0.0628    0.0628  29.4       335.   
#>  6     34.1   34.2   186    0.879 0.121     0.121   12.5        23.5  
#>  7    149.   244.    668    0.865 0.135     0.135   88.5         0.524
#>  8     99.3  103.    314    0.833 0.167     0.167   66.3         2.31 
#>  9     36.9   37.5   192    0.896 0.104     0.104   15.1         4.38 
#> 10     40.7   40.8   185    0.939 0.0612    0.0612  29.3         5.98 
#> # ℹ 90 more rows
# }
```
