# Predictions from evzinb object

Predictions from evzinb object

## Usage

``` r
# S3 method for class 'evzinb'
predict(
  object,
  newdata = NULL,
  type = c("harmonic", "explog", "counts", "pareto_alpha", "zi", "evinf", "count_state",
    "states", "all", "quantile"),
  pred = c("original", "bootstrap_median", "bootstrap_mean"),
  quantile = NULL,
  confint = FALSE,
  conf_level = 0.95,
  multicore = NULL,
  ncores = NULL,
  return_bootstraps = FALSE,
  exclude_degenerate = TRUE,
  ...
)
```

## Arguments

- object:

  An evzinb object for which to produce predicted values

- newdata:

  Optional new data (tibble) to produce predicted values from

- type:

  Character string, 'harmonic' for the harmonic mean and 'explog' for
  exponentiated expected log, 'counts' for predicted count of the
  negative binomial component, 'pareto_alpha' for the predicted pareto
  alpha value, 'states' for the predicted component states (prior),
  'count_state' for predicted probability of the count state, 'evinf'
  for predicted probability of the pareto state,'zi' for the predicted
  probability of the zero state, 'all' for all predicted values, and
  'quantile' for quantile prediction.

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
'pareto_alpha','zi','evinf', 'count_state', and 'quantile' or a tibble
of predicted values for type 'states' and 'all' or if confint=T

## Details

The likelihood, CDF, quantile prediction, residuals and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) all use the
discretised Pareto distribution (the integer-valued distribution the
model is actually fit on). The `'harmonic'` and `'explog'` point
predictions instead use the harmonic and geometric means of the
\*continuous\* Pareto distribution as a closed-form approximation to the
corresponding discretised-Pareto moments; this keeps existing point
predictions unchanged but means they are not computed from exactly the
same distribution as the rest of the model.

`type = 'explog'` is \\C \cdot \exp(1/\alpha\_{PL})\\, which grows
explosively as the fitted Pareto shape \\\alpha\_{PL}\\ approaches 0 –
`evinf_control(alpha_pl_floor = )` (default 0.01) keeps it finite, but
not sane (\\C \cdot e^{100}\\ at the floor). A warning fires whenever
any \\\alpha\_{PL}\\ used in an explog prediction is below 0.1, since
the prediction is effectively undefined well before the floor is reached
(\\e^{10} \approx 2.2 \times 10^4\\); check `glance()$min_alpha_pl`, and
prefer `type = 'harmonic'` when this fires. No separate clamp is applied
beyond `alpha_pl_floor`.

## Parallel processing

Bootstrap fits (and the per-bootstrap work in
[`add_bootstraps`](add_bootstraps.md), [`lr_test`](lr_test.md),
[`compare_models`](compare_models.md), `predict.evzinb` and
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
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 2 of 5 bootstrap replicates; consider widening c.lim.
predict(model)
#> Warning: evinf (explog prediction): 3 fitted Pareto alpha values below 0.1; the geometric-mean prediction C * exp(1 / alpha_pl) is effectively undefined there (e.g. exp(10) ~= 2.2e4). Consider predict(type = "harmonic") instead.
#>   [1]   52.3536768   57.1668381  376.9452687   11.1639892   38.8170530
#>   [6]   35.9055806  126.3183684   90.7615235   37.3853032   34.0306144
#>  [11]    3.8757346    0.7290348  211.7443388   12.7340813   32.1514187
#>  [16]   14.0972436  391.8414182   25.8481011  220.0886230    3.4952654
#>  [21]    9.8384299    5.9820263   12.6066990   12.1342843    0.7762911
#>  [26]   76.3148526   14.1944891   18.5612434    9.6734937   28.8526825
#>  [31]   16.5098346   22.4290440   38.6944355   41.3541612    7.8551612
#>  [36]  114.3658284    3.1848768   22.0641625   10.4536677   71.7961908
#>  [41] 1339.4846821    3.5523495   10.0391497   19.7056928   28.4332799
#>  [46]  102.5571595   57.4677195   38.7979659  100.6043594    2.7228403
#>  [51]   20.6829272   41.1815094   15.9032426   55.6687633   95.2497251
#>  [56]   11.9472583    1.1241349  164.7064469   89.6777188   21.1429786
#>  [61]   18.3951387   52.2695259    3.2066429  119.4794258   66.6140900
#>  [66]   26.6852148   10.1881157   29.5384884  188.2644625   68.4199362
#>  [71]    1.8008339   20.0119216   45.7608924  122.9404738    2.6452961
#>  [76]   17.6827614   24.1374749    4.8513851   12.8528147   53.8795835
#>  [81]   33.1428461   14.5640965    7.6367303   16.1336634    7.4460037
#>  [86]   26.0772200    2.8713642  113.1886227    2.0980822   24.7794615
#>  [91]   11.0102490   64.8139949  752.3960458   74.1029728   38.7799200
#>  [96]   21.6369011  132.7022801    9.7755182   19.1935263   38.9803139
predict(model, type='all', quantile = 0.9) # all available predicted values
#> Warning: evinf (explog prediction): 3 fitted Pareto alpha values below 0.1; the geometric-mean prediction C * exp(1 / alpha_pl) is effectively undefined there (e.g. exp(10) ~= 2.2e4). Consider predict(type = "harmonic") instead.
#> # A tibble: 100 × 10
#>    harmonic explog   q90 pr_zero pr_count pr_evi pr_zc pr_pareto  count
#>       <dbl>  <dbl> <dbl>   <dbl>    <dbl>  <dbl> <dbl>     <dbl>  <dbl>
#>  1     52.4   57.3   216   0.450    0.422 0.128  0.450    0.128   35.9 
#>  2     57.2   57.2   184   0.712    0.255 0.0326 0.712    0.0326 200.  
#>  3    377.   377.   1227   0.500    0.388 0.112  0.500    0.112  917.  
#>  4     11.2   11.8     7   0.674    0.287 0.0385 0.674    0.0385   4.43
#>  5     38.8   38.8   182   0.600    0.344 0.0556 0.600    0.0556  83.0 
#>  6     35.9   35.9   186   0.377    0.488 0.135  0.377    0.135   20.1 
#>  7    126.   194.    486   0.445    0.401 0.154  0.445    0.154  131.  
#>  8     90.8   94.5   278   0.342    0.464 0.193  0.342    0.193   86.4 
#>  9     37.4   38.0   190   0.453    0.431 0.116  0.453    0.116   26.0 
#> 10     34.0   34.2   150   0.649    0.295 0.0559 0.649    0.0559  74.5 
#> # ℹ 90 more rows
#> # ℹ 1 more variable: pareto_alpha <dbl>
# }
```
