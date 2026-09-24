# Predictions from evinb object

Predictions from evinb object

## Usage

``` r
# S3 method for class 'evinb'
predict(
  object,
  newdata = NULL,
  type = c("harmonic", "explog", "counts", "pareto_alpha", "evinf", "count_state",
    "states", "all", "quantile", "distribution", "exceedance", "draws"),
  pred = c("original", "bootstrap_median", "bootstrap_mean"),
  quantile = NULL,
  confint = FALSE,
  conf_level = 0.95,
  multicore = NULL,
  ncores = NULL,
  return_bootstraps = FALSE,
  exclude_degenerate = TRUE,
  support = NULL,
  max_support = 1e+05,
  format = c("long", "matrix"),
  threshold = NULL,
  n_draws = 1000,
  parameter_uncertainty = FALSE,
  seed = NULL,
  keep = NULL,
  clamp_alpha_pl = FALSE,
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
  values, 'quantile' for quantile prediction (scalar or vector, see
  \`quantile\`), 'distribution' for the full predictive distribution,
  'exceedance' for exceedance probabilities (see \`threshold\`), and
  'draws' for predictive draws (see \`n_draws\`).

- pred:

  Type of prediction to be used, defaults to the original prediction
  from the fitted model, with alternatives being the bootstrapped median
  or mean. Note that bootstrap mean may yield infinite values,
  especially when doing quantile prediction

- quantile:

  Quantile(s) for which to produce quantile prediction. A single value
  keeps the existing scalar behavior exactly (a vector, e.g. \`c(.5, .9,
  .99)\` (round10 H.3), returns a tibble with one \`qXX\` column per
  probability (\`qXX_lo\`/\`qXX_hi\` too, with \`confint = TRUE\`)).

- confint:

  Should confidence intervals be made for the predictions? Note: only
  available for vector type predictions and not 'states', 'all' or
  'distribution'.

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

- support:

  (round10 H.2) \`type = 'distribution'\` only: shared support (vector
  of \`y\` values) for every row. \`NULL\` (default) uses \`0:K\`, \`K\`
  the ceiling of the 0.999 mixture quantile of the heaviest-tailed row,
  capped at \`max_support\`.

- max_support:

  (round10 H.2) Upper bound on the default \`support\`'s \`K\` (default
  1e5); ignored when \`support\` is given explicitly.

- format:

  (round10 H.2) \`type = 'distribution'\` only: \`'long'\` (default; a
  tibble with \`.row\`, \`y\`, \`prob\`) or \`'matrix'\` (the raw n x K
  matrix; \`keep\` is ignored).

- threshold:

  (round10 H.4) \`type = 'exceedance'\` only: a vector of thresholds;
  result columns are \`p_ge\_\<threshold\>\`, \`P(Y \>= threshold)\`.

- n_draws:

  (round10 H.5) \`type = 'draws'\` only: number of predictive draws per
  row.

- parameter_uncertainty:

  (round10 H.5) \`type = 'draws'\` only: if \`TRUE\`, each draw uses a
  randomly chosen usable bootstrap replicate's parameters instead of the
  full-sample estimate.

- seed:

  (round10 H.5) \`type = 'draws'\` only: optional seed; the caller's
  \`.Random.seed\` is left untouched either way (same convention as
  \`simulate()\`).

- keep:

  (round10 H.6) Optional character vector of \`newdata\` column names to
  carry into the result (a join key for panel data); supported for
  \`type\` in \`'quantile'\` (vector), \`'distribution'\` (\`format =
  'long'\`), \`'exceedance'\` and \`'draws'\`.

- clamp_alpha_pl:

  (round11 A4) \`type = 'explog'\` (and \`'all'\`) only: \`FALSE\`
  (default) uses the fitted \`alpha_pl\` as-is and warns below 0.1;
  \`TRUE\` clamps at 0.1; a positive number clamps there instead. Exempt
  from the global \`alpha_pl_floor\` (which stays in force for
  \`'harmonic'\`/\`'quantile'\`). The clamp used is recorded as a
  \`"clamp_alpha_pl"\` attribute on the result.

- ...:

  Not used; any name here is an unknown argument and errors (round11
  A5), naming the valid arguments for the requested \`type\`. A known
  argument supplied for a \`type\` it does nothing for (e.g.
  \`threshold\` with \`type = "harmonic"\`) warns instead, since it has
  a named formal below and never reaches here.

## Value

A vector of predicted values for type 'harmonic', 'explog', 'counts',
'pareto_alpha','evinf', 'count_state', and 'quantile' (scalar), or a
tibble for type 'states', 'all', 'quantile' (vector), 'distribution'
(\`format = 'long'\`), 'exceedance', 'draws', or if confint=T; a matrix
for 'distribution' with \`format = 'matrix'\`.

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
explosively as the fitted Pareto shape \\\alpha\_{PL}\\ approaches 0.
Unlike `'harmonic'` and `'quantile'`, this path is exempt from the
global `evinf_control(alpha_pl_floor = )`; `clamp_alpha_pl` is its own,
separate knob (round11 A4). With `clamp_alpha_pl = FALSE` (the default),
the fitted \\\alpha\_{PL}\\ is used as-is and a warning fires whenever
any value is below 0.1, since the prediction is effectively undefined
there (\\e^{10} \approx 2.2 \times 10^4\\, and it may be `Inf`); check
`glance()$min_alpha_pl`, and either pass `clamp_alpha_pl = TRUE` (clamps
at 0.1) or a positive number (clamps there instead), or use
`type = 'harmonic'`. The clamp actually applied (`FALSE`, or the numeric
floor) is recorded as a `"clamp_alpha_pl"` attribute on the result.

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
#> Warning: evinf (explog prediction): 3 fitted Pareto alpha values below 0.1 (smallest: 0.0108); the geometric-mean prediction C * exp(1 / alpha_pl) is effectively undefined there (e.g. exp(10) ~= 2.2e4) and may be Inf. Consider predict(..., type = "explog", clamp_alpha_pl = TRUE).
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
predict(model, type='quantile', quantile = c(.5, .9, .99)) # several quantiles (H.3)
#> # A tibble: 100 × 3
#>      q50   q90   q99
#>    <dbl> <dbl> <dbl>
#>  1     0   229   886
#>  2     0   184   848
#>  3     9   971  6625
#>  4     0     8   347
#>  5     0   184   456
#>  6     0   186   216
#>  7     2   668 26323
#>  8     2   314  1114
#>  9     0   192   346
#> 10     0   185   460
#> # ℹ 90 more rows
predict(model, type='distribution') # the full predictive distribution (H.2)
#> evinf (predict distribution): the default support's upper bound (86294293458, the 0.999 mixture quantile of the heaviest-tailed row) exceeds max_support (1e+05); truncating there. Pass `support = ` explicitly for a wider (or narrower) support.
#> # A tibble: 10,000,100 × 3
#>     .row     y    prob
#>    <int> <int>   <dbl>
#>  1     1     0 0.518  
#>  2     1     1 0.0513 
#>  3     1     2 0.0281 
#>  4     1     3 0.0196 
#>  5     1     4 0.0151 
#>  6     1     5 0.0123 
#>  7     1     6 0.0104 
#>  8     1     7 0.00905
#>  9     1     8 0.00799
#> 10     1     9 0.00716
#> # ℹ 10,000,090 more rows
predict(model, type='exceedance', threshold = c(10, 100)) # exceedance probs (H.4)
#> # A tibble: 100 × 2
#>    p_ge_10 p_ge_100
#>      <dbl>    <dbl>
#>  1  0.321    0.170 
#>  2  0.330    0.166 
#>  3  0.497    0.365 
#>  4  0.0888   0.0408
#>  5  0.302    0.141 
#>  6  0.290    0.151 
#>  7  0.421    0.276 
#>  8  0.427    0.285 
#>  9  0.289    0.142 
#> 10  0.301    0.139 
#> # ℹ 90 more rows
predict(model, type='draws', n_draws = 100) # predictive draws (H.5)
#> # A tibble: 10,000 × 3
#>     .row .draw     y
#>    <int> <int> <dbl>
#>  1     1     1     1
#>  2     2     1   150
#>  3     3     1     0
#>  4     4     1     0
#>  5     5     1     2
#>  6     6     1     0
#>  7     7     1     0
#>  8     8     1     4
#>  9     9     1     0
#> 10    10     1     1
#> # ℹ 9,990 more rows
# }
```
