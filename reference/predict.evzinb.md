# Predictions from evzinb object

Predictions from evzinb object

## Usage

``` r
# S3 method for class 'evzinb'
predict(
  object,
  newdata = NULL,
  type = c("harmonic", "explog", "counts", "pareto_alpha", "zi", "evinf", "count_state",
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
  probability of the zero state, 'all' for all predicted values,
  'quantile' for quantile prediction (scalar or vector, see
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
'pareto_alpha','zi','evinf', 'count_state', and 'quantile' (scalar), or
a tibble for type 'states', 'all', 'quantile' (vector), 'distribution'
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
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
predict(model)
#>   [1]   52.3536771   57.1668381  376.9452691   11.1639895   38.8170530
#>   [6]   35.9055805  126.3183703   90.7615231   37.3853032   34.0306144
#>  [11]    3.8757346    0.7290348  211.7443385   12.7340813   32.1514187
#>  [16]   14.0972436  391.8414865   25.8481011  220.0886401    3.4952654
#>  [21]    9.8384299    5.9820263   12.6066990   12.1342843    0.7762911
#>  [26]   76.3148528   14.1944891   18.5612434    9.6734937   28.8526824
#>  [31]   16.5098346   22.4290440   38.6944355   41.3541611    7.8551612
#>  [36]  114.3658305    3.1848768   22.0641626   10.4536678   71.7961911
#>  [41] 1339.4850102    3.5523495   10.0391498   19.7056931   28.4332808
#>  [46]  102.5571591   57.4677193   38.7979658  100.6043583    2.7228403
#>  [51]   20.6829272   41.1815093   15.9032426   55.6687631   95.2497311
#>  [56]   11.9472584    1.1241349  164.7064467   89.6777185   21.1429786
#>  [61]   18.3951387   52.2695259    3.2066429  119.4794325   66.6140902
#>  [66]   26.6852148   10.1881159   29.5384884  188.2644686   68.4199361
#>  [71]    1.8008339   20.0119216   45.7608922  122.9404737    2.6452961
#>  [76]   17.6827618   24.1374749    4.8513851   12.8528147   53.8795843
#>  [81]   33.1428460   14.5640965    7.6367303   16.1336634    7.4460041
#>  [86]   26.0772201    2.8713642  113.1886226    2.0980822   24.7794615
#>  [91]   11.0102491   64.8139949  752.3961864   74.1029725   38.7799200
#>  [96]   21.6369011  132.7022841    9.7755182   19.1935263   38.9803138
predict(model, type='all', quantile = 0.9) # all available predicted values
#> Warning: evinf (explog prediction): 3 fitted Pareto alpha values below 0.1 (smallest: 0.019); the geometric-mean prediction C * exp(1 / alpha_pl) is effectively undefined there (e.g. exp(10) ~= 2.2e4) and may be Inf. Consider predict(..., type = "explog", clamp_alpha_pl = TRUE).
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
predict(model, type='quantile', quantile = c(.5, .9, .99)) # several quantiles (H.3)
#> # A tibble: 100 × 3
#>      q50   q90   q99
#>    <dbl> <dbl> <dbl>
#>  1     2   216   812
#>  2     0   184   748
#>  3     0  1227  3927
#>  4     0     7   313
#>  5     0   182   343
#>  6     3   186   209
#>  7     9   486 15120
#>  8    26   278   671
#>  9     1   190   321
#> 10     0   150   317
#> # ℹ 90 more rows
predict(model, type='distribution') # the full predictive distribution (H.2)
#> evinf (predict distribution): the default support's upper bound (88558596588672, the 0.999 mixture quantile of the heaviest-tailed row) exceeds max_support (1e+05); truncating there. Pass `support = ` explicitly for a wider (or narrower) support.
#> # A tibble: 10,000,100 × 3
#>     .row     y    prob
#>    <int> <int>   <dbl>
#>  1     1     0 0.476  
#>  2     1     1 0.0182 
#>  3     1     2 0.0152 
#>  4     1     3 0.0134 
#>  5     1     4 0.0121 
#>  6     1     5 0.0112 
#>  7     1     6 0.0104 
#>  8     1     7 0.00979
#>  9     1     8 0.00925
#> 10     1     9 0.00876
#> # ℹ 10,000,090 more rows
predict(model, type='exceedance', threshold = c(10, 100)) # exceedance probs (H.4)
#> # A tibble: 100 × 2
#>    p_ge_10 p_ge_100
#>      <dbl>    <dbl>
#>  1  0.415    0.163 
#>  2  0.262    0.171 
#>  3  0.486    0.431 
#>  4  0.0793   0.0385
#>  5  0.336    0.155 
#>  6  0.404    0.143 
#>  7  0.500    0.324 
#>  8  0.573    0.333 
#>  9  0.380    0.132 
#> 10  0.292    0.132 
#> # ℹ 90 more rows
predict(model, type='draws', n_draws = 100) # predictive draws (H.5)
#> # A tibble: 10,000 × 3
#>     .row .draw     y
#>    <int> <int> <dbl>
#>  1     1     1     0
#>  2     2     1     0
#>  3     3     1  2576
#>  4     4     1     0
#>  5     5     1   137
#>  6     6     1     0
#>  7     7     1     0
#>  8     8     1   310
#>  9     9     1     0
#> 10    10     1    15
#> # ℹ 9,990 more rows
# }
```
