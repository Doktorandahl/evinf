# EVZINB tidy function

EVZINB tidy function

## Usage

``` r
# S3 method for class 'evzinb'
tidy(
  x,
  component = c("all", "count", "zero", "evi", "pareto"),
  coef_type = c("original", "bootstrap_mean", "bootstrap_median"),
  standard_error = TRUE,
  p_value = c("bootstrapped", "approx", "none"),
  confint = c("none", "bootstrapped", "approx"),
  conf_level = 0.95,
  approx_t_value = TRUE,
  symmetric_bootstrap_p = TRUE,
  exclude_degenerate = TRUE,
  ...
)
```

## Arguments

- x:

  An evzinb object

- component:

  Which component should be shown? One of `"count"`, `"zi"`, `"evi"`,
  `"pareto"` or `"all"` (the default).

- coef_type:

  Type of coefficients. Original are the coefficient estimates from the
  non-bootstrapped version of the model. 'bootstrap_mean' are the mean
  coefficients across bootstraps, and 'bootstrap_median' are the median
  coefficients across bootstraps

- standard_error:

  Should bootstrapped standard errors be computed?

- p_value:

  What type of p_values should be computed? 'bootstrapped' are
  bootstrapped p_values through confidence interval inversion. 'approx'
  are p-values based on the t-value produced by dividing the coefficient
  with the standard error.

- confint:

  What type of confidence interval should be computed: `"none"`,
  `"bootstrapped"` (percentile interval from the bootstrap distribution)
  or `"approx"` (`estimate +/- qt() * std.error`).

- conf_level:

  Confidence level for the confidence interval

- approx_t_value:

  Should approximate t-values be returned

- symmetric_bootstrap_p:

  Should bootstrap p-values be computed as symmetric (leaving alpha/2
  percent in each tail)? FALSE gives non-symmetric, but narrower,
  intervals. TRUE corresponds most closely to conventional p-values.

- exclude_degenerate:

  Drop bootstrap replicates flagged degenerate (default TRUE); see the
  alpha_floor argument of evinf_control().

- ...:

  Other arguments passed to the tidy function

## Value

A tibble with one row per coefficient

## Details

When `x` was fitted with `bootstrap = FALSE`, standard errors, p-values,
t-values and confidence intervals are unavailable; the estimate column
is returned on its own and a message is emitted. Requesting bootstrapped
coefficients for such a model is an error.

A bootstrapped `p.value` is never reported below `1 / B`, where `B` is
the number of usable bootstrap replicates: with no draw crossing the
estimate, the true p-value could be anywhere in `[0, 1/B)`, so it is
floored at `1/B` rather than reported as exactly `0`.

With `component = "all"` the output has one row per coefficient *per
component* (a `y.level` column names the component). Like
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html), this
multi-component shape cannot be rendered by
[`modelsummary::modelsummary()`](https://modelsummary.com/man/modelsummary.html)
with its default arguments; pass `shape = term + y.level ~ model` so
each component becomes its own block of rows.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
tidy(model)
#> # A tibble: 16 × 6
#>    y.level term        estimate std.error statistic p.value
#>    <fct>   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
#>  1 zero    (Intercept)    0.511        NA        NA       1
#>  2 zero    x1            -0.674        NA        NA       1
#>  3 zero    x2             0.607        NA        NA       1
#>  4 zero    x3            -0.420        NA        NA       1
#>  5 evi     (Intercept)   -1.73         NA        NA       1
#>  6 evi     x1             0.736        NA        NA       1
#>  7 evi     x2            -0.429        NA        NA       1
#>  8 evi     x3             0.166        NA        NA       1
#>  9 count   (Intercept)    3.13         NA        NA       1
#> 10 count   x1             0.843        NA        NA       1
#> 11 count   x2             0.599        NA        NA       1
#> 12 count   x3             0.241        NA        NA       1
#> 13 pareto  (Intercept)    2.71         NA        NA       1
#> 14 pareto  x1            -2.28         NA        NA       1
#> 15 pareto  x2             1.45         NA        NA       1
#> 16 pareto  x3             1.50         NA        NA       1

# multi-component table with modelsummary
if (requireNamespace("modelsummary", quietly = TRUE) &&
    requireNamespace("broom", quietly = TRUE)) {
  modelsummary::modelsummary(
    model,
    shape = term + y.level ~ model,
    gof_map = gof_map_evinf()
  )
}
#> 
#> +-----------------------+---------+---------+
#> |                       | y.level | (1)     |
#> +=======================+=========+=========+
#> | (Intercept)           | zero    | 0.511   |
#> +-----------------------+---------+---------+
#> |                       | evi     | -1.729  |
#> +-----------------------+---------+---------+
#> |                       | count   | 3.131   |
#> +-----------------------+---------+---------+
#> |                       | pareto  | 2.714   |
#> +-----------------------+---------+---------+
#> | x1                    | zero    | -0.674  |
#> +-----------------------+---------+---------+
#> |                       | evi     | 0.736   |
#> +-----------------------+---------+---------+
#> |                       | count   | 0.843   |
#> +-----------------------+---------+---------+
#> |                       | pareto  | -2.278  |
#> +-----------------------+---------+---------+
#> | x2                    | zero    | 0.607   |
#> +-----------------------+---------+---------+
#> |                       | evi     | -0.429  |
#> +-----------------------+---------+---------+
#> |                       | count   | 0.599   |
#> +-----------------------+---------+---------+
#> |                       | pareto  | 1.449   |
#> +-----------------------+---------+---------+
#> | x3                    | zero    | -0.420  |
#> +-----------------------+---------+---------+
#> |                       | evi     | 0.166   |
#> +-----------------------+---------+---------+
#> |                       | count   | 0.241   |
#> +-----------------------+---------+---------+
#> |                       | pareto  | 1.503   |
#> +-----------------------+---------+---------+
#> | Observations          |         | 100     |
#> +-----------------------+---------+---------+
#> | Parameters            |         | 18      |
#> +-----------------------+---------+---------+
#> | alpha_nb              |         | 1.43    |
#> +-----------------------+---------+---------+
#> | C_EV                  |         | 184     |
#> +-----------------------+---------+---------+
#> | Obs. above C_EV       |         | 10      |
#> +-----------------------+---------+---------+
#> | logLik                |         | -254.03 |
#> +-----------------------+---------+---------+
#> | AIC                   |         | 544.1   |
#> +-----------------------+---------+---------+
#> | BIC                   |         | 591.0   |
#> +-----------------------+---------+---------+
#> | Usable bootstraps     |         | 1       |
#> +-----------------------+---------+---------+
#> | Failed bootstraps     |         | 0       |
#> +-----------------------+---------+---------+
#> | Degenerate bootstraps |         | 4       |
#> +-----------------------+---------+---------+
#> | C_EV on boundary      |         | 0       |
#> +-----------------------+---------+---------+
#> | Converged             |         | TRUE    |
#> +-----------------------+---------+---------+
#> | C_EV profile settled  |         | TRUE    |
#> +-----------------------+---------+---------+ 
# }

# \donttest{
data(hks)
hks_mod <- evzinb(osvAll ~ troopLag + lntpop + brv_AllLag_log,
                  data = hks, n_bootstraps = 2, multicore = FALSE)
#> evinf: using a data-driven candidate range for C_EV: [209, 8586]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (209) in 1 of 2 bootstrap replicates; consider widening c.lim.
tidy(hks_mod)
#> # A tibble: 16 × 6
#>    y.level term           estimate std.error statistic p.value
#>    <fct>   <chr>             <dbl>     <dbl>     <dbl>   <dbl>
#>  1 zero    (Intercept)     6.47      0.119      54.5       0.5
#>  2 zero    troopLag       -0.0491    0.0381     -1.29      0.5
#>  3 zero    lntpop         -0.516     0.0169    -30.5       0.5
#>  4 zero    brv_AllLag_log -0.698     0.00641  -109.        0.5
#>  5 evi     (Intercept)    -8.61      0.0209   -412.        0.5
#>  6 evi     troopLag       -0.0102    0.0165     -0.619     1  
#>  7 evi     lntpop          0.568     0.0141     40.4       0.5
#>  8 evi     brv_AllLag_log  0.115     0.0130      8.83      0.5
#>  9 count   (Intercept)     5.44      0.434      12.5       0.5
#> 10 count   troopLag       -0.00568   0.0132     -0.431     1  
#> 11 count   lntpop         -0.220     0.0376     -5.86      0.5
#> 12 count   brv_AllLag_log  0.0623    0.0472      1.32      0.5
#> 13 pareto  (Intercept)    -1.85      1.70       -1.09      1  
#> 14 pareto  troopLag       -0.0480    0.0875     -0.549     1  
#> 15 pareto  lntpop          0.163     0.131       1.25      1  
#> 16 pareto  brv_AllLag_log  0.0304    0.0792      0.384     1  
# }
```
