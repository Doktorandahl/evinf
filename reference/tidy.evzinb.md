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
#>  1 zero    (Intercept)    0.511    0.245      2.09    0.333
#>  2 zero    x1            -0.674    0.428     -1.57    0.333
#>  3 zero    x2             0.607    0.0781     7.77    0.333
#>  4 zero    x3            -0.420    0.299     -1.40    1    
#>  5 evi     (Intercept)   -1.73     0.285     -6.06    0.333
#>  6 evi     x1             0.736    0.472      1.56    0.333
#>  7 evi     x2            -0.429    0.332     -1.29    0.333
#>  8 evi     x3             0.166    0.267      0.621   0.333
#>  9 count   (Intercept)    3.13     0.354      8.83    0.333
#> 10 count   x1             0.843    0.557      1.51    0.667
#> 11 count   x2             0.599    0.487      1.23    0.333
#> 12 count   x3             0.241    0.442      0.545   0.333
#> 13 pareto  (Intercept)    2.71     0.105     25.8     0.333
#> 14 pareto  x1            -2.28     0.103    -22.1     0.333
#> 15 pareto  x2             1.45     0.0975    14.9     0.333
#> 16 pareto  x3             1.50     0.165      9.13    0.333

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
#> |                       |         | (0.245) |
#> +-----------------------+---------+---------+
#> |                       | evi     | -1.729  |
#> +-----------------------+---------+---------+
#> |                       |         | (0.285) |
#> +-----------------------+---------+---------+
#> |                       | count   | 3.131   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.354) |
#> +-----------------------+---------+---------+
#> |                       | pareto  | 2.714   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.105) |
#> +-----------------------+---------+---------+
#> | x1                    | zero    | -0.674  |
#> +-----------------------+---------+---------+
#> |                       |         | (0.428) |
#> +-----------------------+---------+---------+
#> |                       | evi     | 0.736   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.472) |
#> +-----------------------+---------+---------+
#> |                       | count   | 0.843   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.557) |
#> +-----------------------+---------+---------+
#> |                       | pareto  | -2.278  |
#> +-----------------------+---------+---------+
#> |                       |         | (0.103) |
#> +-----------------------+---------+---------+
#> | x2                    | zero    | 0.607   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.078) |
#> +-----------------------+---------+---------+
#> |                       | evi     | -0.429  |
#> +-----------------------+---------+---------+
#> |                       |         | (0.332) |
#> +-----------------------+---------+---------+
#> |                       | count   | 0.599   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.487) |
#> +-----------------------+---------+---------+
#> |                       | pareto  | 1.449   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.098) |
#> +-----------------------+---------+---------+
#> | x3                    | zero    | -0.420  |
#> +-----------------------+---------+---------+
#> |                       |         | (0.299) |
#> +-----------------------+---------+---------+
#> |                       | evi     | 0.166   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.267) |
#> +-----------------------+---------+---------+
#> |                       | count   | 0.241   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.442) |
#> +-----------------------+---------+---------+
#> |                       | pareto  | 1.503   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.165) |
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
#> | Usable bootstraps     |         | 3       |
#> +-----------------------+---------+---------+
#> | Failed bootstraps     |         | 0       |
#> +-----------------------+---------+---------+
#> | Degenerate bootstraps |         | 2       |
#> +-----------------------+---------+---------+
#> | C_EV on boundary      |         | 0       |
#> +-----------------------+---------+---------+
#> | Converged             |         | TRUE    |
#> +-----------------------+---------+---------+
#> | C_EV profile settled  |         | TRUE    |
#> +-----------------------+---------+---------+ 
# }

if (FALSE) { # \dontrun{
data(hks)
hks_mod <- evzinb(osvAll ~ troopLag + lntpop + brv_AllLag_log,
                  data = hks, n_bootstraps = 5, multicore = FALSE)
tidy(hks_mod)
} # }
```
