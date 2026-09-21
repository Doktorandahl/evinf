# EVINB tidy function

EVINB tidy function

## Usage

``` r
# S3 method for class 'evinb'
tidy(
  x,
  component = c("all", "count", "evi", "pareto"),
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

  An evinb object

- component:

  Which component should be shown? One of `"count"`, `"evi"`, `"pareto"`
  or `"all"` (the default).

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
model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
tidy(model)
#> # A tibble: 12 × 6
#>    y.level term        estimate std.error statistic p.value
#>    <fct>   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
#>  1 evi     (Intercept)   -2.69    0.183     -14.7       0.5
#>  2 evi     x1             0.902   0.213       4.24      0.5
#>  3 evi     x2            -0.580   0.337      -1.72      0.5
#>  4 evi     x3             0.354   0.0746      4.74      0.5
#>  5 count   (Intercept)    2.25    0.197      11.4       0.5
#>  6 count   x1             1.39    0.129      10.8       0.5
#>  7 count   x2             0.192   0.0206      9.32      0.5
#>  8 count   x3             0.386   0.00275   141.        0.5
#>  9 pareto  (Intercept)    2.80    1.46        1.92      0.5
#> 10 pareto  x1            -2.45    0.179     -13.7       0.5
#> 11 pareto  x2             1.57    0.833       1.89      0.5
#> 12 pareto  x3             1.68    1.82        0.920     1  
# }
```
