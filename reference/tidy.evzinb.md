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
#> Error in eval(expr, p) : inv(): matrix is singular
tidy(model)
#> # A tibble: 16 × 6
#>    y.level term        estimate std.error statistic p.value
#>    <fct>   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
#>  1 zero    (Intercept)    0.511    0.204      2.50        0
#>  2 zero    x1            -0.674    0.204     -3.31        0
#>  3 zero    x2             0.607    0.323      1.88        0
#>  4 zero    x3            -0.420    0.0187   -22.5         0
#>  5 evi     (Intercept)   -1.73     0.566     -3.06        0
#>  6 evi     x1             0.736    0.210      3.50        1
#>  7 evi     x2            -0.429    0.792     -0.542       1
#>  8 evi     x3             0.166    0.0134    12.4         0
#>  9 count   (Intercept)    3.13     0.269     11.6         0
#> 10 count   x1             0.843    0.152      5.54        0
#> 11 count   x2             0.599    0.172      3.47        0
#> 12 count   x3             0.241    0.528      0.457       1
#> 13 pareto  (Intercept)    2.71     0.302      8.99        0
#> 14 pareto  x1            -2.28     1.51      -1.50        0
#> 15 pareto  x2             1.45     1.02       1.42        0
#> 16 pareto  x3             1.50     0.772      1.95        0

# multi-component table with modelsummary
if (requireNamespace("modelsummary", quietly = TRUE)) {
  modelsummary::modelsummary(
    model,
    shape = term + y.level ~ model,
    gof_map = gof_map_evinf()
  )
}
#> Error: Package `broom` required for this function to work.
#>   Please install it by running `install.packages("broom")`.
# }

if (FALSE) { # \dontrun{
data(hks)
hks_mod <- evzinb(osvAll ~ troopLag + lntpop + brv_AllLag_log,
                  data = hks, n_bootstraps = 5, multicore = FALSE)
tidy(hks_mod)
} # }
```
