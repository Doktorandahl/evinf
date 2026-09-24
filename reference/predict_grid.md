# Build a prediction grid and evaluate a model over it

Varies one covariate over a range (or set of values) while holding the
others at representative values, and returns the model predictions in
long form - the raw material for an effect plot.

## Usage

``` r
predict_grid(
  object,
  variable,
  values = NULL,
  n = 50,
  at = list(),
  fixed = c("mean", "median"),
  type = c("states", "harmonic", "explog", "quantile", "counts", "pareto_alpha"),
  quantile = NULL,
  confint = FALSE,
  conf_level = 0.9,
  clamp_alpha_pl = FALSE
)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` model.

- variable:

  Name of the covariate to vary.

- values:

  Optional explicit values for `variable`; otherwise an evenly spaced
  sequence (numeric) or all levels (factor).

- n:

  Number of grid points for a numeric `variable`.

- at:

  Named list of values at which to hold specific other covariates.

- fixed:

  How to hold the remaining numeric covariates: `"mean"` or `"median"`.

- type:

  Prediction type: `"states"`, `"harmonic"`, `"explog"`, `"quantile"`,
  `"counts"` or `"pareto_alpha"`.

- quantile:

  Quantile for `type = "quantile"`.

- confint:

  Add bootstrap confidence intervals (not for `"states"`).

- conf_level:

  Confidence level.

- clamp_alpha_pl:

  (round11 A4) `type = "explog"` only; forwarded to
  [`predict()`](https://rdrr.io/r/stats/predict.html). See
  [`?predict.evzinb`](predict.evzinb.md).

## Value

A tibble with columns `variable`, `value`, `type`, `estimate`,
`conf.low`, `conf.high` (the last two `NA` when `confint = FALSE`).

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
predict_grid(model, "x1", type = "harmonic")
#> # A tibble: 50 × 6
#>    variable value type     estimate conf.low conf.high
#>    <chr>    <dbl> <chr>       <dbl>    <dbl>     <dbl>
#>  1 x1       -3.01 harmonic    0.385       NA        NA
#>  2 x1       -2.91 harmonic    0.443       NA        NA
#>  3 x1       -2.81 harmonic    0.510       NA        NA
#>  4 x1       -2.71 harmonic    0.587       NA        NA
#>  5 x1       -2.61 harmonic    0.675       NA        NA
#>  6 x1       -2.51 harmonic    0.777       NA        NA
#>  7 x1       -2.41 harmonic    0.893       NA        NA
#>  8 x1       -2.31 harmonic    1.03        NA        NA
#>  9 x1       -2.21 harmonic    1.18        NA        NA
#> 10 x1       -2.10 harmonic    1.35        NA        NA
#> # ℹ 40 more rows
# }

# \donttest{
data(hks)
hks_mod <- evzinb(osvAll ~ troopLag + lntpop + brv_AllLag_log,
                  data = hks, n_bootstraps = 2, multicore = FALSE)
#> evinf: using a data-driven candidate range for C_EV: [209, 8586]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
predict_grid(hks_mod, "troopLag", type = "states")
#> # A tibble: 150 × 6
#>    variable value type     estimate conf.low conf.high
#>    <chr>    <dbl> <chr>       <dbl>    <dbl>     <dbl>
#>  1 troopLag 0     pr_zero    0.742        NA        NA
#>  2 troopLag 0     pr_count   0.248        NA        NA
#>  3 troopLag 0     pr_evi     0.0100       NA        NA
#>  4 troopLag 0.596 pr_zero    0.736        NA        NA
#>  5 troopLag 0.596 pr_count   0.254        NA        NA
#>  6 troopLag 0.596 pr_evi     0.0102       NA        NA
#>  7 troopLag 1.19  pr_zero    0.730        NA        NA
#>  8 troopLag 1.19  pr_count   0.259        NA        NA
#>  9 troopLag 1.19  pr_evi     0.0103       NA        NA
#> 10 troopLag 1.79  pr_zero    0.725        NA        NA
#> # ℹ 140 more rows
# }
```
