# Augment data with fitted values, residuals and state classification

Augment data with fitted values, residuals and state classification

## Usage

``` r
# S3 method for class 'evinb'
augment(x, data = NULL, newdata = NULL, seed = NULL, ...)
```

## Arguments

- x:

  A fitted `evinb` model.

- data, newdata:

  Data to augment. `newdata` takes precedence over `data` when both are
  given (the broom convention); both default to `NULL`, in which case
  the model's own estimation data is used.

- seed:

  Optional seed for the randomized quantile residual (`.resid`); see
  `residuals(type = "quantile")`.

- ...:

  Unused.

## Value

See [`augment.evzinb`](augment.evzinb.md); `evinb` has no zero state, so
`.prob_zero`/`.post_zero` are not produced.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evinb(y~x1+x2+x3, data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
augment(model)
#> # A tibble: 100 × 12
#>        y      x1     x2      x3 .fitted .prob_count .prob_evi .state .resid
#>    <dbl>   <dbl>  <dbl>   <dbl>   <dbl>       <dbl>     <dbl> <fct>   <dbl>
#>  1     2  0.702  -0.134 -0.253     53.5       0.888    0.112  count   0.204
#>  2     0  0.654   2.29   1.03      60.1       0.955    0.0447 count  -1.85 
#>  3     0  1.94    2.81   1.56     410.        0.883    0.117  count  -1.04 
#>  4     0 -0.772  -1.15  -1.25      12.4       0.959    0.0408 count   0.540
#>  5   184  0.332   1.23   1.13      39.2       0.937    0.0628 count   1.34 
#>  6    23  0.0440 -0.625  0.861     34.1       0.879    0.121  count   0.734
#>  7     0  1.62    0.813 -0.453    149.        0.865    0.135  count  -0.457
#>  8   207  1.25    0.304  0.372     99.3       0.833    0.167  count   0.837
#>  9     2  0.366  -0.324  0.0459    36.9       0.896    0.104  count   0.321
#> 10     0  0.787   1.03  -0.425     40.7       0.939    0.0612 count  -0.266
#> # ℹ 90 more rows
#> # ℹ 3 more variables: .post_count <dbl>, .post_evi <dbl>, .post_state <fct>
# }
```
