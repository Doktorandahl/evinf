# Augment data with fitted values, residuals and state classification

Augment data with fitted values, residuals and state classification

## Usage

``` r
# S3 method for class 'evzinb'
augment(x, data = NULL, newdata = NULL, seed = NULL, ...)
```

## Arguments

- x:

  A fitted `evzinb` model.

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

`data`/`newdata` (or the estimation data) with columns appended:
`.fitted` (the `predict(type = "harmonic")` point prediction),
`.prob_zero`/`.prob_count`/`.prob_evi` (prior state probabilities –
[`augment.evinb()`](augment.evinb.md) has no zero state, so `.prob_zero`
is dropped there) and `.state` (the MAP prior state, from
[`classify_states()`](classify_states.md)). When the response column is
present, three more columns are added: `.resid` (a seeded randomized
quantile residual, `residuals(type = "quantile")`),
`.post_zero`/`.post_count`/`.post_evi` (posterior state probabilities)
and `.post_state` (the MAP posterior state).

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3, data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
augment(model)
#> # A tibble: 100 × 14
#>        y      x1     x2      x3 .fitted .prob_zero .prob_count .prob_evi .state
#>    <dbl>   <dbl>  <dbl>   <dbl>   <dbl>      <dbl>       <dbl>     <dbl> <fct> 
#>  1     2  0.702  -0.134 -0.253     52.4      0.450       0.422    0.128  zero  
#>  2     0  0.654   2.29   1.03      57.2      0.712       0.255    0.0326 zero  
#>  3     0  1.94    2.81   1.56     377.       0.500       0.388    0.112  zero  
#>  4     0 -0.772  -1.15  -1.25      11.2      0.674       0.287    0.0385 zero  
#>  5   184  0.332   1.23   1.13      38.8      0.600       0.344    0.0556 zero  
#>  6    23  0.0440 -0.625  0.861     35.9      0.377       0.488    0.135  count 
#>  7     0  1.62    0.813 -0.453    126.       0.445       0.401    0.154  zero  
#>  8   207  1.25    0.304  0.372     90.8      0.342       0.464    0.193  count 
#>  9     2  0.366  -0.324  0.0459    37.4      0.453       0.431    0.116  zero  
#> 10     0  0.787   1.03  -0.425     34.0      0.649       0.295    0.0559 zero  
#> # ℹ 90 more rows
#> # ℹ 5 more variables: .resid <dbl>, .post_zero <dbl>, .post_count <dbl>,
#> #   .post_evi <dbl>, .post_state <fct>
# }
```
