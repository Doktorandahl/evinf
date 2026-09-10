# Prior and posterior state classification

For each observation, returns the prior state probabilities (from the
inflation components) and, for the estimation data, the posterior state
probabilities (responsibilities), together with the maximum-a-posteriori
state under each.

## Usage

``` r
classify_states(
  object,
  rule = c("map", "threshold"),
  threshold = 0.5,
  newdata = NULL
)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` model.

- rule:

  `"map"` (assign to the most probable state) or `"threshold"` (assign
  only when the top probability exceeds `threshold`, otherwise `NA`).

- threshold:

  Threshold for `rule = "threshold"`.

- newdata:

  Optional data. Prior probabilities are always available; posterior
  probabilities need the response column and are `NA` otherwise.

## Value

A tibble with `prior_zero` / `prior_count` / `prior_evi`, `posterior_*`,
`map_prior`, `map_posterior` (factors with levels `zero`, `count`,
`evi`) and `y`.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
classify_states(model)
#> # A tibble: 100 × 9
#>    prior_zero prior_count prior_evi posterior_zero posterior_count posterior_evi
#>         <dbl>       <dbl>     <dbl>          <dbl>           <dbl>         <dbl>
#>  1      0.451       0.421    0.128           0             1               0    
#>  2      0.712       0.256    0.0325          0.993         0.00670         0    
#>  3      0.499       0.390    0.111           0.995         0.00502         0    
#>  4      0.676       0.285    0.0386          0.906         0.0941          0    
#>  5      0.600       0.344    0.0555          0             0.0104          0.990
#>  6      0.378       0.487    0.135           0             1               0    
#>  7      0.445       0.401    0.153           0.978         0.0223          0    
#>  8      0.343       0.464    0.193           0             0.219           0.781
#>  9      0.454       0.430    0.116           0             1               0    
#> 10      0.649       0.295    0.0559          0.983         0.0167          0    
#> # ℹ 90 more rows
#> # ℹ 3 more variables: map_prior <fct>, map_posterior <fct>, y <dbl>
# }
```
