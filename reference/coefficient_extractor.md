# Bootstrap coefficient extractor

Bootstrap coefficient extractor

## Usage

``` r
coefficient_extractor(object, ...)

# S3 method for class 'evzinb'
coefficient_extractor(
  object,
  component = c("all", "count", "zero", "evi", "pareto"),
  exclude_degenerate = TRUE,
  ...
)

# S3 method for class 'evinb'
coefficient_extractor(
  object,
  component = c("all", "count", "evi", "pareto"),
  exclude_degenerate = TRUE,
  ...
)

# S3 method for class 'zinbboot'
coefficient_extractor(object, component = c("all", "count", "zero"), ...)

# S3 method for class 'nbboot'
coefficient_extractor(object, ...)
```

## Arguments

- object:

  a fitted model with bootstraps of class evzinb, evinb, nbboot, or
  zinbboot

- ...:

  Arguments passed to methods, in particular `component` (not for
  nbboot): one of `"count"`, `"zero"`, `"evi"`, `"pareto"` or `"all"`.
  The pre-0.9.4 names `"nb"`, `"zi"` and `"evinf"` are still accepted
  with a deprecation warning.

- component:

  Which component should be extracted

- exclude_degenerate:

  Drop bootstrap replicates flagged degenerate (default `TRUE`); see
  [`evinf_control`](evinf_control.md).

## Value

A tibble with coefficient values, one row per bootstrap and component

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
coefficient_extractor(model, component = 'all')
#> # A tibble: 12 × 5
#>    `(Intercept)`      x1     x2      x3 .component
#>            <dbl>   <dbl>  <dbl>   <dbl> <chr>     
#>  1         2.74   0.571   0.518  1.16   count     
#>  2         3.62   0.0664  0.899  0.104  count     
#>  3         3.02   0.650   0.564 -0.245  count     
#>  4         0.390 -0.410   1.19  -0.769  zero      
#>  5         1.10  -1.07    0.725 -1.03   zero      
#>  6         0.902 -0.272   0.490 -0.148  zero      
#>  7        -1.40   0.501   0.137 -0.187  evi       
#>  8        -0.732  0.146   0.238  0.0272 evi       
#>  9        -1.32   1.47   -1.06   0.676  evi       
#> 10         2.71  -2.26    1.44   1.50   pareto    
#> 11         2.81  -2.66    1.69   1.74   pareto    
#> 12         1.98  -1.44    0.698  1.37   pareto    
# }
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
zinb_comp <- compare_models(model)
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: alternation limit reached
#> Warning: alternation limit reached
coefficient_extractor(zinb_comp$zinb)
#> # A tibble: 10 × 5
#>    `(Intercept)`     x1      x2      x3 .component
#>            <dbl>  <dbl>   <dbl>   <dbl> <chr>     
#>  1        3.23    1.08   0.318   0.183  count     
#>  2        4.07    1.35   0.212  -0.903  count     
#>  3        4.15    0.626 -0.0286  0.0312 count     
#>  4        3.82    1.02  -0.537   0.351  count     
#>  5        2.81    1.58   0.594   0.601  count     
#>  6        0.539  -0.779  0.772  -0.240  zero      
#>  7        0.499  -0.636  0.160  -0.255  zero      
#>  8       -0.0832 -1.35   0.215   0.401  zero      
#>  9        0.613  -1.65   1.08   -1.35   zero      
#> 10       -0.131  -0.894  1.51   -1.14   zero      
# }
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
zinb_comp <- compare_models(model)
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: alternation limit reached
coefficient_extractor(zinb_comp$nb)
#> # A tibble: 5 × 4
#>   `(Intercept)`    x1      x2     x3
#>           <dbl> <dbl>   <dbl>  <dbl>
#> 1          2.75 1.76  -1.18   0.775 
#> 2          3.36 1.32  -0.636  0.944 
#> 3          3.14 1.75  -0.915  0.0786
#> 4          2.67 1.41   0.0792 0.332 
#> 5          3.35 0.892 -0.139  0.450 
# }
```
