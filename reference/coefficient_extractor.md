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
#> # A tibble: 16 × 5
#>    `(Intercept)`     x1     x2      x3 .component
#>            <dbl>  <dbl>  <dbl>   <dbl> <chr>     
#>  1        3.11    0.604 0.548   0.882  count     
#>  2        3.29    0.184 1.06    0.0132 count     
#>  3        2.92    0.947 0.499   0.433  count     
#>  4        3.07    0.722 0.355   0.448  count     
#>  5        1.09   -1.46  1.17   -1.01   zero      
#>  6        0.0606 -0.868 0.934  -0.0557 zero      
#>  7        0.233  -1.12  0.659  -0.525  zero      
#>  8        0.695  -0.752 0.552  -0.640  zero      
#>  9       -1.31   -0.158 0.943  -0.660  evi       
#> 10       -1.48    0.217 0.316   0.443  evi       
#> 11       -1.44    0.451 0.0772  0.272  evi       
#> 12       -1.28    1.11  0.261  -0.274  evi       
#> 13        3.30   -3.09  1.12    2.54   pareto    
#> 14        2.72   -2.40  1.51    1.62   pareto    
#> 15        3.00   -3.43  2.22    2.46   pareto    
#> 16        3.18   -2.94  1.22    2.34   pareto    
# }
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
#> Error in eval(expr, p) : inv(): matrix is singular
zinb_comp <- compare_models(model)
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Error in -bootstrap$boot_id : invalid argument to unary operator
coefficient_extractor(zinb_comp$zinb)
#> # A tibble: 6 × 5
#>   `(Intercept)`      x1     x2     x3 .component
#>           <dbl>   <dbl>  <dbl>  <dbl> <chr>     
#> 1         3.91   0.979   0.559 -0.341 count     
#> 2         3.94   0.943  -0.433  0.125 count     
#> 3         3.54  -0.0195  0.310  0.380 count     
#> 4         0.176 -0.988   1.28  -0.815 zero      
#> 5         0.284 -1.91    1.67  -1.70  zero      
#> 6        -1.03  -1.22    1.38   0.310 zero      
# }
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
#> Error in eval(expr, p) : inv(): matrix is singular
zinb_comp <- compare_models(model)
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: alternation limit reached
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
coefficient_extractor(zinb_comp$nb)
#> # A tibble: 3 × 4
#>   `(Intercept)`    x1     x2    x3
#>           <dbl> <dbl>  <dbl> <dbl>
#> 1          2.97  1.37 -0.304 1.19 
#> 2          2.95  1.26 -0.113 0.199
#> 3          3.27  1.01 -0.178 0.288
# }
```
