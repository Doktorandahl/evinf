# EVINB print function

EVINB print function

## Usage

``` r
# S3 method for class 'evinb'
print(x, ...)
```

## Arguments

- x:

  A fitted evinb model

- ...:

  Not used

## Value

`x`, invisibly.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
#> Error in eval(expr, p) : inv(): matrix is singular
print(model)
#> 
#>  Fitted EVINB model with formulas: 
#>  NB:     y ~ x1 + x2 + x3
#>  EVI:    y ~ x1 + x2 + x3 
#>  Pareto: y ~ x1 + x2 + x3 
#>  ______ 
#>  Converged:                      TRUE 
#>  C_EV:                           184 
#>  Candidate range for C_EV:       [173, 263]  (data-driven) 
#>  Observations at or above C_EV:  10 
#>  Parameters:                     14 
#>  Bootstraps (failed, degenerate): 3 (2, 0)
#>    (call failed_bootstraps() for the details)
# }
```
