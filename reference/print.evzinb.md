# EVZINB print function

EVZINB print function

## Usage

``` r
# S3 method for class 'evzinb'
print(x, ...)
```

## Arguments

- x:

  A fitted evzinb model

- ...:

  Not used

## Value

`x`, invisibly.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
print(model)
#> 
#>  Fitted EVZINB model with formulas: 
#>  NB:     y ~ x1 + x2 + x3
#>  ZI:     y ~ x1 + x2 + x3
#>  EVI:    y ~ x1 + x2 + x3 
#>  Pareto: y ~ x1 + x2 + x3 
#>  ______ 
#>  Converged:                      TRUE 
#>  C_EV:                           184 
#>  Candidate range for C_EV:       [173, 263]  (data-driven) 
#>  Observations at or above C_EV:  10 
#>  Parameters:                     18 
#>  Bootstraps (failed, degenerate): 2 (0, 3)
#>    (call failed_bootstraps() for the details)
# }
```
