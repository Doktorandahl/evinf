# Cross-tabulate prior vs posterior state classification

Cross-tabulate prior vs posterior state classification

## Usage

``` r
state_table(object)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` model.

## Value

An object of class `evinf_state_table`: the `map_prior` x
`map_posterior` contingency table with row percentages.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
state_table(model)
#> Prior x posterior state classification (counts)
#>        posterior
#> prior   zero count evi
#>   zero    51    18   3
#>   count   10    12   6
#>   evi      0     0   0
#> 
#> Row percentages
#>        posterior
#> prior   zero count  evi
#>   zero  70.8  25.0  4.2
#>   count 35.7  42.9 21.4
#>   evi                  
# }
```
