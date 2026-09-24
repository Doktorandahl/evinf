# Random draws from a fitted evinb model

**Superseded:** kept for backwards compatibility, with no plans for
removal, but new code should use [`simulate()`](evinf-s3-predict.md),
which returns a tidy data frame (and a reproducible `"seed"` attribute)
instead of a bare vector/list.

## Usage

``` r
revinb_fit(object, newdata = NULL, n_draws = 1)
```

## Arguments

- object:

  A fitted EVINB object

- newdata:

  Optional newdata

- n_draws:

  Number of random draws to make

## Value

A vector of randomly drawn values from the fitted evinb if n_draws == 1,
or a list of length n_draws with random drawn values if n_draws \> 1

## See also

[`simulate.evinb`](evinf-s3-predict.md)

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evinb(y~x1+x2+x3, data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 2 of 5 bootstrap replicates; consider widening c.lim.
revinb_fit(model)
#>   [1] 1.000000e+00 0.000000e+00 1.860000e+02 0.000000e+00 1.000000e+00
#>   [6] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.880000e+02
#>  [11] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [16] 4.000000e+00 1.702000e+03 2.000000e+00 2.380000e+02 0.000000e+00
#>  [21] 0.000000e+00 0.000000e+00 2.000000e+00 0.000000e+00 0.000000e+00
#>  [26] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#>  [31] 0.000000e+00 0.000000e+00 8.400000e+01 1.840000e+02 0.000000e+00
#>  [36] 1.000000e+00 0.000000e+00 0.000000e+00 3.000000e+00 1.100000e+01
#>  [41] 0.000000e+00 3.000000e+00 0.000000e+00 5.000000e+01 0.000000e+00
#>  [46] 1.860000e+02 0.000000e+00 0.000000e+00 6.000000e+01 0.000000e+00
#>  [51] 0.000000e+00 0.000000e+00 0.000000e+00 2.800000e+02 0.000000e+00
#>  [56] 0.000000e+00 0.000000e+00 2.110000e+02 3.580000e+02 0.000000e+00
#>  [61] 2.000000e+01 0.000000e+00 0.000000e+00 1.900000e+01 1.810000e+02
#>  [66] 0.000000e+00 0.000000e+00 1.880000e+02 0.000000e+00 3.370000e+02
#>  [71] 0.000000e+00 0.000000e+00 0.000000e+00 2.025000e+03 0.000000e+00
#>  [76] 0.000000e+00 0.000000e+00 2.000000e+00 1.850000e+02 0.000000e+00
#>  [81] 3.900000e+01 1.840000e+02 0.000000e+00 0.000000e+00 1.900000e+01
#>  [86] 0.000000e+00 0.000000e+00 2.190000e+02 0.000000e+00 1.500000e+01
#>  [91] 0.000000e+00 1.840000e+02 3.129716e+25 0.000000e+00 0.000000e+00
#>  [96] 0.000000e+00 3.916000e+03 0.000000e+00 0.000000e+00 3.000000e+00
# }
```
