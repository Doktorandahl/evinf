# Log-likelihood profile over the candidate values of C_EV

Returns the log-likelihood evaluated at each candidate value of
\\C\_{EV}\\ from the final EM iteration. Use it to judge whether
\\C\_{EV}\\ is well identified and whether the candidate range (`c.lim`)
is wide enough.

## Usage

``` r
c_profile(object)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` object.

## Value

A tibble with columns `c` and `loglik`.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
c_profile(model)
#> # A tibble: 9 × 2
#>       c loglik
#>   <dbl>  <dbl>
#> 1   173  -259.
#> 2   184  -254.
#> 3   190  -257.
#> 4   191  -260.
#> 5   207  -262.
#> 6   210  -263.
#> 7   216  -266.
#> 8   237  -318.
#> 9   263  -321.
# }
```
