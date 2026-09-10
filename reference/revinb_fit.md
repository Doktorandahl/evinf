# Random draws from a fitted evinb model

\`revinb_fit()\` is superseded by \[simulate()\]\[simulate.evinb\],
which returns a tidy data frame; it is kept for backwards compatibility.

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

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evinb(y~x1+x2+x3, data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
revinb_fit(model)
#>   [1]    350     13      9      1      0      0    220      0      9      0
#>  [11]      0      0    188      7      0     20      0      0      0      0
#>  [21]     11      0     11      0      0      0      4      4      1      0
#>  [31]      0    125    184      0      0      0      1      0      0    379
#>  [41]      0      1      0      0     18      0      3     79    250      0
#>  [51]      0      0      0      0      2      0      0      0      0      0
#>  [61]      0    220      0      0      0      1      6      0 103194      0
#>  [71]      0      0      2      4      0      0      0      0      2     37
#>  [81]      0    184      0      0      0      0      0   1940      0      1
#>  [91]      0    184      0      0    118      0     33      0      1    181
# }
```
