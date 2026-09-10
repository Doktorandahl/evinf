# Random draws from a fitted evzinb model

\`revzinb_fit()\` is superseded by \[simulate()\]\[simulate.evzinb\],
which returns a tidy data frame; it is kept for backwards compatibility.

## Usage

``` r
revzinb_fit(object, newdata = NULL, n_draws = 1)
```

## Arguments

- object:

  A fitted EVZINB object

- newdata:

  Optional newdata

- n_draws:

  Number of random draws to make

## Value

A vector of randomly drawn values from the fitted evzinb if n_draws ==
1, or a list of length n_draws with random drawn values if n_draws \> 1

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3, data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
revzinb_fit(model)
#>   [1]    10     0  1100     0    45     6     0   106    26   278     0     0
#>  [13]   138     1     0     0     0     0    15     0     5     0     0     0
#>  [25]     3     0     8     2     0     9     0     0     0    34     8    23
#>  [37]     0     3     0    81    26     0     0     0     8    11     0     0
#>  [49]    50     1    23     9     0   146     3     0     0 23061    69    82
#>  [61]     0     0     0     0     0     0     0     0     0   220     0     0
#>  [73]   264     0     0     0    11     0     0   193     0    25     0    23
#>  [85]     0     0     0     0     0    24     0    62     0   188     0     6
#>  [97]     5     0     0     0
# }
```
