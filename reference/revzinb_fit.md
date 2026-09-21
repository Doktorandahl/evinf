# Random draws from a fitted evzinb model

**Superseded:** kept for backwards compatibility, with no plans for
removal, but new code should use [`simulate()`](evinf-s3-predict.md),
which returns a tidy data frame (and a reproducible `"seed"` attribute)
instead of a bare vector/list.

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

## See also

[`simulate.evzinb`](evinf-s3-predict.md)

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3, data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
revzinb_fit(model)
#>   [1]    0    0    0    1    0    0   24  201    5    0    0    0  567    5   29
#>  [16]    0    5   33   61   25    0    0    0    0    0    8    7    6   46    0
#>  [31]   53    7    0   39    0    0    0    0    0    6    0    0    0    7    0
#>  [46]  187    0   35  211    0   39    0   61    0    3    0    0 1108   88   12
#>  [61]    0    0    5   24    8    0    0    0    0    5    0    0   50   68    0
#>  [76]    0   92    0    0    0    0    7    0    0    0    0    0    0    0    2
#>  [91]    0    0    0  106   95   43  419    0    3    0
# }
```
