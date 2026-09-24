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
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
revzinb_fit(model)
#>   [1]   5   0   0   0   0   0   0   0   0  52   0   0 155   0   0   6   0   0
#>  [19]   8   0   0   0   0   0   3   8   0   0   0   0   0   0   0   0   9   0
#>  [37]   0  32   0 365   0   0   0   9   0   0 190  18 206 184   0 149   0   0
#>  [55]   8   0   0   0  86 427   0   0   0 276 699   0   0   2 192  85   0   0
#>  [73]   1 219   0  43  27   0   0   4 190 186   0   0   0   0   0 374   0   0
#>  [91]   0 235  44   6   0  60  37   8   0   0
# }
```
