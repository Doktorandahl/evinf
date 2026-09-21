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
revinb_fit(model)
#>   [1]       1      14     184       0      46       0      39     755       0
#>  [10]       0       0       0      90       0       0       0      28       0
#>  [19] 1008840       6       0       3       0       7       0       1     184
#>  [28]       0       0       0     109       1       0       0       0       6
#>  [37]       0       0     200       0      14       0       0       0      10
#>  [46]     196       0       0     277       0       0     184      10       0
#>  [55]       0       0       0       0     298       0       0     184       0
#>  [64]      22       0       0       1       0       0       3       0       1
#>  [73]      11      29       0       0       0       0       0       0      17
#>  [82]       0       0       1       0       3       0       1       9       0
#>  [91]       0     187      37      10      52       0     526       0       0
#> [100]       0
# }
```
