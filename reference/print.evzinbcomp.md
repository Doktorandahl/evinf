# Print method for compare_models() output

Shows the compared models, a compact
[`glance()`](https://generics.r-lib.org/reference/glance.html) table
(one row per slot, key fit statistics only), and the paired bootstrap
comparison from [`compare_fit()`](compare_fit.md) (round10 I.4, audit
§5.9) – previously only the compared-model names and bootstrap count.

## Usage

``` r
# S3 method for class 'evzinbcomp'
print(x, metrics = c("aic", "bic", "rmse", "rmsle"), ...)
```

## Arguments

- x:

  An `evzinbcomp` object returned by
  [`compare_models`](compare_models.md).

- metrics:

  Metrics for the [`compare_fit()`](compare_fit.md) table; see there.

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
print(compare_models(model))
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: alternation limit reached
#> 
#> Model comparison of evzinb
#>  Compared models: nb, zinb
#>  Number of compared models: 2
#>  Number of bootstraps:5
#> 
#> Fit summary
#>  model nobs npar logLik   aic   bic
#>  model  100   18 -254.0 544.1 591.0
#>     nb  100    4 -273.2 556.4 569.4
#>   zinb  100    9 -266.1 550.3 573.7
#> 
#> Paired bootstrap fit comparison (evinf - compared)
#>   negative median favours the extreme-value model
#> 
#>  model metric median_difference prop_evinf_better n_pairs
#>     nb    aic        -3807.0000               1.0       2
#>     nb    bic        -3773.0000               0.5       2
#>     nb   rmse         1127.0000               0.0       2
#>     nb  rmsle            0.6354               0.0       2
#>   zinb    aic          -16.3500               1.0       2
#>   zinb    bic            7.0950               0.5       2
#>   zinb   rmse         1135.0000               0.0       2
#>   zinb  rmsle            0.6976               0.0       2
# }
```
