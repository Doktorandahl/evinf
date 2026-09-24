# zipboot and poissonboot glance functions

zipboot and poissonboot glance functions

## Usage

``` r
# S3 method for class 'zipboot'
glance(x, ...)

# S3 method for class 'poissonboot'
glance(x, ...)
```

## Arguments

- x:

  A poissonboot or zipboot object

- ...:

  Further arguments to be passed to glance()

## Value

A one-row tibble of goodness-of-fit statistics, including whether the
full-sample model converged and the number of (failed) bootstraps.
`alpha` is `NA`: neither baseline has a dispersion parameter (round9
E.3).

## See also

[`glance`](https://generics.r-lib.org/reference/glance.html)

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, family = "poisson", n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
zip_comp <- compare_models(model)
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
glance(zip_comp$zip)
#> # A tibble: 1 × 9
#>    nobs  npar alpha   aic   bic logLik converged n_bootstraps
#>   <int> <int> <dbl> <dbl> <dbl>  <dbl> <lgl>            <int>
#> 1   100     8    NA 4648. 4669. -2316. TRUE                 5
#> # ℹ 1 more variable: n_failed_bootstraps <int>
# }
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, family = "poisson", n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 3 of 5 bootstrap replicates; consider widening c.lim.
zip_comp <- compare_models(model)
#> Warning: glm.fit: algorithm did not converge
glance(zip_comp$poisson)
#> # A tibble: 1 × 9
#>    nobs  npar alpha   aic   bic logLik converged n_bootstraps
#>   <int> <int> <dbl> <dbl> <dbl>  <dbl> <lgl>            <int>
#> 1   100     4    NA 9576. 9586. -4784. TRUE                 5
#> # ℹ 1 more variable: n_failed_bootstraps <int>
# }
```
