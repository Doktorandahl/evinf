# EVZINB and EVINB glance functions

EVZINB and EVINB glance functions

## Usage

``` r
# S3 method for class 'evzinb'
glance(x, ...)
```

## Arguments

- x:

  An EVZINB or EVINB object

- ...:

  Further arguments to be passed to glance()

## Value

A one-row tibble of goodness-of-fit statistics: number of observations
and parameters, alpha_NB, C_EV, AIC, BIC, log-likelihood, whether the EM
algorithm converged, the number of EM steps, the number of observations
at or above C_EV, and the number of usable, failed and degenerate
bootstraps (`NA` when the model was fitted without bootstrapping; see
[`evinf_control`](evinf_control.md) for "degenerate").

## See also

[`glance`](https://generics.r-lib.org/reference/glance.html)

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
glance(model)
#> # A tibble: 1 × 13
#>    nobs  npar alpha parameter   aic   bic logLik converged n_above_c n_em_steps
#>   <int> <int> <dbl>     <dbl> <dbl> <dbl>  <dbl> <lgl>         <int>      <int>
#> 1   100    18  1.43       184  544.  591.  -254. TRUE             10         16
#> # ℹ 3 more variables: n_bootstraps <int>, n_failed_bootstraps <int>,
#> #   n_degenerate_bootstraps <int>
# }
```
