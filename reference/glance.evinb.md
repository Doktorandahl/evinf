# EVZINB and EVINB glance functions

EVZINB and EVINB glance functions

## Usage

``` r
# S3 method for class 'evinb'
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
at or above C_EV, and, for a bootstrapped model, the bootstrap replicate
counts (`NA` otherwise): `n_bootstraps` (usable), `n_failed_bootstraps`,
`n_degenerate_bootstraps` — these three partition the number of
replicates requested (see [`evinf_control`](evinf_control.md) for
"degenerate") — and `n_c_on_boundary`, the number of replicates whose
C_EV landed on a candidate-grid endpoint (informational, not counted as
degenerate).

## See also

[`glance`](https://generics.r-lib.org/reference/glance.html)

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
#> Warning: C_EV reached the boundary of the candidate range in 2 of 5 bootstrap replicates; consider widening c.lim.
glance(model)
#> # A tibble: 1 × 14
#>    nobs  npar alpha parameter   aic   bic logLik converged n_above_c n_em_steps
#>   <int> <int> <dbl>     <dbl> <dbl> <dbl>  <dbl> <lgl>         <int>      <int>
#> 1   100    14  10.0       184  550.  587.  -261. TRUE             10         26
#> # ℹ 4 more variables: n_bootstraps <int>, n_failed_bootstraps <int>,
#> #   n_degenerate_bootstraps <int>, n_c_on_boundary <int>
# }
```
