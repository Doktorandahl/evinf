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
algorithm converged (`converged`) and whether the C_EV profile settled
within `max.c.iter` (`c_converged`; `NA` for a model fitted before this
field existed), the number of EM steps, the number of observations at or
above C_EV, and, for a bootstrapped model, the bootstrap replicate
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
#> Warning: C_EV equalled the lower endpoint (173) in 2 of 5 bootstrap replicates; consider widening c.lim.
glance(model)
#> # A tibble: 1 × 15
#>    nobs  npar alpha parameter   aic   bic logLik converged c_converged n_above_c
#>   <int> <int> <dbl>     <dbl> <dbl> <dbl>  <dbl> <lgl>     <lgl>           <int>
#> 1   100    14  10.0       184  550.  587.  -261. TRUE      TRUE               10
#> # ℹ 5 more variables: n_em_steps <int>, n_bootstraps <int>,
#> #   n_failed_bootstraps <int>, n_degenerate_bootstraps <int>,
#> #   n_c_on_boundary <int>
# }
```
