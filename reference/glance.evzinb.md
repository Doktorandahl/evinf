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
and parameters, family (round9 E.0, e.g. "nbinom/mixture" or
"poisson/mixture"), alpha_NB (NA when the family has none, round9 E.1),
C_EV, AIC, BIC, log-likelihood, whether the EM algorithm converged
(`converged`) and whether the C_EV profile settled within `max.c.iter`
(`c_converged`; `NA` for a model fitted before this field existed), the
number of EM steps, the number of observations at or above C_EV, the
smallest fitted Pareto shape (`min_alpha_pl`; see `alpha_pl_floor` in
[`evinf_control`](evinf_control.md)), `sum_weights` (round9 D.2: the sum
of `weights =`, equal to `nobs` for an unweighted fit – this, not
`nobs`, is what `aic`/`bic` and the model's degrees of freedom are
computed from), and, for a bootstrapped model, the bootstrap replicate
counts (`NA` otherwise): `n_bootstraps` (usable), `n_failed_bootstraps`,
`n_degenerate_bootstraps` — these three partition the number of
replicates requested (see [`evinf_control`](evinf_control.md) for
"degenerate") — and `n_c_on_boundary`, the number of replicates whose
C_EV landed on a candidate-grid endpoint (informational, not counted as
degenerate); `oob_fraction_mean`/`oob_fraction_min`/ `oob_fraction_max`
(round10 0.9), the mean/min/max of each usable replicate's out-of-bag
row fraction (`NA` without bootstraps) – most informative for the block
schemes, where overlapping blocks change how much of the data a
replicate leaves out; `n_starts`/ `n_starts_at_best` (round10 G.1, `NA`
for `n_starts <= 1`), how many perturbed starts were tried and how many
of them reached the winning log-likelihood (within 1e-4);
`median_boot_em_steps`/ `n_boot_c_capped` (round10 G.3, `NA` without
bootstraps), the median EM step count across usable replicates and how
many had their C_EV profile hit `max.c.iter` without settling in either
phase.

## See also

[`glance`](https://generics.r-lib.org/reference/glance.html)

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
glance(model)
#> # A tibble: 1 × 25
#>    nobs sum_weights  npar family    alpha parameter   aic   bic logLik converged
#>   <int>       <dbl> <int> <chr>     <dbl>     <dbl> <dbl> <dbl>  <dbl> <lgl>    
#> 1   100         100    18 nbinom/m…  1.43       184  544.  591.  -254. TRUE     
#> # ℹ 15 more variables: c_converged <lgl>, n_above_c <int>, n_em_steps <int>,
#> #   min_alpha_pl <dbl>, n_bootstraps <int>, n_failed_bootstraps <int>,
#> #   n_degenerate_bootstraps <int>, n_c_on_boundary <int>,
#> #   oob_fraction_mean <dbl>, oob_fraction_min <dbl>, oob_fraction_max <dbl>,
#> #   n_starts <int>, n_starts_at_best <int>, median_boot_em_steps <dbl>,
#> #   n_boot_c_capped <int>
# }
```
