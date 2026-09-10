# Compare the extreme-value model against its plainer counterparts

For each compared model (`nb`, `zinb`, and any winsorized / razorized
variants) and each metric, computes the paired bootstrap difference
`compared - evinf` (so a negative median favours the extreme-value
model), the proportion of bootstraps in which the evinf model is better,
and the number of bootstrap pairs where both fits succeeded.

## Usage

``` r
compare_fit(comp, metrics = c("aic", "bic", "rmse", "rmsle"), ...)

# S3 method for class 'evzinbcomp'
plot(x, metrics = c("aic", "bic", "rmse", "rmsle"), ...)
```

## Arguments

- comp:

  An `evzinbcomp` object from [`compare_models()`](compare_models.md).

- metrics:

  Metrics to show.

- ...:

  Unused.

- x:

  An `evzinbcomp` object.

## Value

A tibble of class `evinf_compare_fit` with columns `model`, `metric`,
`median_difference`, `prop_evinf_better`, `n_pairs`.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 10)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
#> Error in eval(expr, p) : inv(): matrix is singular
#> Warning: C_EV reached the boundary of the candidate range in 1 of 10 bootstrap replicates; consider widening c.lim.
compare_fit(compare_models(model))
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Error in -bootstrap$boot_id : invalid argument to unary operator
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: alternation limit reached
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
#> Error in bootstrap$boot_id : $ operator is invalid for atomic vectors
#> Error in bootstrap$boot_id : $ operator is invalid for atomic vectors
#> Error in bootstrap$boot_id : $ operator is invalid for atomic vectors
#> Error in bootstrap$boot_id : $ operator is invalid for atomic vectors
#> Paired bootstrap fit comparison (compared - evinf)
#>   negative median favours the extreme-value model
#> 
#>  model metric median_difference prop_evinf_better n_pairs
#>     nb    aic           45.0400             1.000       8
#>     nb    bic           11.1700             0.750       8
#>     nb   rmse          -77.6200             0.250       8
#>     nb  rmsle           -0.3111             0.125       8
#>   zinb    aic           25.8100             0.875       8
#>   zinb    bic            2.3680             0.625       8
#>   zinb   rmse          -91.8300             0.000       8
#>   zinb  rmsle           -0.4376             0.000       8
# }
```
