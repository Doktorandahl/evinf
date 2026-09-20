# Compare the extreme-value model against its plainer counterparts

For each compared model (`nb`, `zinb`, and any winsorized / razorized
variants) and each metric, computes the paired bootstrap difference
`evinf - compared` (so a negative median favours the extreme-value
model, matching Table B3 of the ISQ appendix), the proportion of
bootstraps in which the evinf model is better, and the number of
bootstrap pairs where both fits succeeded.

## Usage

``` r
compare_fit(
  comp,
  metrics = c("aic", "bic", "rmse", "rmsle"),
  exclude_degenerate = TRUE,
  ...
)

# S3 method for class 'evzinbcomp'
plot(
  x,
  metrics = c("aic", "bic", "rmse", "rmsle"),
  exclude_degenerate = TRUE,
  ...
)
```

## Arguments

- comp:

  An `evzinbcomp` object from [`compare_models()`](compare_models.md).

- metrics:

  Metrics to show.

- exclude_degenerate:

  Drop bootstrap replicates of the evinf model flagged degenerate
  (default `TRUE`); see the `alpha_floor` argument of
  [`evinf_control`](evinf_control.md). Excluded replicates are treated
  as missing rather than dropped, so the pairing with the compared
  models' replicates is preserved.

- ...:

  Unused.

- x:

  An `evzinbcomp` object.

## Value

A tibble of class `evinf_compare_fit` with columns `model`, `metric`,
`median_difference`, `prop_evinf_better`, `n_pairs`. `median_difference`
/ `prop_evinf_better` are `NA` for the `aic` / `bic` rows of a `*_razor`
or `*_winsor` slot (see Details); `n_pairs` is left as computed.

## Details

The `aic` / `bic` rows are `NA` for a `*_razor` or `*_winsor` slot
(round8 0.5): a razorised fit is estimated on fewer observations than
the evinf model, and a winsorised fit is estimated on a different
outcome, so their AIC/BIC are not on the same scale as the evinf model's
and a difference between them is not meaningful. The RMSE / RMSLE rows
for those slots remain comparable, because out-of-bag error is always
computed against the raw (un-winsorised, un-razorised) outcome.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 10)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 2 of 10 bootstrap replicates; consider widening c.lim.
compare_fit(compare_models(model))
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
#> Warning: step size truncated due to divergence
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
#> Paired bootstrap fit comparison (evinf - compared)
#>   negative median favours the extreme-value model
#> 
#>  model metric median_difference prop_evinf_better n_pairs
#>     nb    aic          -76.0200             1.000       3
#>     nb    bic          -42.1500             0.667       3
#>     nb   rmse           33.4700             0.000       3
#>     nb  rmsle            0.2116             0.000       3
#>   zinb    aic          -26.5600             1.000       3
#>   zinb    bic           -3.1160             0.667       3
#>   zinb   rmse           62.5200             0.000       3
#>   zinb  rmsle            0.2452             0.000       3
# }
```
