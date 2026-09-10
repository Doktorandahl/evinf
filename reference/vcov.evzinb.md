# Bootstrap covariance matrix of an evzinb / evinb model

Bootstrap covariance matrix of an evzinb / evinb model

## Usage

``` r
# S3 method for class 'evzinb'
vcov(object, exclude_degenerate = TRUE, ...)

# S3 method for class 'evinb'
vcov(object, exclude_degenerate = TRUE, ...)
```

## Arguments

- object:

  A fitted model with bootstraps.

- exclude_degenerate:

  Drop bootstrap replicates flagged degenerate (default `TRUE`); see
  [`evinf_control`](evinf_control.md).

- ...:

  Unused.

## Value

The covariance matrix of the bootstrap coefficient draws, with the same
names and order as `coef(object)` (includes `alpha_nb` and `c_ev`).

## Details

`c_ev` is estimated on the grid of unique observed values of the
response, not by a smooth optimiser, so it sits on a discrete scale. Its
bootstrap variance (and a percentile
[`confint()`](https://rdrr.io/r/stats/confint.html) on it) is
meaningful, but delta-method standard errors that perturb `c_ev`
continuously – as
[`marginaleffects::avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html)
etc. do – should be read with that in mind. For uncertainty on covariate
*effects*, prefer the bootstrap route,
[`marginal_effects()`](marginal_effects.md).
