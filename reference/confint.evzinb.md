# Confidence intervals for an evzinb / evinb model

Confidence intervals for an evzinb / evinb model

## Usage

``` r
# S3 method for class 'evzinb'
confint(
  object,
  parm,
  level = 0.95,
  type = c("percentile", "approx"),
  exclude_degenerate = TRUE,
  ...
)

# S3 method for class 'evinb'
confint(
  object,
  parm,
  level = 0.95,
  type = c("percentile", "approx"),
  exclude_degenerate = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted model with bootstraps.

- parm:

  Which parameters (names as in `coef(object)`); default all.

- level:

  Confidence level.

- type:

  `"percentile"` (bootstrap percentile intervals, the default) or
  `"approx"` (`estimate +/- qnorm() * sqrt(diag(vcov))`).

- exclude_degenerate:

  Drop bootstrap replicates flagged degenerate (default `TRUE`); see
  [`evinf_control`](evinf_control.md).

- ...:

  Unused.

## Value

A two-column matrix.

## Examples

``` r
# \donttest{
data(genevzinb2)
m <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 25)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the upper candidate endpoint (263) in 1 and the lower endpoint (173) in 4 of 25 bootstrap replicates; consider widening c.lim.
# the Pareto shape (alpha_nb) and the threshold (c_ev) are in coef()/confint()
confint(m, parm = c("count_x1", "alpha_nb", "c_ev"))
#>                2.5 %     97.5 %
#> count_x1   0.3056582   1.832540
#> alpha_nb   1.0615163   4.466584
#> c_ev     184.0000000 190.775000
# }
```
