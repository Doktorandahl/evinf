# EVZINB summary function

EVZINB summary function

## Usage

``` r
# S3 method for class 'evzinb'
summary(
  object,
  coef = c("original", "bootstrapped_mean", "bootstrapped_median"),
  standard_error = TRUE,
  p_value = c("bootstrapped", "approx", "both", "none"),
  bootstrapped_props = c("none", "mean", "median"),
  approx_t_value = TRUE,
  symmetric_bootstrap_p = TRUE,
  exclude_degenerate = TRUE,
  ...
)
```

## Arguments

- object:

  an EVZINB object

- coef:

  Type of coefficients. Original are the coefficient estimates from the
  non-bootstrapped version of the model. 'bootstrapped_mean' are the
  mean coefficients across bootstraps, and 'bootstrapped_median' are the
  median coefficients across bootstraps

- standard_error:

  Should standard errors be computed?

- p_value:

  What type of p_values should be computed? 'bootstrapped' are
  bootstrapped p_values through confidence interval inversion. 'approx'
  are p-values based on the t-value produced by dividing the coefficient
  with the standard error. 'both' returns both.

- bootstrapped_props:

  Type of bootstrapped proportions of component proportions to be
  returned

- approx_t_value:

  Should approximate t-values be returned

- symmetric_bootstrap_p:

  Should bootstrap p-values be computed as symmetric (leaving alpha/2
  percent in each tail)? FALSE gives non-symmetric, but narrower,
  intervals. TRUE corresponds most closely to conventional p-values.

- exclude_degenerate:

  Drop bootstrap replicates flagged degenerate (default TRUE); see the
  alpha_floor argument of evinf_control().

- ...:

  Additional arguments passed to the summary function

## Value

An EVZINB summary object

## Details

When `object` was fitted with `bootstrap = FALSE`, the bootstrap-based
quantities (standard errors, p-values, approximate t-values and
bootstrapped proportions) are unavailable; the summary then reports the
point estimates only, `n_failed_bootstraps` is `NA`, and a message is
emitted. Requesting bootstrapped coefficients for such a model is an
error.

A bootstrapped p-value is never reported below `1 / B`, where `B` is the
number of usable bootstrap replicates: with no draw crossing the
estimate, the true p-value could be anywhere in `[0, 1/B)`, so it is
floored at `1/B` rather than reported as exactly `0`.
[`print.summary.evzinb()`](print.summary.evzinb.md) /
[`print.summary.evinb()`](print.summary.evzinb.md) show this as
`"< 1/B"` (e.g. `"<0.01"` for 100 usable bootstraps) rather than the
default, misleadingly precise `"<2e-16"`.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
summary(model)
#> EVZINB model summary
#> ====================
#> 
#> Count component (negative binomial)
#>             Estimate Std. Error approx t Pr(boot)
#> (Intercept)   3.1305     0.3708    8.442     <0.2
#> x1            0.8432     0.4595    1.835     <0.2
#> x2            0.5987     0.3678    1.628     <0.2
#> x3            0.2411     0.4552    0.530      0.4
#> 
#> Zero-inflation component
#>             Estimate Std. Error approx t Pr(boot)
#> (Intercept)   0.5109     0.5085    1.005      0.4
#> x1           -0.6738     0.4244   -1.587     <0.2
#> x2            0.6071     0.3350    1.812     <0.2
#> x3           -0.4197     0.3242   -1.295      0.4
#> 
#> Extreme-value inflation component
#>             Estimate Std. Error approx t Pr(boot)
#> (Intercept)  -1.7286     3.7981   -0.455     <0.2
#> x1            0.7356     3.8051    0.193     <0.2
#> x2           -0.4289     2.2152   -0.194     <0.2
#> x3            0.1662     0.6213    0.268     <0.2
#> 
#> Pareto (extreme value) component
#>             Estimate Std. Error approx t Pr(boot)
#> (Intercept)    2.714      3.584    0.757     <0.2
#> x1            -2.278      3.831   -0.595     <0.2
#> x2             1.449      1.771    0.818     <0.2
#> x3             1.503      2.696    0.558     <0.2
#> 
#> ----------------------------------------
#> alpha_NB: 1.431   C_EV: 184
#> Observations at or above C_EV: 10
#> Mean state proportions:  zero = 0.575   count = 0.341   evi = 0.084
#> Observations: 100   Parameters: 18   df: 82
#> logLik: -254   AIC: 544.1   BIC: 591   Converged: TRUE
#> Bootstraps: failed = 0, degenerate = 0
# }
```
