# EVINB summary function

EVINB summary function

## Usage

``` r
# S3 method for class 'evinb'
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

  an EVINB object

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

An EVINB summary object

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
model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
summary(model)
#> EVINB model summary
#> ===================
#> 
#> Count component (negative binomial)
#>             Estimate Std. Error approx t Pr(boot)
#> (Intercept)   2.2479     0.5733    3.921     <0.3
#> x1            1.3922     0.2197    6.337     <0.3
#> x2            0.1915     0.2479    0.773        1
#> x3            0.3861     0.2857    1.351     <0.3
#> 
#> Extreme-value inflation component
#>             Estimate Std. Error approx t Pr(boot)
#> (Intercept)  -2.6905     0.2304  -11.675     <0.3
#> x1            0.9018     0.2826    3.192     <0.3
#> x2           -0.5799     0.5431   -1.068     <0.3
#> x3            0.3535     0.1207    2.929     <0.3
#> 
#> Pareto (extreme value) component
#>             Estimate Std. Error approx t Pr(boot)
#> (Intercept)  2.80276    0.09925   28.240     <0.3
#> x1          -2.44928    0.16461  -14.879     <0.3
#> x2           1.57041    0.09964   15.762     <0.3
#> x3           1.67644    0.47380    3.538     <0.3
#> 
#> ----------------------------------------
#> alpha_NB: 10.05   C_EV: 184
#> Observations at or above C_EV: 10
#> Mean state proportions:  count = 0.920   evi = 0.080
#> Observations: 100   Parameters: 14   df: 86
#> logLik: -261.1   AIC: 550.2   BIC: 586.7   Converged: TRUE
#> Bootstraps: failed = 0, degenerate = 2
# }
```
