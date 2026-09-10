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

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
summary(model)
#> EVINB model summary
#> ===================
#> 
#> Count component (negative binomial)
#>             Estimate Std. Error approx t Pr(boot)    
#> (Intercept)   2.2479     0.1112   20.219   <2e-16 ***
#> x1            1.3922     1.1606    1.200   <2e-16 ***
#> x2            0.1915     0.9433    0.203    0.667    
#> x3            0.3861     0.3836    1.007   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Extreme-value inflation component
#>             Estimate Std. Error approx t Pr(boot)    
#> (Intercept)  -2.6905     0.3615   -7.443   <2e-16 ***
#> x1            0.9018     0.5262    1.714   <2e-16 ***
#> x2           -0.5799     0.4139   -1.401   <2e-16 ***
#> x3            0.3535     0.1724    2.050   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Pareto (extreme value) component
#>             Estimate Std. Error approx t Pr(boot)    
#> (Intercept)   2.8028     0.2285   12.269   <2e-16 ***
#> x1           -2.4493     0.8955   -2.735   <2e-16 ***
#> x2            1.5704     0.6219    2.525   <2e-16 ***
#> x3            1.6764     0.5166    3.245   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ----------------------------------------
#> alpha_NB: 10.05   C_EV: 184
#> Observations at or above C_EV: 10
#> Mean state proportions:  count = 0.920   evi = 0.080
#> Observations: 100   Parameters: 14   df: 86
#> logLik: -261.1   AIC: 550.2   BIC: 586.7   Converged: TRUE
#> Bootstraps: failed = 1, degenerate = 1
# }
```
