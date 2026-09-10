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

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
summary(model)
#> EVZINB model summary
#> ====================
#> 
#> Count component (negative binomial)
#>             Estimate Std. Error approx t Pr(boot)    
#> (Intercept)   3.1305     0.4151    7.541   <2e-16 ***
#> x1            0.8432     0.4614    1.828   <2e-16 ***
#> x2            0.5987     0.3929    1.524   <2e-16 ***
#> x3            0.2411     0.1893    1.274      0.5    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Zero-inflation component
#>             Estimate Std. Error approx t Pr(boot)    
#> (Intercept)   0.5109     0.4517    1.131      0.5    
#> x1           -0.6738     0.2344   -2.874   <2e-16 ***
#> x2            0.6071     0.3814    1.592   <2e-16 ***
#> x3           -0.4197     0.2767   -1.517   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Extreme-value inflation component
#>             Estimate Std. Error approx t Pr(boot)    
#> (Intercept)  -1.7286     0.6390  -2.7051   <2e-16 ***
#> x1            0.7356     0.6055   1.2149   <2e-16 ***
#> x2           -0.4289     0.5719  -0.7499   <2e-16 ***
#> x3            0.1662     0.1819   0.9136   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Pareto (extreme value) component
#>             Estimate Std. Error approx t Pr(boot)    
#> (Intercept)   2.7136     0.1481   18.327   <2e-16 ***
#> x1           -2.2778     0.3571   -6.378   <2e-16 ***
#> x2            1.4488     0.2457    5.898   <2e-16 ***
#> x3            1.5035     0.6040    2.489   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> ----------------------------------------
#> alpha_NB: 1.431   C_EV: 184
#> Observations at or above C_EV: 10
#> Mean state proportions:  zero = 0.575   count = 0.340   evi = 0.084
#> Observations: 100   Parameters: 18   df: 82
#> logLik: -254   AIC: 544.1   BIC: 591   Converged: TRUE
#> Bootstraps: failed = 1, degenerate = 0
# }
```
