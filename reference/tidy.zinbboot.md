# Tidy function for zinbboot

Tidy function for zinbboot

## Usage

``` r
# S3 method for class 'zinbboot'
tidy(
  x,
  component = c("all", "count", "zero"),
  coef_type = c("original", "bootstrap_mean", "bootstrap_median"),
  standard_error = TRUE,
  p_value = c("bootstrapped", "approx", "none"),
  confint = c("none", "bootstrapped", "approx"),
  conf_level = 0.95,
  approx_t_value = TRUE,
  symmetric_bootstrap_p = TRUE,
  ...
)
```

## Arguments

- x:

  A fitted bootstrapped zero-inflated model

- component:

  Which component should be shown?

- coef_type:

  What type of coefficient should be reported, original, bootstrapped
  mean, or bootstrapped median

- standard_error:

  Should bootstrapped standard errors be reported?

- p_value:

  What type of p-value should be reported? Bootstrapped p_values,
  approximate p-values, or none?

- confint:

  What type of confidence intervals should be reported? Bootstrapped
  p_values, approximate p-values, or none?

- conf_level:

  Confidence level for confidence intervals

- approx_t_value:

  Should approximate t_values be reported

- symmetric_bootstrap_p:

  Should bootstrap p-values be computed as symmetric (leaving alpha/2
  percent in each tail)? FALSE gives non-symmetric, but narrower,
  intervals. TRUE corresponds most closely to conventional p-values.

- ...:

  Other arguments to be passed to tidy

## Value

A tidy function for a bootstrapped zinb model

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
zinb_comp <- compare_models(model)
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
#> Warning: glm.fit: algorithm did not converge
tidy(zinb_comp$zinb)
#> # A tibble: 8 × 6
#>   y.level term        estimate std.error statistic p.value
#>   <chr>   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
#> 1 zero    (Intercept)   0.288      2.84     0.101      1  
#> 2 zero    x1           -0.835      1.09    -0.767      0.2
#> 3 zero    x2            0.603      0.736    0.819      0.4
#> 4 zero    x3           -0.493      1.94    -0.254      0.2
#> 5 count   (Intercept)   4.07       0.352   11.6        0.2
#> 6 count   x1            0.767      0.491    1.56       0.2
#> 7 count   x2           -0.151      0.210   -0.718      0.4
#> 8 count   x3            0.0373     0.529    0.0705     0.8
# }
```
