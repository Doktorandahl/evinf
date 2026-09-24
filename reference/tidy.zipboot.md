# Tidy function for zipboot

Tidy function for zipboot

## Usage

``` r
# S3 method for class 'zipboot'
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

  A fitted bootstrapped zero-inflated Poisson model

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

A tidy function for a bootstrapped ZIP model

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, family = "poisson", n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 4 of 5 bootstrap replicates; consider widening c.lim.
zip_comp <- compare_models(model)
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
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: alternation limit reached
tidy(zip_comp$zip)
#> # A tibble: 8 × 6
#>   y.level term        estimate std.error statistic p.value
#>   <chr>   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
#> 1 zero    (Intercept)   0.508      0.402     1.26      0.2
#> 2 zero    x1           -0.855      0.217    -3.94      0.2
#> 3 zero    x2            0.575      0.103     5.59      0.2
#> 4 zero    x3           -0.455      0.251    -1.81      0.2
#> 5 count   (Intercept)   4.24       0.334    12.7       0.2
#> 6 count   x1            0.599      0.422     1.42      0.2
#> 7 count   x2           -0.0685     0.362    -0.189     0.4
#> 8 count   x3            0.0700     0.405     0.173     0.4
# }
```
