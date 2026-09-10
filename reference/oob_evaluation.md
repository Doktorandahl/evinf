# Out of bag predictive performance of EVZINB and EVINB models

Out of bag predictive performance of EVZINB and EVINB models

## Usage

``` r
oob_evaluation(
  object,
  predict_type = c("harmonic", "explog"),
  metric = c("rmsle", "rmse", "mse", "mae")
)
```

## Arguments

- object:

  A fitted evzinb / evinb model with bootstraps, or an `evzinbcomp`
  object from [`compare_models()`](compare_models.md).

- predict_type:

  What type of prediction should be made? Harmonic mean, or
  exp(log(prediction))?

- metric:

  What metric should be used for the out of bag evaluation? Default
  options include rmsle, rmse, mse, and mae. Can also take a user
  supplied function of the form function(y_pred,y_true) which returns a
  single value

## Value

For a single model, a vector of length `n_bootstraps`. For an
`evzinbcomp` object, a tibble with one column per compared model
(`evinf`, `nb`, `zinb`, ...) and one row per bootstrap.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
oob_evaluation(model)
#> Error in bootstrap$boot_id : $ operator is invalid for atomic vectors
#> [1] 3.053584 8.654060       NA 3.339729 4.003586
# }
```
