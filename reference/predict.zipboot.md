# Prediction for zipboot

Prediction for zipboot

## Usage

``` r
# S3 method for class 'zipboot'
predict(
  object,
  newdata = NULL,
  type = c("predicted", "counts", "zi", "count_state", "states", "all", "quantile"),
  pred = c("original", "bootstrap_median", "bootstrap_mean"),
  quantile = NULL,
  confint = FALSE,
  conf_level = 0.9,
  ...
)
```

## Arguments

- object:

  a fitted zipboot object

- newdata:

  Data to make predictions on

- type:

  What prediction should be computed? One of `"predicted"`, `"counts"`,
  `"zi"`, `"count_state"`, `"states"`, `"all"` or `"quantile"`. (A
  zero-inflated Poisson has no extreme-value state, so `"evinf"` is not
  accepted.)

- pred:

  Prediction type, 'original', 'bootstra_median', or 'bootstrap_mean'

- quantile:

  Quantile for quantile prediction

- confint:

  Should confidence intervals be created?

- conf_level:

  Confidence level when predicting with CIs

- ...:

  Not used

## Value

Predictions from zipboot
