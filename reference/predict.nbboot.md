# Prediction for nbboot

Prediction for nbboot

## Usage

``` r
# S3 method for class 'nbboot'
predict(
  object,
  newdata = NULL,
  type = c("predicted", "all", "quantile"),
  pred = c("original", "bootstrap_median", "bootstrap_mean"),
  quantile = NULL,
  confint = FALSE,
  conf_level = 0.9,
  ...
)
```

## Arguments

- object:

  a fitted nbboot object

- newdata:

  Data to make predictions on

- type:

  What prediction should be computed?

- pred:

  Prediction type, 'original', 'bootstrap_median', or 'bootstrap_mean'

- quantile:

  Quantile for quantile prediction

- confint:

  Should confidence intervals be created?

- conf_level:

  Confidence level when predicting with CIs

- ...:

  Not used

## Value

Predictions from nbboot
