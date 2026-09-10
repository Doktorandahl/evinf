# logLik / nobs / model.frame / formula / terms for evzinb / evinb models

logLik / nobs / model.frame / formula / terms for evzinb / evinb models

## Usage

``` r
# S3 method for class 'evzinb'
logLik(object, ...)

# S3 method for class 'evinb'
logLik(object, ...)

# S3 method for class 'evzinb'
nobs(object, ...)

# S3 method for class 'evinb'
nobs(object, ...)

# S3 method for class 'evzinb'
model.frame(formula, ...)

# S3 method for class 'evinb'
model.frame(formula, ...)

# S3 method for class 'evzinb'
formula(x, component = "count", ...)

# S3 method for class 'evinb'
formula(x, component = "count", ...)

# S3 method for class 'evzinb'
terms(x, component = "count", ...)

# S3 method for class 'evinb'
terms(x, component = "count", ...)
```

## Arguments

- object:

  A fitted model.

- ...:

  Unused.

- formula:

  A fitted model (for the `formula` method).

- x:

  A fitted model (for the `terms` method).

- component:

  Which component (`"count"` by default) for
  [`formula()`](https://rdrr.io/r/stats/formula.html) /
  [`terms()`](https://rdrr.io/r/stats/terms.html).

## Value

As for the corresponding stats generic.
[`logLik()`](https://rdrr.io/r/stats/logLik.html) carries
`df = length(object$par.all)` and `nobs`, so
[`AIC()`](https://rdrr.io/r/stats/AIC.html) /
[`BIC()`](https://rdrr.io/r/stats/AIC.html) work and match `object$AIC`
/ `object$BIC`.
