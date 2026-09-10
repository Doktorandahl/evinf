# Fitted values, residuals and simulations from an evzinb / evinb model

Fitted values, residuals and simulations from an evzinb / evinb model

## Usage

``` r
# S3 method for class 'evzinb'
fitted(object, type = c("harmonic", "explog", "counts", "pareto_alpha"), ...)

# S3 method for class 'evinb'
fitted(object, type = c("harmonic", "explog", "counts", "pareto_alpha"), ...)

# S3 method for class 'evzinb'
residuals(object, type = c("response", "quantile"), seed = NULL, ...)

# S3 method for class 'evinb'
residuals(object, type = c("response", "quantile"), seed = NULL, ...)

# S3 method for class 'evzinb'
simulate(object, nsim = 1, seed = NULL, newdata = NULL, ...)

# S3 method for class 'evinb'
simulate(object, nsim = 1, seed = NULL, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted model.

- type:

  For [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), one of
  `"harmonic"`, `"explog"`, `"counts"`, `"pareto_alpha"`. For
  [`residuals()`](https://rdrr.io/r/stats/residuals.html), `"response"`
  (\\y - \\ harmonic prediction) or `"quantile"` (randomized quantile
  residuals from the mixture CDF).

- ...:

  Unused.

- seed:

  Optional RNG seed for the randomized quantile residuals /
  [`simulate()`](https://rdrr.io/r/stats/simulate.html).

- nsim:

  Number of simulated response vectors.

- newdata:

  Optional data to simulate for.

## Value

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) /
[`residuals()`](https://rdrr.io/r/stats/residuals.html) return a numeric
vector; [`simulate()`](https://rdrr.io/r/stats/simulate.html) a data
frame with columns `sim_1`, `sim_2`, ...
