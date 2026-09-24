# Fitted values, residuals and simulations from an evzinb / evinb model

Fitted values, residuals and simulations from an evzinb / evinb model

## Usage

``` r
# S3 method for class 'evzinb'
fitted(object, type = c("harmonic", "explog", "counts", "pareto_alpha"), ...)

# S3 method for class 'evinb'
fitted(object, type = c("harmonic", "explog", "counts", "pareto_alpha"), ...)

# S3 method for class 'evzinb'
residuals(
  object,
  type = c("response", "quantile"),
  seed = NULL,
  newdata = NULL,
  ...
)

# S3 method for class 'evinb'
residuals(
  object,
  type = c("response", "quantile"),
  seed = NULL,
  newdata = NULL,
  ...
)

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

  For [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  forwarded to [`predict()`](https://rdrr.io/r/stats/predict.html)
  (round11 A4: e.g. `clamp_alpha_pl` for `type = "explog"`); unused
  otherwise.

- seed:

  Optional RNG seed for the randomized quantile residuals /
  [`simulate()`](https://rdrr.io/r/stats/simulate.html). If given, the
  caller's RNG state is restored on exit (the seed only affects this
  call's draws).

- newdata:

  Optional data to simulate for, or (round10 I.2, for
  [`residuals()`](https://rdrr.io/r/stats/residuals.html)) to compute
  residuals for instead of the estimation data; the response column
  (named by the count formula's left-hand side) must be present.

- nsim:

  Number of simulated response vectors.

## Value

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) /
[`residuals()`](https://rdrr.io/r/stats/residuals.html) return a numeric
vector; [`simulate()`](https://rdrr.io/r/stats/simulate.html) a data
frame with columns `sim_1`, `sim_2`, ..., with a `"seed"` attribute
following the [`simulate`](https://rdrr.io/r/stats/simulate.html)
convention (the given `seed`, or, if none was given, the RNG state
before the draws were made).
