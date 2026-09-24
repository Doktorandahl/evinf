# The predictive distribution of a fitted evzinb / evinb model

`evinf_pmf()` / `evinf_cdf()` are the package's one definition of the
predictive distribution – family-aware across all four (count x zero)
combinations – that
`predict(type = "distribution"/"quantile"/ "exceedance"/"draws")`,
`residuals(type = "quantile")` and
[`classify_states()`](classify_states.md) all build on.

## Usage

``` r
evinf_pmf(
  object,
  newdata = NULL,
  y = NULL,
  support = NULL,
  is_zinb = inherits(object, "evzinb")
)

evinf_cdf(
  object,
  newdata = NULL,
  y = NULL,
  support = NULL,
  is_zinb = inherits(object, "evzinb")
)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` model.

- newdata:

  Optional new data; `NULL` uses the estimation data.

- y:

  Optional numeric vector, one value per row of `newdata` (or the
  estimation data): the pmf/cdf at each row's own `y`. Exactly one of
  `y` or `support` must be given.

- support:

  Optional numeric vector, the same support for every row: returns an
  \\n \times K\\ matrix. Exactly one of `y` or `support` must be given.

- is_zinb:

  Optional override for whether `object` is a `evzinb`-shaped model (a
  zero state present). Defaults to `inherits(object, "evzinb")`, correct
  for a full-sample fitted object; a caller iterating over
  `object$bootstraps` (which does not reliably carry that class –
  round10 H.4/H.5) must pass this explicitly.

## Value

A numeric vector (with `y`) or an \\n \times K\\ numeric matrix (with
`support`).
