# Coefficients of an evzinb / evinb model

Coefficients of an evzinb / evinb model

## Usage

``` r
# S3 method for class 'evzinb'
coef(object, component = "all", ...)

# S3 method for class 'evinb'
coef(object, component = "all", ...)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` model.

- component:

  One of `"all"` (the default), `"count"`, `"zero"`, `"evi"`,
  `"pareto"`. Deprecated aliases `"nb"`, `"zi"`, `"evinf"` are accepted.

- ...:

  Unused.

## Value

For `"all"`, a named numeric vector `<component>_<term>`, `alpha_nb`,
`c_ev`; for a single component the plain named vector.
