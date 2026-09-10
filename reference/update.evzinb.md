# Update and re-fit an evzinb / evinb model

Update and re-fit an evzinb / evinb model

## Usage

``` r
# S3 method for class 'evzinb'
update(
  object,
  formula_nb.,
  formula_zi.,
  formula_evi.,
  formula_pareto.,
  ...,
  evaluate = TRUE
)

# S3 method for class 'evinb'
update(
  object,
  formula_nb.,
  formula_evi.,
  formula_pareto.,
  ...,
  evaluate = TRUE
)
```

## Arguments

- object:

  A fitted model (must carry `object$call`).

- formula_nb., formula_zi., formula_evi., formula_pareto.:

  Optional
  [`update.formula`](https://rdrr.io/r/stats/update.formula.html)-style
  changes per component, e.g. `formula_pareto. = . ~ . - x3`.

- ...:

  Other arguments of [`evzinb`](evzinb.md) / [`evinb`](evinb.md) to
  change.

- evaluate:

  If `TRUE` (default) re-fit; otherwise return the updated call.

## Value

The updated fit, or the call.
