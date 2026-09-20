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
  changes per component, e.g. `formula_pareto. = . ~ . - x3`. Each
  component is updated against its *own current* formula (the one the
  fitted object actually used, whether the user supplied it or it was
  inherited from `formula_nb` when `NULL`). Changing `formula_nb.` does
  not propagate to a component formula that was originally left `NULL`
  and inherited from it – that component keeps whatever formula it was
  fitted with; pass that component's own `formula_*.` explicitly to
  change it too.

- ...:

  Other arguments of [`evzinb`](evzinb.md) / [`evinb`](evinb.md) to
  change.

- evaluate:

  If `TRUE` (default) re-fit; otherwise return the updated call.

## Value

The updated fit, or the call.

## Details

If `object`'s candidate range for \\C\_{EV}\\ was itself data-driven
(`control$c.lim` was `NULL` at the original fit; see
`object$c_lim_default`), and `data` is one of the arguments being
changed, `control$c.lim` and `control$init.C` are reset to `NULL` so
they are re-resolved from the new data (with the usual message) instead
of silently reusing the range chosen for the old data. A `c.lim` the
user pinned explicitly is always kept as-is.
