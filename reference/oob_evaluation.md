# Out of bag predictive performance of EVZINB and EVINB models

Out of bag predictive performance of EVZINB and EVINB models

## Usage

``` r
oob_evaluation(
  object,
  predict_type = c("harmonic", "explog"),
  metric = c("rmsle", "rmse", "mse", "mae"),
  exclude_degenerate = TRUE
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

- exclude_degenerate:

  For a single evinf model (or its `$model` slot inside an
  `evzinbcomp`), return `NA` for bootstrap replicates flagged degenerate
  (default `TRUE`) instead of their out-of-bag error, so positions still
  line up with the compared models' replicates; see the `alpha_floor`
  argument of [`evinf_control`](evinf_control.md).

## Value

For a single model, a vector of length `n_bootstraps`. For an
`evzinbcomp` object, a tibble with one column per compared model
(`evinf`, `nb`, `zinb`, ...) and one row per bootstrap. A replicate that
is `NA` in any column (a degenerate evinf replicate, or a compared
model's fit that errored) is `NA` in *every* column, so a column-wise
`na.rm = TRUE` summary (e.g.
`summarize(oob, across(everything(), median, na.rm = TRUE))`) compares
the same set of replicates across models rather than silently different
ones; the row mask itself is available as `attr(out, "excluded")`
(round8 0.4, review §5).

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
oob_evaluation(model)
#> Warning: evinf (harmonic prediction): 2 fitted Pareto alpha values below the floor (0.01); clamping to the floor. Set `alpha_pl_floor` in evinf_control() to change this.
#> Warning: evinf (explog prediction): 2 fitted Pareto alpha values below the floor (0.01); clamping to the floor. Set `alpha_pl_floor` in evinf_control() to change this.
#> [1] 3.325499 2.603393       NA 3.127722       NA
# }
```
