# `marginaleffects` support for evzinb / evinb models

On load the package registers the `"evzinb"` / `"evinb"` classes with
marginaleffects and provides `get_coef()`, `set_coef()`, `get_vcov()`
and `get_predict()` methods, so
[`marginaleffects::avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html),
`predictions()` etc. work. Standard errors come from the delta method
applied to the bootstrap covariance matrix
([`vcov.evzinb`](vcov.evzinb.md)).

## Details

The coefficient vector includes `c_ev`, the extreme-value threshold,
which is estimated on the grid of unique observed response values rather
than by a smooth optimiser. marginaleffects perturbs every coefficient
continuously when it builds the delta-method Jacobian, so standard
errors for quantities that depend strongly on `c_ev` should be treated
as approximate. The recommended route for uncertainty on covariate
effects is the bootstrap-based
[`marginal_effects()`](marginal_effects.md).
