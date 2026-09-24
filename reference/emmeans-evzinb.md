# emmeans support for evzinb / evinb models

On load, if emmeans is installed, the package registers
`recover_data.evzinb()`/`emm_basis.evzinb()` (and the `evinb`
equivalents) with emmeans via
[`emmeans::.emm_register()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.html),
so `emmeans::emmeans(model, ...)` works for the count component's linear
predictor (a log link; `type = "response"`/`regrid()` back-transforms it
the usual emmeans way). Standard errors come from the count-component
block of the bootstrap covariance matrix
([`vcov.evzinb`](vcov.evzinb.md)).

## Details

Only the count component is supported (`component = "count"`, the
default and, for now, the only accepted value) – see the source comments
in `R/emmeans_evzinb.R` for why a state-probability (zero/count/evi)
multinomial basis is not implemented.
`predict(type = "states", confint = TRUE)` covers that case instead,
with bootstrap rather than delta-method intervals.
