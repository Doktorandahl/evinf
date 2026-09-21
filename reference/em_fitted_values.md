# Fitted quantities from a fitted EVZINB / EVINB parameter set

Computes, for every observation, the count-component mean \\\mu\_{NB}\\,
the Pareto shape \\\alpha\_{PL}\\, several summaries of the Pareto tail
(mean, median, \\\exp E\[\log y\]\\, \\1/E\[1/y\]\\) and the
corresponding state-probability-weighted predictions \\\hat y\\. This is
the vectorised equivalent of the per-observation loop that used to sit
at the end of the EM driver; it is identical for the EVZINB and EVINB
models.

## Usage

``` r
em_fitted_values(
  x_obj,
  par,
  props,
  c_ev,
  model = c("evzinb", "evinb"),
  floor = 0.01
)
```

## Arguments

- x_obj:

  List of component design matrices (see
  [`em_extend_design`](em_extend_design.md)); its `offset.nb` element
  (or a zero offset) enters \\\mu\_{NB}\\.

- par:

  List of estimated parameters with elements `Beta.NB` and `Beta.PL`.

- props:

  Numeric matrix \\n \times 3\\ of prior state probabilities (zero,
  count, extreme-value), as returned in `par.mat$Props`.

- c_ev:

  The estimated extreme-value threshold \\C\_{EV}\\.

- model:

  Either `"evzinb"` or `"evinb"`; accepted for interface stability, but
  the computation is currently the same for both.

- floor:

  Floor for the Pareto shape used only when computing the derived tail
  summaries below (round9 0.1, review §2): `exp.E.log.y` and
  `median.pl.vec` involve `exp(1 / alpha.pl.vec)`, which is `Inf` (or
  absurdly large) once `alpha.pl.vec` collapses toward 0. The returned
  `alpha.pl.vec` itself is left unclamped, so callers (e.g.
  [`glance()`](https://generics.r-lib.org/reference/glance.html)'s
  `min_alpha_pl`) still see the true fitted shape. This runs on every EM
  fit, including every bootstrap replicate, so it clamps silently
  (`warn = FALSE`); the collapse is still visible via
  [`glance()`](https://generics.r-lib.org/reference/glance.html)/[`print()`](https://rdrr.io/r/base/print.html)
  without warning on every fit.

## Value

A named list with `mu.nb.vec`, `alpha.pl.vec`, `exp.E.log.y`, `E.inv.y`,
`mean.pl.vec`, `median.pl.vec`, `y.hat.plmedian`, `y.hat.plmean`,
`y.hat.plexpElogy` and `y.hat.pl.E.inv.y`, each a length-\\n\\ numeric
vector.

## See also

[`evzinb()`](evzinb.md), [`evinb()`](evinb.md)
