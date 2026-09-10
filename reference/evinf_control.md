# Control settings for evzinb() / evinb()

Bundles the ~20 tuning arguments of the EM fitting routine into a single
object, in the style of
[`stats::glm.control()`](https://rdrr.io/r/stats/glm.control.html) or
[`pscl::zeroinfl.control()`](https://rdrr.io/pkg/pscl/man/zeroinfl.control.html).
Pass the result as the `control` argument of [`evzinb()`](evzinb.md) /
[`evinb()`](evinb.md).

## Usage

``` r
evinf_control(
  max.diff.par = 0.01,
  max.no.em.steps = 500,
  max.no.em.steps.warmup = 5,
  c.lim = NULL,
  prune.c.range = FALSE,
  max.upd.par.zc.multinomial = 0.5,
  max.upd.par.pl.multinomial = 0.5,
  max.upd.par.nb = 0.5,
  max.upd.par.pl = 0.5,
  no.m.bfgs.steps.multinomial = 3,
  no.m.bfgs.steps.nb = 3,
  no.m.bfgs.steps.pl = 3,
  pdf.pl.type = c("approx", "exact"),
  eta.int = c(-1, 1),
  init.Beta.multinom.ZC = NULL,
  init.Beta.multinom.PL = NULL,
  init.Beta.NB = NULL,
  init.Beta.PL = NULL,
  init.Alpha.NB = 0.01,
  init.C = NULL,
  alpha_floor = 0.001,
  coef_limit = 50
)
```

## Arguments

- max.diff.par:

  EM convergence tolerance: the algorithm has converged when the maximum
  absolute change in the parameter estimates falls below this. The EM
  log-likelihood can have close local optima; a tighter tolerance
  reduces (without eliminating) the chance that two future backends run
  with the same `boot_seed` settle in different ones (see the
  **Reproducibility** section of `?`[`evzinb`](evzinb.md)).

- max.no.em.steps:

  Maximum number of EM steps.

- max.no.em.steps.warmup:

  Number of EM steps in each warm-up round.

- c.lim:

  `NULL` or a numeric vector of length 2. The candidate set for
  \\C\_{EV}\\ is the unique observed response values within this range.
  `NULL` (the default) uses a data-driven range (see Details).

- prune.c.range:

  `FALSE`, or a number in \\0, 1\]: thin the candidate set to about
  `length(c.lim) * (1 - prune.c.range)` values.

- max.upd.par.zc.multinomial, max.upd.par.pl.multinomial,
  max.upd.par.nb, max.upd.par.pl:

  Maximum parameter-change step sizes for the zero-inflation,
  extreme-value inflation, count and Pareto components.

- no.m.bfgs.steps.multinomial, no.m.bfgs.steps.nb, no.m.bfgs.steps.pl:

  Number of BFGS steps per M-step for the multinomial, count and Pareto
  blocks.

- pdf.pl.type:

  Pareto density approximation: `"approx"` or `"exact"`.

- eta.int:

  Interval for the eta line search, a numeric vector of length 2.

- init.Beta.multinom.ZC, init.Beta.multinom.PL, init.Beta.NB,
  init.Beta.PL:

  Optional starting values for the component coefficient vectors (`NULL`
  starts at zero).

- init.Alpha.NB:

  Starting value for the negative-binomial dispersion.

- init.C:

  `NULL` or a starting value for \\C\_{EV}\\ within `c.lim`. `NULL` uses
  the median of the candidate set.

- alpha_floor, coef_limit:

  Thresholds for flagging a bootstrap replicate as *degenerate*
  (`$degenerate`, `$degenerate_reason`), so it is excluded from
  bootstrap summaries by default (see `exclude_degenerate`) and counted
  in [`glance()`](https://generics.r-lib.org/reference/glance.html) /
  [`failed_bootstraps()`](failed_bootstraps.md). A replicate is
  degenerate if (in this order): its EM did not converge; any fitted
  coefficient in a linear predictor (`Beta.*`) or the negative-binomial
  dispersion is non-finite or exceeds `coef_limit` in absolute value; or
  its smallest fitted Pareto shape is non-finite or below `alpha_floor`.
  `alpha_floor` (default `0.001`) is deliberately low — it catches an
  outright tail collapse (shape near 0), not a merely heavy tail.
  `coef_limit` (default `50`) is on the linear-predictor scale, where
  `50` is already extreme.

## Value

A list of class `"evinf_control"`.

## Details

When `c.lim = NULL`, [`evzinb()`](evzinb.md) / [`evinb()`](evinb.md)
choose a range from the data: the lower bound is the 90th percentile of
the positive response values (rounded down to an observed value), the
upper bound is the third-largest unique response value; if fewer than 10
unique positive values lie in that range the lower bound drops to the
75th percentile. The chosen range is printed with a message and stored
in `object$control$c.lim`.

## Examples

``` r
evinf_control(max.no.em.steps = 300, c.lim = c(20, 500))
#> <evinf_control>
#>   max.diff.par             0.01
#>   max.no.em.steps          300
#>   max.no.em.steps.warmup   5
#>   c.lim                    20, 500
#>   prune.c.range            FALSE
#>   pdf.pl.type              approx
#>   init.Alpha.NB            0.01
#>   init.C                   NULL (data-driven)
```
