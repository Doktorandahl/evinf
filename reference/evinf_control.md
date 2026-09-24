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
  coef_limit = 50,
  max.c.iter = 50,
  alpha_pl_floor = 0.01,
  n_starts = 1L,
  start_jitter = 0.5,
  chunk_size = NULL
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

  `FALSE`, or a number in \\0, 1): thin the candidate set to about
  `length(c.lim) * (1 - prune.c.range)` values.

- max.upd.par.zc.multinomial, max.upd.par.pl.multinomial,
  max.upd.par.nb, max.upd.par.pl:

  Maximum parameter-change step sizes for the zero-inflation,
  extreme-value inflation, count and Pareto components.

- no.m.bfgs.steps.multinomial, no.m.bfgs.steps.nb, no.m.bfgs.steps.pl:

  Number of BFGS steps per M-step for the multinomial, count and Pareto
  blocks.

- pdf.pl.type:

  Which Pareto-block derivatives the M-step uses for \\\beta\_{PL}\\:
  `"approx"` (the default) uses the gradient/Hessian of the *continuous*
  Pareto log-density; `"exact"` uses the gradient/Hessian of the
  *discretised* Pareto log-pmf that the likelihood itself always uses
  (this only changes the Newton step taken each M-step, not what is
  being maximised). The two typically converge to nearly the same
  estimates; `"exact"` can help when `"approx"`'s steps are poorly
  scaled for a heavy-tailed fit.

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

- max.c.iter:

  Maximum number of outer ECME iterations (each one re-profiles
  \\C\_{EV}\\ over the candidate grid) per phase (warm-up, convergence).
  Guards against the profile oscillating between two candidate values
  forever. If the convergence-phase loop hits this cap, the fit's
  `converge` is set to `FALSE` (see `$c_converged` to tell this apart
  from the EM inner loop not converging) and, for a full-sample fit, a
  [`warning()`](https://rdrr.io/r/base/warning.html) names the last two
  \\C\_{EV}\\ values visited. The warm-up phase is capped independently
  and recorded in `$c_warmup_capped` (with its own
  [`warning()`](https://rdrr.io/r/base/warning.html) for a full-sample
  fit); it does not affect `converge`, since warm-up is a short
  exploratory phase and not settling there is not by itself a sign the
  fit failed.

- alpha_pl_floor:

  Floor for the fitted Pareto shape \\\alpha\_{PL}\\ (round9 0.1). Some
  observations' fitted \\\alpha\_{PL}\\ can collapse toward 0 (heaviest
  possible tail); several downstream quantities involve
  \\\exp(1/\alpha\_{PL})\\ or \\1/\alpha\_{PL}\\ and silently return
  `Inf` or an astronomically large finite number when that happens
  (`predict(type = "harmonic")`, `predict(type = "explog")`, the
  continuous-Pareto mixture quantile, and the `$fitted` tail summaries).
  Any fitted \\\alpha\_{PL}\\ below `alpha_pl_floor` is clamped to it
  before these calculations, with a
  [`warning()`](https://rdrr.io/r/base/warning.html) naming how many
  observations were clamped; the unclamped values are still what
  [`glance()`](https://generics.r-lib.org/reference/glance.html)'s
  `min_alpha_pl` and [`print()`](https://rdrr.io/r/base/print.html)'s
  note report, so a collapsed EV shape stays visible.

- n_starts:

  Number of starting points for the full-sample fit (round10 G.1, audit
  §5.10). `1` (the default) fits once from the default start, exactly as
  before this argument existed. With `n_starts > 1`, that default start
  plus `n_starts - 1` perturbed starts each run through the EM (in
  parallel via [`evinf_pmap`](evinf_pmap.md), seeded from `start_seed`),
  and the replicate with the highest final log-likelihood is kept – the
  EM can land in different local optima depending on floating-point
  accumulation order, and multiple starts are the principled response.
  Every replicate's outcome is recorded in `object$starts`; bootstrap
  replicates warm-start from the winning solution, same as they would
  from a single-start fit.

- start_jitter:

  Standard deviation of the `N(0, start_jitter^2)` perturbation applied
  to each component's starting coefficients for the `n_starts - 1`
  perturbed starts; their starting \\C\_{EV}\\ is drawn uniformly from
  the candidate grid instead. Unused when `n_starts <= 1`.

- chunk_size:

  Number of bootstrap replicates (or starts) handed to a worker at a
  time in [`evinf_pmap`](evinf_pmap.md) (round10 J.2, audit §5.11).
  `NULL` (the default) picks
  `max(1, ceiling(B / (4 * future::nbrOfWorkers())))` at dispatch time
  (`B` the number of replicates being dispatched) – four chunks per
  worker, balancing per-task scheduling overhead (worse for `B` far
  larger than the worker count and `chunk_size = 1`) against a straggler
  chunk leaving workers idle near the end (worse for a `chunk_size`
  close to `B / nbrOfWorkers()`). Chunk size only changes how replicates
  are grouped for dispatch, never which `boot_id`s they draw (each
  replicate's L'Ecuyer stream comes from its position in the full
  sequence, not from the chunk it lands in), so bootstrap output is
  identical across chunk sizes for the same seed. Set explicitly to
  override, e.g. `1L` to restore the pre-J.2 default.

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
