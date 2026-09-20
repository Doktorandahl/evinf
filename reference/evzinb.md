# Running an extreme value and zero inflated negative binomial model with bootstrapping

Running an extreme value and zero inflated negative binomial model with
bootstrapping

## Usage

``` r
evzinb(
  formula_nb,
  formula_zi = NULL,
  formula_evi = NULL,
  formula_pareto = NULL,
  data,
  bootstrap = TRUE,
  n_bootstraps = 100,
  multicore = NULL,
  ncores = NULL,
  block = NULL,
  weights = NULL,
  family = evinf_family(),
  boot_seed = NULL,
  control = evinf_control(),
  max.diff.par,
  max.no.em.steps,
  max.no.em.steps.warmup,
  c.lim,
  prune.c.range,
  max.upd.par.zc.multinomial,
  max.upd.par.pl.multinomial,
  max.upd.par.nb,
  max.upd.par.pl,
  no.m.bfgs.steps.multinomial,
  no.m.bfgs.steps.nb,
  no.m.bfgs.steps.pl,
  pdf.pl.type,
  eta.int,
  init.Beta.multinom.ZC,
  init.Beta.multinom.PL,
  init.Beta.NB,
  init.Beta.PL,
  init.Alpha.NB,
  init.C,
  verbose = FALSE
)
```

## Arguments

- formula_nb:

  Formula for the negative binomial (count) component of the model. May
  include an [`offset()`](https://rdrr.io/r/stats/offset.html) term
  (\\\mu\_{NB} = \exp(x'\beta + offset)\\).

- formula_zi:

  Formula for the zero-inflation component of the model. If NULL taken
  as the same formula as nb, with any
  [`offset()`](https://rdrr.io/r/stats/offset.html) term stripped (an
  offset applies only where it is written explicitly, never by
  inheritance). May include its own
  [`offset()`](https://rdrr.io/r/stats/offset.html) term: since the
  zero-inflation logit is the log-odds of the zero state *against* the
  count state, an offset there shifts that log-odds, e.g.
  `offset(log(exposure))` makes a larger exposure relatively less likely
  to land in the structural-zero state.

- formula_evi:

  Formula for the extreme-value inflation component of the model. If
  NULL taken as the same formula as nb, offset stripped as above. May
  include its own [`offset()`](https://rdrr.io/r/stats/offset.html) term
  (e.g. `offset(log(population))` for a probability of an extreme event
  that scales with exposure), read the same way: it shifts the EVI
  log-odds against the count state.

- formula_pareto:

  Formula for the pareto (extreme value) component of the model. If NULL
  taken as the same formula as nb, offset stripped as above.
  [`offset()`](https://rdrr.io/r/stats/offset.html) is **not** supported
  here (errors if present): an offset on a shape parameter has no clear
  reading.

- data:

  data to run the model on

- bootstrap:

  Should bootstrapping be performed. Needed to obtain standard errors
  and p-values

- n_bootstraps:

  Number of bootstraps to run. For use of bootstrapped p-values, at
  least 1,000 bootstraps are recommended. For approximate p-values, a
  lower number can be sufficient

- multicore:

  Bootstrap parallelisation shortcut. The default (`NULL`) respects
  whatever `future` plan is currently set (see the **Parallel
  processing** section). `TRUE` sets a temporary
  [`multisession`](https://future.futureverse.org/reference/multisession.html)
  plan for the duration of the call; `FALSE` forces sequential
  execution. The previous plan is always restored on exit.

- ncores:

  Number of workers when `multicore = TRUE`. Default (`NULL`) is one
  less than the number of available cores. Ignored when `multicore` is
  `NULL` or `FALSE`.

- block:

  Optional case-identifier column for block bootstrapping, given either
  as a bare column name (`block = id`) or a string (`block = "id"`).
  Note that the bundled [`hks`](hks.md) data contain no conflict
  identifier, so the conflict-level cluster bootstrap in Randahl and
  Vegelius (2024) cannot be reproduced from them directly (see
  [`?hks`](hks.md)).

- weights:

  Optional observation weights (round9 D.2), given as a bare column name
  (`weights = wt`), a string naming a column (`weights = "wt"`), or a
  numeric vector. Must be positive and finite; need not be integers
  (analytic weights are allowed, not just frequency counts).
  **Frequency-weight semantics**: every observation's contribution to
  the log-likelihood and to the EM/M-step accumulations is multiplied by
  its weight, and [`nobs()`](https://rdrr.io/r/stats/nobs.html) – and
  therefore `AIC`, `BIC` and the approximate t-based p-values – use
  `sum(weights)`, not the row count (see `sum_weights` in
  [`glance.evzinb`](glance.evzinb.md)). This interpretation of AIC/BIC
  assumes the weights really are frequency weights (repeat-count
  equivalents); for analytic weights the information-criterion values
  are still computed this way but their usual interpretation is weaker.
  **`weights` does not give design-based standard errors for survey
  data** – a sampling weight changes the point estimate, not the
  variance under the sampling design; for that, resample primary
  sampling units with `block =` instead. The bootstrap resamples rows
  exactly as without weights and carries each drawn row's weight along
  (the resampling probabilities themselves are not reweighted).

- family:

  An [`evinf_family()`](evinf_family.md) object, or a string as
  shorthand for its `count` argument (round9 E.0/E.1), e.g.
  `family = "poisson"`. The default reproduces today's negative-binomial
  count state exactly. `count = "poisson"` drops `Alpha.NB` entirely
  (not merely fixes it): it is absent from `par.all`,
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html), and shown
  as absent (not `NA`) in
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`glance()`](https://generics.r-lib.org/reference/glance.html).

- boot_seed:

  Optional bootstrap seed for reproducibility. When supplied it is used
  as-is; when `NULL` a seed is drawn and recorded, so
  `object$boot_seeds` is always populated for a bootstrapped model (see
  [`add_bootstraps`](add_bootstraps.md)).

- control:

  An [`evinf_control()`](evinf_control.md) object holding the EM tuning
  settings (tolerances, candidate range for \\C\_{EV}\\, BFGS steps,
  starting values, ...).

- max.diff.par, max.no.em.steps, max.no.em.steps.warmup, c.lim,
  prune.c.range, max.upd.par.zc.multinomial, max.upd.par.pl.multinomial,
  max.upd.par.nb, max.upd.par.pl, no.m.bfgs.steps.multinomial,
  no.m.bfgs.steps.nb, no.m.bfgs.steps.pl, pdf.pl.type, eta.int,
  init.Beta.multinom.ZC, init.Beta.multinom.PL, init.Beta.NB,
  init.Beta.PL, init.Alpha.NB, init.C:

  **Deprecated.** These EM tuning arguments still work but should be
  passed through `control = evinf_control(...)`; supplying one directly
  overrides the corresponding `control` element and emits a warning. See
  [`evinf_control`](evinf_control.md) for their meaning. Will be removed
  in 0.11.0.

- verbose:

  Logical: should progress of the full run of the model be tracked?

## Value

An object of class 'evzinb'

## Parallel processing

Bootstrap fits (and the per-bootstrap work in
[`add_bootstraps`](add_bootstraps.md), [`lr_test`](lr_test.md),
[`compare_models`](compare_models.md),
[`predict.evzinb`](predict.evzinb.md) and
[`marginal_effects`](marginal_effects.md)) are dispatched with furrr on
top of a future plan. Set the plan once for your session and leave
`multicore` at its default:


    future::plan(future::multisession, workers = 8)
    progressr::handlers(global = TRUE)   # opt in to progress bars
    model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 1000)

Any `future` backend works (`multisession`, `cluster`, `callr`, a HPC
batchtools plan, ...). As a convenience `multicore = TRUE` sets a
temporary `multisession` plan for the single call and restores the
previous plan on exit. Large models may need
`options(future.globals.maxSize = <bytes>)` to raise the default limit
on the data shipped to each worker.

## Reproducibility

The bootstrap **resample indices** (`object$bootstraps[[i]]$boot_id`)
are fully determined by `boot_seed` and are independent of the future
backend, the number of workers, and chunking — a fixed seed always draws
the same resamples.

The **fitted coefficients** are reproducible only to within
floating-point noise. BLAS operations are not bit-reproducible across
processes, so a sequential run and a `multisession` run of the same seed
can differ in the last few digits; and because the EM log-likelihood can
have close local optima, an occasional replicate converges to a
different one under a different backend. For replication material,
record the `future` plan and
[`utils::sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
alongside `boot_seed`, and produce the canonical set of estimates under
a single fixed plan.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
# }

# Peacekeeping and one-sided violence, with an in-formula log1p() transform
# for the Pareto (extreme-value) component (see `?hks`).
if (FALSE) { # \dontrun{
data(hks)
hks_mod <- evzinb(
  osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
    lntpop + brv_AllLag_log + osvAllLagDum + incomp,
  formula_pareto = ~ log1p(troopLag),
  data = hks, n_bootstraps = 5, multicore = FALSE
)
summary(hks_mod)
glance(hks_mod)
predict(hks_mod, type = "harmonic")
} # }
```
