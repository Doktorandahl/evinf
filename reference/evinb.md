# Running an extreme value inflated negative binomial model with bootstrapping

Running an extreme value inflated negative binomial model with
bootstrapping

## Usage

``` r
evinb(
  formula_nb,
  formula_evi = NULL,
  formula_pareto = NULL,
  data,
  bootstrap = TRUE,
  n_bootstraps = 100,
  multicore = NULL,
  ncores = NULL,
  block = NULL,
  boot_seed = NULL,
  control = evinf_control(),
  max.diff.par,
  max.no.em.steps,
  max.no.em.steps.warmup,
  c.lim,
  prune.c.range,
  max.upd.par.pl.multinomial,
  max.upd.par.nb,
  max.upd.par.pl,
  no.m.bfgs.steps.multinomial,
  no.m.bfgs.steps.nb,
  no.m.bfgs.steps.pl,
  pdf.pl.type,
  eta.int,
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

  Formula for the negative binomial (count) component of the model

- formula_evi:

  Formula for the extreme-value inflation component of the model. If
  NULL taken as the same formula as nb

- formula_pareto:

  Formula for the pareto (extreme value) component of the model. If NULL
  taken as the same formula as nb

- data:

  Data to run the model on

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

- boot_seed:

  Optional bootstrap seed for reproducibility. When supplied it is used
  as-is; when `NULL` a seed is drawn and recorded, so
  `object$boot_seeds` is always populated for a bootstrapped model (see
  [`add_bootstraps`](add_bootstraps.md)).

- control:

  An [`evinf_control()`](evinf_control.md) object holding the EM tuning
  settings.

- max.diff.par, max.no.em.steps, max.no.em.steps.warmup, c.lim,
  prune.c.range, max.upd.par.pl.multinomial, max.upd.par.nb,
  max.upd.par.pl, no.m.bfgs.steps.multinomial, no.m.bfgs.steps.nb,
  no.m.bfgs.steps.pl, pdf.pl.type, eta.int, init.Beta.multinom.PL,
  init.Beta.NB, init.Beta.PL, init.Alpha.NB, init.C:

  **Deprecated.** Pass these through `control = evinf_control(...)`;
  supplying one directly overrides the corresponding `control` element
  and emits a warning. See [`evinf_control`](evinf_control.md).

- verbose:

  Should progress be printed for the first run of evinb

## Value

An object of class 'evinb'

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
model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
# }

if (FALSE) { # \dontrun{
data(hks)
hks_mod <- evinb(
  osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
    lntpop + brv_AllLag_log + osvAllLagDum + incomp,
  formula_pareto = ~ log1p(troopLag),
  data = hks, n_bootstraps = 5, multicore = FALSE
)
} # }
```
