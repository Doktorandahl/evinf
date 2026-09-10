# Add more bootstrap replicates to a fitted model

Fits `n` additional bootstrap replicates with the same block structure
and appends them, so you don't have to refit everything when the
original number turns out to be too few.

## Usage

``` r
add_bootstraps(
  object,
  n,
  boot_seed = NULL,
  multicore = NULL,
  ncores = NULL,
  verbose = FALSE
)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` model (with bootstraps).

- n:

  Number of additional bootstrap replicates.

- boot_seed:

  RNG seed for this batch. Each batch needs its own seed: the bootstrap
  stream is fully determined by the seed, so reusing one that the model
  (or an earlier `add_bootstraps()` call) already used would silently
  duplicate those draws and is an error. When `NULL` a fresh seed is
  drawn and recorded in `object$boot_seeds`.

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

- verbose:

  Show a progress bar.

## Value

The model with `n` more replicates in `object$bootstraps` (names
continue `bootstrap_<k>`).

## Parallel processing

Bootstrap fits (and the per-bootstrap work in `add_bootstraps`,
[`lr_test`](lr_test.md), [`compare_models`](compare_models.md),
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
previous plan on exit. Bootstrap draws are reproducible from `boot_seed`
and do not depend on the number of workers. Large models may need
`options(future.globals.maxSize = <bytes>)` to raise the default limit
on the data shipped to each worker.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
model <- add_bootstraps(model, 5)
#> Error in eval(expr, p) : inv(): matrix is singular
#> Error in eval(expr, p) : inv(): matrix is singular
# }
```
