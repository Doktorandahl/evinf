# Likelihood ratio test for individual variables of evzinb

Likelihood ratio test for individual variables of evzinb

## Usage

``` r
lr_test(
  object,
  vars,
  single = TRUE,
  bootstrap = FALSE,
  multicore = NULL,
  ncores = NULL,
  verbose = FALSE
)
```

## Arguments

- object:

  EVZINB or EVINB object to perform likelihood ratio test on

- vars:

  Either a list of character vectors with variable names which to be
  restricted in the LR test or a character vector of variable names. If
  a list, each character vector of the list will be run separately,
  allowing for multiple variables to be restricted as once. If a
  character vector, parameter 'single' can be used to determine whether
  all variables in the vector should be restricted at once (single =
  FALSE) or if the variables should be restricted one by one (single =
  TRUE)

- single:

  Logical. Determining whether variables in 'vars' should be restricted
  individually (single = TRUE) or all at once (single = FALSE)

- bootstrap:

  Should LR tests be conducted on each bootstrapped sample or only on
  the original sample.

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

  Logical. Should the function be verbose?

## Value

A tibble with one row per performed LR test, or, when
`bootstrap = TRUE`, a list with the summary tibble (`results`) and the
per-bootstrap statistics (`boot_results`).

## Details

The likelihood ratio statistic is \\2(\ell\_{full} -
\ell\_{restricted})\\ and is compared to a chi-square distribution with
degrees of freedom equal to the number of design-matrix columns removed
across all model components (so a restricted four-level factor
contributes three degrees of freedom).

## Parallel processing

Bootstrap fits (and the per-bootstrap work in
[`add_bootstraps`](add_bootstraps.md), `lr_test`,
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
#> Error in eval(expr, p) : inv(): matrix is singular
lr_test(model,'x1')
#> Warning: The following arguments of evzinb() are deprecated; pass them through `control = evinf_control()`: max.diff.par, max.no.em.steps, max.no.em.steps.warmup, c.lim, prune.c.range, max.upd.par.zc.multinomial, max.upd.par.pl.multinomial, max.upd.par.nb, max.upd.par.pl, no.m.bfgs.steps.multinomial, no.m.bfgs.steps.nb, no.m.bfgs.steps.pl, pdf.pl.type, eta.int, init.Beta.multinom.ZC, init.Beta.multinom.PL, init.Beta.NB, init.Beta.PL, init.Alpha.NB, init.C.
#> # A tibble: 1 × 6
#>   vars  loglik_full loglik_restricted    df statistic      prob
#>   <chr>       <dbl>             <dbl> <int>     <dbl>     <dbl>
#> 1 x1          -254.             -266.     4      24.5 0.0000641
# }
```
