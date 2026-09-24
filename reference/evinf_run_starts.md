# Run em_fit() from one or more starting points and keep the best (round10 G.1)

With `control$n_starts <= 1` (the default), this is a single,
unperturbed call to [`em_fit`](em_fit.md) – bit-identical to before
`n_starts` existed. With `control$n_starts > 1`, the default start
(`ini.val` as given) plus `control$n_starts - 1` perturbed starts
(`Beta.*` jittered by `N(0, start_jitter^2)`, `C` drawn uniformly from
the candidate grid) each run through [`em_fit()`](em_fit.md), in
parallel via [`evinf_pmap`](evinf_pmap.md) seeded from `start_seed`, and
the replicate with the highest final log-likelihood is kept.

## Usage

``` r
evinf_run_starts(
  y,
  x.obj,
  ini.val,
  control,
  model,
  family,
  verbose = FALSE,
  start_seed = NULL
)
```

## Arguments

- y, x.obj, ini.val, control, model, family:

  As for [`em_fit`](em_fit.md).

- verbose:

  Controls [`evinf_pmap()`](evinf_pmap.md)'s progress bar for the
  multi-start pass (each start's own EM console output is always
  suppressed there, regardless of `verbose` – see
  `evinf_em_fit_silent()`). For `n_starts <= 1`, controls the single
  [`em_fit()`](em_fit.md) call's own output, exactly as before.

- start_seed:

  Integer seed for the perturbed starts; drawn if `NULL` and
  `n_starts > 1`.

## Value

A list with `best` (the winning [`em_fit()`](em_fit.md) return list),
`starts` (a tibble with one row per start – `start`, `loglik`, `C`,
`converged`, `n_em_steps`, and list-columns `loglik_trace`/`c_trace`;
`NULL` when `n_starts <= 1`, to keep single-start objects small) and
`start_seed` (`NULL` when `n_starts <= 1`).
