# Evaluate an expression under a temporary `future::plan()`

Evaluate an expression under a temporary
[`future::plan()`](https://future.futureverse.org/reference/plan.html)

## Usage

``` r
evinf_with_plan(multicore = NULL, ncores = NULL, expr)
```

## Arguments

- multicore:

  `NULL` (default) leaves the current
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  untouched; `TRUE` sets
  [`future::multisession`](https://future.futureverse.org/reference/multisession.html)
  with `ncores` workers; `FALSE` sets `sequential`. The plan is always
  restored on exit.

- ncores:

  Number of workers for `multicore = TRUE`; default is one less than
  [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html).

- expr:

  The expression to evaluate (lazily, after the plan is set).

## Value

The value of `expr`.
