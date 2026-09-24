# Run a function over a list in parallel, with a progress bar

The single point through which every parallel loop in the package goes.
Wraps
[`furrr::future_map()`](https://furrr.futureverse.org/reference/future_map.html)
with reproducible per-element RNG streams (`seed`) and a progressr
progress bar that updates once per element.

## Usage

``` r
evinf_pmap(
  .x,
  .f,
  ...,
  seed,
  label = "bootstrap",
  verbose = FALSE,
  chunk_size = NULL
)
```

## Arguments

- .x:

  A list or vector to iterate over.

- .f:

  A function applied to each element of `.x` (plus `...`).

- ...:

  Additional arguments passed on to `.f` on every call.

- seed:

  An integer seed. `furrr` derives an independent L'Ecuyer stream for
  each element from it, identical regardless of the number of workers or
  the chunk size.

- label:

  Progress-bar message.

- verbose:

  If `TRUE`, always show a progress bar. Otherwise one is shown only
  when the user has enabled a global progressr handler.

- chunk_size:

  Number of elements handed to a worker at a time. `NULL` (the default,
  round10 J.2, audit §5.11) picks
  `max(1, ceiling(length(.x) / (4 * future::nbrOfWorkers())))` – four
  chunks per worker, balancing per-task scheduling overhead (a
  `chunk_size` of 1 when `length(.x)` is far larger than the worker
  count) against a straggler chunk leaving workers idle near the end (a
  `chunk_size` close to `length(.x) / nbrOfWorkers()`). This only
  changes how elements are grouped for dispatch, never each element's
  L'Ecuyer stream (verified empirically: derived from the element's
  position in `.x`, not the chunk it lands in), so results are identical
  across chunk sizes for the same `seed`.

## Value

A list, the element-wise results of `.f`.
