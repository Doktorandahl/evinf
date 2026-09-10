# The bootstrap replicates usable for a summary

Drops the `try-error` replicates and, unless
`exclude_degenerate = FALSE`, the ones flagged degenerate.

## Usage

``` r
evinf_usable_bootstraps(object, exclude_degenerate = TRUE)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` model.

- exclude_degenerate:

  Also drop replicates flagged degenerate (default `TRUE`).

## Value

A (possibly empty) list of bootstrap fits.
