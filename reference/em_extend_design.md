# Prepend the intercept column to each component design matrix

The C++ estimation routines (`log_lik_fun()`, `update_bfgs_fun()`)
expect each component's design matrix to include a leading column of
ones. This helper turns the raw `x_obj` (design matrices without
intercept, as built by `evinf_design()`) into the extended matrices,
supplying a ones-only matrix for any component with no covariates and a
zero offset when none was given.

## Usage

``` r
em_extend_design(x_obj, n)
```

## Arguments

- x_obj:

  List with elements `X.multinom.ZC`, `X.multinom.PL`, `X.NB`, `X.PL`
  (each a numeric matrix \\n \times p\\ without an intercept column, or
  `NULL`) and optionally `offset.nb`, `offset.zc`, `offset.pl_mult`
  (round9 D.1) and `weights` (round9 D.2).

- n:

  Number of observations.

## Value

A list with elements `zc`, `pl_mult`, `nb`, `pl` (the extended numeric
matrices), `offset`, `offset_zc`, `offset_pl_mult` (length-`n` numeric,
zero when not supplied) and `weights` (length-`n` numeric, one when not
supplied).

## See also

[`evzinb()`](evzinb.md), [`evinb()`](evinb.md)
