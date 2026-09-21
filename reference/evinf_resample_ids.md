# Bootstrap resample indices under a chosen scheme

The single generator behind `bootstrap_scheme = ` in
[`evzinb()`](evzinb.md) / [`evinb()`](evinb.md) /
[`add_bootstraps()`](add_bootstraps.md) (round9 F, audit §5.7).

## Usage

``` r
evinf_resample_ids(
  n,
  scheme = c("iid", "cluster", "moving_block", "stationary"),
  block_vec = NULL,
  time_vec = NULL,
  block_length = NULL
)
```

## Arguments

- n:

  Number of rows in the estimation data.

- scheme:

  One of `"iid"` (plain row resampling), `"cluster"` (resample whole
  units named by `block_vec`, keeping every row of a drawn unit – what
  `block = ` has always done), `"moving_block"` or `"stationary"`
  (block-resample each unit's own time series, see Kunsch 1989 / Politis
  and Romano 1994).

- block_vec:

  Optional unit identifier, one value per row. Required for `"cluster"`.
  For `"moving_block"` / `"stationary"`, treated as a single unit (the
  whole data) when `NULL`.

- time_vec:

  Required for `"moving_block"` / `"stationary"`: a numeric/orderable
  time index, one value per row, that must already be strictly
  increasing within every unit's rows as they appear in the data – this
  function never sorts it for you (see Details).

- block_length:

  Block length \\L\\ for `"moving_block"` / `"stationary"`. `NULL` (the
  default) uses \\\lceil T^{1/3} \rceil\\ for each unit's own length
  \\T\\ (a message names the value(s) used); a user-supplied
  `block_length` that exceeds some unit's \\T\\ is an error rather than
  a silent clamp.

## Value

An integer vector of length `n`: row indices into the estimation data,
with replacement, resampled under `scheme`.

## Details

With overlapping blocks (`"moving_block"` / `"stationary"`) the set of
rows never drawn in a given replicate is both smaller and more
temporally correlated than under i.i.d. resampling, so an out-of-bag
error estimated from it (see [`oob_evaluation`](oob_evaluation.md)) is
optimistic relative to genuine forecasting performance.
