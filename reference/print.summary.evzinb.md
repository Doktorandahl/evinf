# Print methods for evzinb / evinb summaries

Print methods for evzinb / evinb summaries

## Usage

``` r
# S3 method for class 'summary.evzinb'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)

# S3 method for class 'summary.evinb'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- x:

  A `summary.evzinb` or `summary.evinb` object.

- digits:

  Number of significant digits for the coefficient tables.

- signif.stars:

  Logical; show significance stars.

- ...:

  Not used.

## Value

`x`, invisibly.
