# Stacked tidy / glance for a model comparison

Stacked tidy / glance for a model comparison

## Usage

``` r
# S3 method for class 'evzinbcomp'
tidy(x, ...)

# S3 method for class 'evzinbcomp'
glance(x, ...)
```

## Arguments

- x:

  An `evzinbcomp` object.

- ...:

  Passed to the per-model
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html) /
  [`glance()`](https://generics.r-lib.org/reference/glance.html).

## Value

A tibble with a leading `model` column.
