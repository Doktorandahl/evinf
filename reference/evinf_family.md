# Model family for evzinb() / evinb() (round9 E.0)

Specifies the distributional family of the count state (`count`) and of
the zero state (`zero`). One object, so later additions (a
truncated-count state, a geometric count) do not add more arguments to
[`evzinb()`](evzinb.md) / [`evinb()`](evinb.md).

## Usage

``` r
evinf_family(count = c("nbinom", "poisson"), zero = c("mixture", "hurdle"))
```

## Arguments

- count:

  `"nbinom"` (the default; a negative-binomial count state) or
  `"poisson"` (a Poisson count state – `Alpha.NB` is dropped entirely:
  not in `par.all`, not in
  [`coef()`](https://rdrr.io/r/stats/coef.html)/
  [`vcov()`](https://rdrr.io/r/stats/vcov.html)/[`confint()`](https://rdrr.io/r/stats/confint.html)/[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
  shown as absent (not `NA`) in
  [`summary()`](https://rdrr.io/r/base/summary.html)/[`glance()`](https://generics.r-lib.org/reference/glance.html)).

- zero:

  `"mixture"` (the default; the existing zero-inflation mixture, where a
  zero can come from either the zero state or the count state) or
  `"hurdle"` (the zero state owns every zero and the count state is
  zero-truncated; [`evzinb()`](evzinb.md) only – [`evinb()`](evinb.md)
  has no zero state to hurdle over and errors if asked for one).

## Value

An object of class `"evinf_family"`: a list with elements `count` and
`zero`.

## Examples

``` r
evinf_family()
#> <evinf_family>
#>   count: nbinom
#>   zero:  mixture
evinf_family(count = "poisson")
#> <evinf_family>
#>   count: poisson
#>   zero:  mixture
```
