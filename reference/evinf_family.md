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

## Details

**Sign convention against `pscl`** (round10 0.9, review §7): for
`zero = "hurdle"`,
[`hurdle`](https://rdrr.io/pkg/pscl/man/hurdle.html)'s zero-hurdle
component models \\P(Y \> 0)\\, while evinf's zero-state coefficients
model \\P(\text{zero state})\\ directly – so a fitted
`evzinb(family = evinf_family(zero = "hurdle"))` model's zero
coefficients (`Beta.multinom.ZC`) are the exact negatives of the
corresponding `pscl::hurdle(..., dist = "negbin")` zero coefficients,
everything else about the fit being equal (verified in `test-hurdle.R`).
This is a clean reparameterisation because a hurdle's zero state owns
every zero outright. `zero = "mixture"` (the default) has no such simple
relationship to
[`zeroinfl`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)'s
zero-inflation coefficients: a mixture's zeros can also arise from the
count state, so the two models' zero components are not a
reparameterisation of each other.

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
