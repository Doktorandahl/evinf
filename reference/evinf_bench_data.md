# Synthetic data from an EVZINB model with known true parameters

Generates `n` observations from an EVZINB mixture with fixed,
hand-chosen true parameters (see Details) –
`formula_nb = y ~ x1 + x2 + x3`, `formula_zi`/`formula_evi` sharing the
same three covariates, `formula_pareto` depending on `x1` alone. Meant
for benchmarking (`inst/bench/bench_evinf.R`) and scale testing, not as
a realistic applied example – see [`genevzinb2`](genevzinb2.md) for
that.

## Usage

``` r
evinf_bench_data(n = 1000, seed = NULL)
```

## Arguments

- n:

  Number of observations.

- seed:

  Optional seed; the caller's `.Random.seed` is left untouched either
  way (the same convention as
  [`simulate()`](https://rdrr.io/r/stats/simulate.html)/
  `predict(type = "draws")`).

## Value

A tibble with columns `y`, `x1`, `x2`, `x3`.

## Details

True parameters (on each component's own scale, the same
parameterisation [`coef()`](https://rdrr.io/r/stats/coef.html) returns):

- Zero-state log-odds against the count-state baseline: \\-0.5 + 0.8
  x_1 - 0.5 x_2\\.

- Evi-state log-odds against the count-state baseline: \\-2.5 + 0.5
  x_1 + 0.3 x_3\\ (a deliberately small share, matching an
  extreme-value-inflation process being the rare state).

- Count state (negative binomial, log link): \\\mu = \exp(1.5 + 0.4
  x_1 - 0.3 x_2 + 0.2 x_3)\\, dispersion \\\alpha\_{NB} = 0.8\\.

- Evi state (the discretised Pareto tail the package fits throughout):
  threshold \\C = 50\\, shape \\\alpha\_{PL} = \exp(0.6 + 0.15 x_1)\\
  (kept well above 0 for any simulated `x1`, so the tail never
  collapses).

## Examples

``` r
d <- evinf_bench_data(200, seed = 1)
# \donttest{
m <- evzinb(y ~ x1 + x2 + x3, data = d, bootstrap = FALSE)
#> evinf: using a data-driven candidate range for C_EV: [10, 93]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
# }
```
