# Goodness-of-fit map for `modelsummary`

Returns the `gof_map` that
[`modelsummary::modelsummary()`](https://modelsummary.com/man/modelsummary.html)
needs to label and format the goodness-of-fit rows produced by
[`glance.evzinb()`](glance.evzinb.md) /
[`glance.evinb()`](glance.evinb.md). The bundled data object
[`gm_evzinb`](gm_evzinb.md) is simply `gof_map_evinf()`.

## Usage

``` r
gof_map_evinf(extra = NULL)
```

## Arguments

- extra:

  Optional data frame / tibble of additional rows (columns `raw`,
  `clean`, `fmt`) to append, e.g. for a statistic added by a custom
  [`glance()`](https://generics.r-lib.org/reference/glance.html) method.

## Value

A tibble with columns `raw` (the
[`glance()`](https://generics.r-lib.org/reference/glance.html) column
name, unchanged), `clean` (the title-case row label in the table) and
`fmt` (default number of decimals). Rows: observations, parameters,
alpha_NB, C_EV, observations at or above C_EV, log-likelihood, AIC, BIC,
the usable / failed / degenerate bootstrap counts, the number of
bootstrap replicates with C_EV on the candidate-grid boundary, and
convergence.

## Details

The coefficient table produced by
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) with
`component = "all"` has one row per coefficient *per component*.
[`modelsummary::modelsummary()`](https://modelsummary.com/man/modelsummary.html)
cannot render that with its default arguments (the same limitation it
has for
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)); pass
`shape = term + y.level ~ model` together with
`gof_map = gof_map_evinf()`.

## Examples

``` r
gof_map_evinf()
#> # A tibble: 13 × 3
#>    raw                     clean                   fmt
#>    <chr>                   <chr>                 <dbl>
#>  1 nobs                    Observations              0
#>  2 npar                    Parameters                0
#>  3 alpha                   alpha_nb                  2
#>  4 parameter               C_EV                      0
#>  5 n_above_c               Obs. above C_EV           0
#>  6 logLik                  logLik                    2
#>  7 aic                     AIC                       1
#>  8 bic                     BIC                       1
#>  9 n_bootstraps            Usable bootstraps         0
#> 10 n_failed_bootstraps     Failed bootstraps         0
#> 11 n_degenerate_bootstraps Degenerate bootstraps     0
#> 12 n_c_on_boundary         C_EV on boundary          0
#> 13 converged               Converged                 0
gof_map_evinf(extra = data.frame(raw = "n_em_steps", clean = "EM steps",
                                 fmt = 0))
#> # A tibble: 14 × 3
#>    raw                     clean                   fmt
#>    <chr>                   <chr>                 <dbl>
#>  1 nobs                    Observations              0
#>  2 npar                    Parameters                0
#>  3 alpha                   alpha_nb                  2
#>  4 parameter               C_EV                      0
#>  5 n_above_c               Obs. above C_EV           0
#>  6 logLik                  logLik                    2
#>  7 aic                     AIC                       1
#>  8 bic                     BIC                       1
#>  9 n_bootstraps            Usable bootstraps         0
#> 10 n_failed_bootstraps     Failed bootstraps         0
#> 11 n_degenerate_bootstraps Degenerate bootstraps     0
#> 12 n_c_on_boundary         C_EV on boundary          0
#> 13 converged               Converged                 0
#> 14 n_em_steps              EM steps                  0

# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
#> Warning: C_EV reached the boundary of the candidate range in 1 of 5 bootstrap replicates; consider widening c.lim.
if (requireNamespace("modelsummary", quietly = TRUE)) {
  modelsummary::modelsummary(
    model,
    shape = term + y.level ~ model,
    gof_map = gof_map_evinf()
  )
}
#> Error: Package `broom` required for this function to work.
#>   Please install it by running `install.packages("broom")`.
# }
```
