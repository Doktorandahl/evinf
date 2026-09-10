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
name), `clean` (the row label in the table) and `fmt` (default number of
decimals). Rows: number of observations, number of parameters, alpha_NB,
C_EV, observations at or above C_EV, log-likelihood, AIC, BIC, number of
successful and failed bootstraps, and convergence.

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
#> # A tibble: 11 × 3
#>    raw                 clean                 fmt
#>    <chr>               <chr>               <dbl>
#>  1 nobs                obs                     0
#>  2 npar                par                     0
#>  3 alpha               alpha_nb                2
#>  4 parameter           C_EV                    0
#>  5 n_above_c           obs_above_c_ev          0
#>  6 logLik              logLik                  2
#>  7 aic                 AIC                     1
#>  8 bic                 BIC                     1
#>  9 n_bootstraps        n_bootstraps            0
#> 10 n_failed_bootstraps n_failed_bootstraps     0
#> 11 converged           converged               0
gof_map_evinf(extra = data.frame(raw = "n_em_steps", clean = "EM steps",
                                 fmt = 0))
#> # A tibble: 12 × 3
#>    raw                 clean                 fmt
#>    <chr>               <chr>               <dbl>
#>  1 nobs                obs                     0
#>  2 npar                par                     0
#>  3 alpha               alpha_nb                2
#>  4 parameter           C_EV                    0
#>  5 n_above_c           obs_above_c_ev          0
#>  6 logLik              logLik                  2
#>  7 aic                 AIC                     1
#>  8 bic                 BIC                     1
#>  9 n_bootstraps        n_bootstraps            0
#> 10 n_failed_bootstraps n_failed_bootstraps     0
#> 11 converged           converged               0
#> 12 n_em_steps          EM steps                0

# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Error in eval(expr, p) : inv(): matrix is singular
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
