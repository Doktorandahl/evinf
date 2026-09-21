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
convergence (both of the EM loop and of the C_EV profile).

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
#> 14 c_converged             C_EV profile settled      0
gof_map_evinf(extra = data.frame(raw = "n_em_steps", clean = "EM steps",
                                 fmt = 0))
#> # A tibble: 15 × 3
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
#> 14 c_converged             C_EV profile settled      0
#> 15 n_em_steps              EM steps                  0

# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
if (requireNamespace("modelsummary", quietly = TRUE) &&
    requireNamespace("broom", quietly = TRUE)) {
  modelsummary::modelsummary(
    model,
    shape = term + y.level ~ model,
    gof_map = gof_map_evinf()
  )
}
#> 
#> +-----------------------+---------+---------+
#> |                       | y.level | (1)     |
#> +=======================+=========+=========+
#> | (Intercept)           | zero    | 0.511   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.459) |
#> +-----------------------+---------+---------+
#> |                       | evi     | -1.729  |
#> +-----------------------+---------+---------+
#> |                       |         | (1.607) |
#> +-----------------------+---------+---------+
#> |                       | count   | 3.131   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.153) |
#> +-----------------------+---------+---------+
#> |                       | pareto  | 2.714   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.194) |
#> +-----------------------+---------+---------+
#> | x1                    | zero    | -0.674  |
#> +-----------------------+---------+---------+
#> |                       |         | (0.390) |
#> +-----------------------+---------+---------+
#> |                       | evi     | 0.736   |
#> +-----------------------+---------+---------+
#> |                       |         | (1.463) |
#> +-----------------------+---------+---------+
#> |                       | count   | 0.843   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.752) |
#> +-----------------------+---------+---------+
#> |                       | pareto  | -2.278  |
#> +-----------------------+---------+---------+
#> |                       |         | (0.011) |
#> +-----------------------+---------+---------+
#> | x2                    | zero    | 0.607   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.663) |
#> +-----------------------+---------+---------+
#> |                       | evi     | -0.429  |
#> +-----------------------+---------+---------+
#> |                       |         | (1.625) |
#> +-----------------------+---------+---------+
#> |                       | count   | 0.599   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.212) |
#> +-----------------------+---------+---------+
#> |                       | pareto  | 1.449   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.097) |
#> +-----------------------+---------+---------+
#> | x3                    | zero    | -0.420  |
#> +-----------------------+---------+---------+
#> |                       |         | (0.270) |
#> +-----------------------+---------+---------+
#> |                       | evi     | 0.166   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.897) |
#> +-----------------------+---------+---------+
#> |                       | count   | 0.241   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.209) |
#> +-----------------------+---------+---------+
#> |                       | pareto  | 1.503   |
#> +-----------------------+---------+---------+
#> |                       |         | (0.233) |
#> +-----------------------+---------+---------+
#> | Observations          |         | 100     |
#> +-----------------------+---------+---------+
#> | Parameters            |         | 18      |
#> +-----------------------+---------+---------+
#> | alpha_nb              |         | 1.43    |
#> +-----------------------+---------+---------+
#> | C_EV                  |         | 184     |
#> +-----------------------+---------+---------+
#> | Obs. above C_EV       |         | 10      |
#> +-----------------------+---------+---------+
#> | logLik                |         | -254.03 |
#> +-----------------------+---------+---------+
#> | AIC                   |         | 544.1   |
#> +-----------------------+---------+---------+
#> | BIC                   |         | 591.0   |
#> +-----------------------+---------+---------+
#> | Usable bootstraps     |         | 2       |
#> +-----------------------+---------+---------+
#> | Failed bootstraps     |         | 0       |
#> +-----------------------+---------+---------+
#> | Degenerate bootstraps |         | 3       |
#> +-----------------------+---------+---------+
#> | C_EV on boundary      |         | 1       |
#> +-----------------------+---------+---------+
#> | Converged             |         | TRUE    |
#> +-----------------------+---------+---------+
#> | C_EV profile settled  |         | TRUE    |
#> +-----------------------+---------+---------+ 
# }
```
