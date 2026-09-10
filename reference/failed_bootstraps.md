# Inspect the bootstrap replicates that failed or came out degenerate

Bootstrap fits that error are kept as `try-error` objects rather than
discarded; fits that converge but whose extreme-value tail is so heavy
that summaries from them are effectively unbounded are flagged
*degenerate* (see `alpha_floor` in [`evinf_control`](evinf_control.md)).
This returns both, which are useful diagnostics (e.g. a resample with a
singular design, or one whose Pareto shape collapsed).

## Usage

``` r
failed_bootstraps(object)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` model.

## Value

A tibble with columns `id`, `type` (`"error"` or `"degenerate"`) and
`message` (the error text, or the reason the replicate is degenerate);
zero rows when every replicate is usable.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
failed_bootstraps(model)
#> # A tibble: 4 × 3
#>   id          type       message                                                
#>   <chr>       <chr>      <chr>                                                  
#> 1 bootstrap_1 degenerate smallest fitted Pareto shape 4.49e-11 < alpha_floor (0…
#> 2 bootstrap_2 degenerate smallest fitted Pareto shape 3.76e-42 < alpha_floor (0…
#> 3 bootstrap_4 degenerate smallest fitted Pareto shape 7.92e-12 < alpha_floor (0…
#> 4 bootstrap_5 degenerate smallest fitted Pareto shape 6.04e-13 < alpha_floor (0…
# }
```
