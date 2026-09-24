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

A tibble with columns `id`, `type` (`"error"`, `"degenerate"` or
`"not_converged"` – round10 G.3: a replicate that ran without erroring
and isn't degenerate, but whose `converge` is `FALSE`, was previously
invisible here), `message` (the error text, the degeneracy reason, or
`NA` for `"not_converged"` – see `n_em_steps`/`c_converged`/
`c_warmup_capped` there instead) and, for `"not_converged"` rows,
`n_em_steps`, `c_converged`, `c_warmup_capped` (`NA` for
`"error"`/`"degenerate"` rows); zero rows when every replicate is usable
and converged.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
failed_bootstraps(model)
#> # A tibble: 3 × 6
#>   id          type       message          n_em_steps c_converged c_warmup_capped
#>   <chr>       <chr>      <chr>                 <int> <lgl>       <lgl>          
#> 1 bootstrap_3 degenerate smallest fitted…         NA NA          NA             
#> 2 bootstrap_4 degenerate smallest fitted…         NA NA          NA             
#> 3 bootstrap_5 degenerate coefficient (In…         NA NA          NA             
# }
```
