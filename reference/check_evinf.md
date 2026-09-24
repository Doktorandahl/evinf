# Diagnostic report for a fitted evzinb / evinb model

A single report gathering the convergence and health checks otherwise
spread across [`print()`](https://rdrr.io/r/base/print.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html) and
[`failed_bootstraps()`](failed_bootstraps.md): EM and
\\C\_{EV}\\-profile convergence, whether \\C\_{EV}\\ landed on a
candidate-grid boundary (full sample and across bootstrap replicates),
the smallest fitted Pareto shape against `alpha_pl_floor`, failed and
degenerate bootstrap replicates, start agreement (see `n_starts` in
[`evinf_control`](evinf_control.md)), and, for the block bootstrap
schemes, the out-of-bag row fraction.

## Usage

``` r
check_evinf(object)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` model.

## Value

A classed tibble (`"evinf_check"`) with columns `check`, `status`
(`"ok"`, `"note"` or `"warning"`) and `detail` (a one-line description,
with a suggestion for anything not `"ok"`).

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
check_evinf(model)
#> evinf diagnostic report
#>   6 checks, 0 warnings, 1 note
#> 
#> [ok]     EM converged                                      
#>          The EM algorithm converged.
#> [ok]     C_EV profile converged (convergence phase)        
#>          The C_EV profile settled within max.c.iter.
#> [ok]     C_EV profile converged (warm-up phase)            
#>          The warm-up phase settled within max.c.iter.
#> [ok]     C_EV on candidate-grid boundary                   
#>          The fitted C_EV lies inside the candidate range.
#> [ok]     Fitted Pareto shape (alpha_pl) not collapsed      
#>          min(alpha_pl) = 0.019, above the floor.
#> [note]   Bootstrap replicates usable                       
#>          3 usable of 5 (0 failed, 2 degenerate) -- see failed_bootstraps() for the details.
# }
```
