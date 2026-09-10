# Function to compare evzinb or evinb models with zinb and nb models

Function to compare evzinb or evinb models with zinb and nb models

## Usage

``` r
compare_models(
  object,
  nb_comparison = TRUE,
  zinb_comparison = TRUE,
  winsorize = FALSE,
  razorize = FALSE,
  cutoff_value = 10,
  init_theta = NULL,
  multicore = NULL,
  ncores = NULL
)
```

## Arguments

- object:

  A fitted evzinb or evinb model object

- nb_comparison:

  Should comparison be made with a negative binomial model?

- zinb_comparison:

  Should comparisons be made with the zinb model? Not available for
  `evinb` objects (there is no zero-inflation component); it defaults to
  `FALSE` for those, with a message, and errors if set to `TRUE`
  explicitly.

- winsorize:

  Should winsorizing be done in the comparisons?

- razorize:

  Should razorizing (trimming) be done in the comparisons?

- cutoff_value:

  Integer: Which observation should be used as a basis for
  winsorizing/razorising. E.g. 10 means that everything larger than the
  10th observation will be winsorized/razorised

- init_theta:

  Optional initial value for theta in the NB specification

- multicore:

  Bootstrap parallelisation shortcut. The default (`NULL`) respects
  whatever `future` plan is currently set (see the **Parallel
  processing** section). `TRUE` sets a temporary
  [`multisession`](https://future.futureverse.org/reference/multisession.html)
  plan for the duration of the call; `FALSE` forces sequential
  execution. The previous plan is always restored on exit.

- ncores:

  Number of workers when `multicore = TRUE`. Default (`NULL`) is one
  less than the number of available cores. Ignored when `multicore` is
  `NULL` or `FALSE`.

## Value

An object of class `evzinbcomp`: a list whose first element `model` is
the original evzinb/evinb model (also available as `evzinb` for
backwards compatibility), followed by the compared `nb` / `zinb` models
(and their winsorized/razorized variants when requested).

## Parallel processing

Bootstrap fits (and the per-bootstrap work in
[`add_bootstraps`](add_bootstraps.md), [`lr_test`](lr_test.md),
`compare_models`, [`predict.evzinb`](predict.evzinb.md) and
[`marginal_effects`](marginal_effects.md)) are dispatched with furrr on
top of a future plan. Set the plan once for your session and leave
`multicore` at its default:


    future::plan(future::multisession, workers = 8)
    progressr::handlers(global = TRUE)   # opt in to progress bars
    model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 1000)

Any `future` backend works (`multisession`, `cluster`, `callr`, a HPC
batchtools plan, ...). As a convenience `multicore = TRUE` sets a
temporary `multisession` plan for the single call and restores the
previous plan on exit. Bootstrap draws are reproducible from `boot_seed`
and do not depend on the number of workers. Large models may need
`options(future.globals.maxSize = <bytes>)` to raise the default limit
on the data shipped to each worker.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
compare_models(model)
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: iteration limit reached
#> Warning: NaNs produced
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: alternation limit reached
#> 
#>  Model comparison of  evzinb 
#>   Compared models:  nb, zinb 
#>   Number of compared models:  2 
#>  Number of bootstraps: 5 
# }

if (FALSE) { # \dontrun{
data(hks)
hks_mod <- evzinb(osvAll ~ troopLag + lntpop + brv_AllLag_log + osvAllLagDum,
                  formula_pareto = ~ log1p(troopLag),
                  data = hks, n_bootstraps = 5, multicore = FALSE)
cmp <- compare_models(hks_mod)
compare_fit(cmp)
} # }
```
