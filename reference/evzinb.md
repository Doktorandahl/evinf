# Running an extreme value and zero inflated negative binomial model with bootstrapping

Running an extreme value and zero inflated negative binomial model with
bootstrapping

## Usage

``` r
evzinb(
  formula_nb,
  formula_zi = NULL,
  formula_evi = NULL,
  formula_pareto = NULL,
  data,
  bootstrap = TRUE,
  n_bootstraps = 100,
  multicore = NULL,
  ncores = NULL,
  block = NULL,
  weights = NULL,
  time = NULL,
  bootstrap_scheme = NULL,
  block_length = NULL,
  family = evinf_family(),
  boot_seed = NULL,
  start_seed = NULL,
  control = evinf_control(),
  max.diff.par,
  max.no.em.steps,
  max.no.em.steps.warmup,
  c.lim,
  prune.c.range,
  max.upd.par.zc.multinomial,
  max.upd.par.pl.multinomial,
  max.upd.par.nb,
  max.upd.par.pl,
  no.m.bfgs.steps.multinomial,
  no.m.bfgs.steps.nb,
  no.m.bfgs.steps.pl,
  pdf.pl.type,
  eta.int,
  init.Beta.multinom.ZC,
  init.Beta.multinom.PL,
  init.Beta.NB,
  init.Beta.PL,
  init.Alpha.NB,
  init.C,
  verbose = FALSE
)
```

## Arguments

- formula_nb:

  Formula for the negative binomial (count) component of the model. May
  include an [`offset()`](https://rdrr.io/r/stats/offset.html) term
  (\\\mu\_{NB} = \exp(x'\beta + offset)\\).

- formula_zi:

  Formula for the zero-inflation component of the model. If NULL taken
  as the same formula as nb, with any
  [`offset()`](https://rdrr.io/r/stats/offset.html) term stripped (an
  offset applies only where it is written explicitly, never by
  inheritance). May include its own
  [`offset()`](https://rdrr.io/r/stats/offset.html) term: since the
  zero-inflation logit is the log-odds of the zero state *against* the
  count state, an offset there shifts that log-odds, e.g.
  `offset(log(exposure))` makes a larger exposure relatively less likely
  to land in the structural-zero state.

- formula_evi:

  Formula for the extreme-value inflation component of the model. If
  NULL taken as the same formula as nb, offset stripped as above. May
  include its own [`offset()`](https://rdrr.io/r/stats/offset.html) term
  (e.g. `offset(log(population))` for a probability of an extreme event
  that scales with exposure), read the same way: it shifts the EVI
  log-odds against the count state.

- formula_pareto:

  Formula for the pareto (extreme value) component of the model. If NULL
  taken as the same formula as nb, offset stripped as above.
  [`offset()`](https://rdrr.io/r/stats/offset.html) is **not** supported
  here (errors if present): an offset on a shape parameter has no clear
  reading.

- data:

  data to run the model on

- bootstrap:

  Should bootstrapping be performed. Needed to obtain standard errors
  and p-values

- n_bootstraps:

  Number of bootstraps to run. For use of bootstrapped p-values, at
  least 1,000 bootstraps are recommended. For approximate p-values, a
  lower number can be sufficient

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

- block:

  Optional case-identifier column for block bootstrapping, given either
  as a bare column name (`block = id`) or a string (`block = "id"`).
  Note that the bundled [`hks`](hks.md) data contain no conflict
  identifier, so the conflict-level cluster bootstrap in Randahl and
  Vegelius (2024) cannot be reproduced from them directly (see
  [`?hks`](hks.md)).

- weights:

  Optional observation weights (round9 D.2), given as a bare column name
  (`weights = wt`), a string naming a column (`weights = "wt"`), or a
  numeric vector. Must be positive and finite; need not be integers
  (analytic weights are allowed, not just frequency counts).
  **Frequency-weight semantics**: every observation's contribution to
  the log-likelihood and to the EM/M-step accumulations is multiplied by
  its weight, and [`nobs()`](https://rdrr.io/r/stats/nobs.html) – and
  therefore `AIC`, `BIC` and the approximate t-based p-values – use
  `sum(weights)`, not the row count (see `sum_weights` in
  [`glance.evzinb`](glance.evzinb.md)). This interpretation of AIC/BIC
  assumes the weights really are frequency weights (repeat-count
  equivalents); for analytic weights the information-criterion values
  are still computed this way but their usual interpretation is weaker.
  **`weights` does not give design-based standard errors for survey
  data** – a sampling weight changes the point estimate, not the
  variance under the sampling design; for that, resample primary
  sampling units with `block =` instead. The bootstrap resamples rows
  exactly as without weights and carries each drawn row's weight along
  (the resampling probabilities themselves are not reweighted).

- time:

  Optional time index for panel/time-series bootstrap resampling (round9
  F), given as a bare column name (`time = t`) or a string
  (`time = "t"`); required for
  `bootstrap_scheme %in% c("moving_block", "stationary")`. Must be
  strictly increasing within every `block` unit's rows as they already
  appear in the data – this is never sorted for you; sort `data` by
  `(block, time)` first if it isn't already.

- bootstrap_scheme:

  One of `"iid"` (plain row resampling, the default when `block` is not
  given), `"cluster"` (resample whole `block` units, the default when
  `block` is given – what `block =` has always done), `"moving_block"`
  or `"stationary"` (block-resample each unit's own time series; see
  Kunsch 1989 / Politis and Romano 1994). The two block schemes need
  `time`; `block` is optional for them (the whole data is treated as one
  unit when omitted). With overlapping blocks the out-of-bag set is
  smaller and more temporally correlated than under iid resampling, so
  out-of-bag error ([`oob_evaluation`](oob_evaluation.md)) is optimistic
  relative to genuine forecasting performance; each bootstrap
  replicate's realised out-of-bag fraction is stored as `$oob_fraction`
  next to `$boot_id`.

- block_length:

  Block length for
  `bootstrap_scheme %in% c("moving_block", "stationary")`; `NULL` (the
  default) uses `ceiling(T^(1/3))` for each unit's own length \\T\\ (a
  message names the value(s) used, once, at the original fit – not on
  every bootstrap replicate).

- family:

  An [`evinf_family()`](evinf_family.md) object, or a string as
  shorthand for its `count` argument (round9 E.0/E.1/E.2), e.g.
  `family = "poisson"`. The default reproduces today's
  negative-binomial, mixture-zero model exactly. `count = "poisson"`
  drops `Alpha.NB` entirely (not merely fixes it): it is absent from
  `par.all`, [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html), and shown
  as absent (not `NA`) in
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`glance()`](https://generics.r-lib.org/reference/glance.html).
  `zero = "hurdle"` makes the zero state own every zero (rather than
  competing with the count state for them) and zero-truncates the count
  state; verified to match
  [`pscl::hurdle()`](https://rdrr.io/pkg/pscl/man/hurdle.html)'s
  coefficients and log-likelihood to numerical precision when the
  extreme-value state is unreachable.

- boot_seed:

  Optional bootstrap seed for reproducibility. When supplied it is used
  as-is; when `NULL` a seed is drawn and recorded, so
  `object$boot_seeds` is always populated for a bootstrapped model (see
  [`add_bootstraps`](add_bootstraps.md)).

- start_seed:

  Optional seed for the perturbed starts when `control$n_starts > 1`
  (round10 G.1); recorded as `object$start_seed` (`NULL` for the default
  `n_starts = 1`, same as `boot_seed` for an unbootstrapped fit). Unused
  otherwise.

- control:

  An [`evinf_control()`](evinf_control.md) object holding the EM tuning
  settings (tolerances, candidate range for \\C\_{EV}\\, BFGS steps,
  starting values, ...).

- max.diff.par, max.no.em.steps, max.no.em.steps.warmup, c.lim,
  prune.c.range, max.upd.par.zc.multinomial, max.upd.par.pl.multinomial,
  max.upd.par.nb, max.upd.par.pl, no.m.bfgs.steps.multinomial,
  no.m.bfgs.steps.nb, no.m.bfgs.steps.pl, pdf.pl.type, eta.int,
  init.Beta.multinom.ZC, init.Beta.multinom.PL, init.Beta.NB,
  init.Beta.PL, init.Alpha.NB, init.C:

  **Deprecated.** These EM tuning arguments still work but should be
  passed through `control = evinf_control(...)`; supplying one directly
  overrides the corresponding `control` element and emits a warning. See
  [`evinf_control`](evinf_control.md) for their meaning. Will be removed
  in a future release.

- verbose:

  Logical: should progress of the full run of the model be tracked?

## Value

An object of class 'evzinb'

## Parallel processing

Bootstrap fits (and the per-bootstrap work in
[`add_bootstraps`](add_bootstraps.md), [`lr_test`](lr_test.md),
[`compare_models`](compare_models.md),
[`predict.evzinb`](predict.evzinb.md) and
[`marginal_effects`](marginal_effects.md)) are dispatched with furrr on
top of a future plan. Set the plan once for your session and leave
`multicore` at its default:


    future::plan(future::multisession, workers = 8)
    progressr::handlers(global = TRUE)   # opt in to progress bars
    model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 1000)

Any `future` backend works (`multisession`, `cluster`, `callr`, a HPC
batchtools plan, ...). As a convenience `multicore = TRUE` sets a
temporary `multisession` plan for the single call and restores the
previous plan on exit. Large models may need
`options(future.globals.maxSize = <bytes>)` to raise the default limit
on the data shipped to each worker.

## Reproducibility

The bootstrap **resample indices** (`object$bootstraps[[i]]$boot_id`)
are fully determined by `boot_seed` and are independent of the future
backend, the number of workers, and chunking — a fixed seed always draws
the same resamples.

The **fitted coefficients** are reproducible only to within
floating-point noise. BLAS operations are not bit-reproducible across
processes, so a sequential run and a `multisession` run of the same seed
can differ in the last few digits; and because the EM log-likelihood can
have close local optima, an occasional replicate converges to a
different one under a different backend. For replication material,
record the `future` plan and
[`utils::sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
alongside `boot_seed`, and produce the canonical set of estimates under
a single fixed plan.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
#> Warning: C_EV equalled the lower endpoint (173) in 1 of 5 bootstrap replicates; consider widening c.lim.
# }

# Peacekeeping and one-sided violence, with an in-formula log1p() transform
# for the Pareto (extreme-value) component (see `?hks`).
# \donttest{
data(hks)
hks_mod <- evzinb(
  osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
    lntpop + brv_AllLag_log + osvAllLagDum + incomp,
  formula_pareto = ~ log1p(troopLag),
  data = hks, n_bootstraps = 2, multicore = FALSE
)
#> evinf: using a data-driven candidate range for C_EV: [209, 8586]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
summary(hks_mod)
#> EVZINB model summary
#> ====================
#> 
#> Count component (negative binomial)
#>                       Estimate Std. Error approx t Pr(boot)
#> (Intercept)           2.602327                           <1
#> troopLag             -0.150160                           <1
#> policeLag            -2.661202                           <1
#> militaryobserversLag  5.586607                           <1
#> epduration           -0.002241                           <1
#> lntpop               -0.115017                           <1
#> brv_AllLag_log        0.036009                           <1
#> osvAllLagDum          0.392168                           <1
#> incomp                0.965755                           <1
#> 
#> Zero-inflation component
#>                       Estimate Std. Error approx t Pr(boot)
#> (Intercept)           9.318688                           <1
#> troopLag              0.167433                           <1
#> policeLag             0.045592                           <1
#> militaryobserversLag -3.973772                           <1
#> epduration           -0.005826                           <1
#> lntpop               -0.527963                           <1
#> brv_AllLag_log       -0.376696                           <1
#> osvAllLagDum         -3.292438                           <1
#> incomp               -0.946799                           <1
#> 
#> Extreme-value inflation component
#>                        Estimate Std. Error approx t Pr(boot)
#> (Intercept)          -23.999287                           <1
#> troopLag               0.205188                           <1
#> policeLag             -4.440042                           <1
#> militaryobserversLag  -5.952385                           <1
#> epduration            -0.006611                           <1
#> lntpop                 1.478051                           <1
#> brv_AllLag_log         0.043702                           <1
#> osvAllLagDum          -0.035886                           <1
#> incomp                 3.668050                           <1
#> 
#> Pareto (extreme value) component
#>                 Estimate Std. Error approx t Pr(boot)
#> (Intercept)     -0.05172                           <1
#> log1p(troopLag) -0.86289                           <1
#> 
#> ----------------------------------------
#> alpha_NB: 2.671   C_EV: 209
#> Observations at or above C_EV: 89
#> Mean state proportions:  zero = 0.712   count = 0.272   evi = 0.015
#> Observations: 3746   Parameters: 31   df: 3715
#> logLik: -5580   AIC: 11220   BIC: 11410   Converged: TRUE
#> Bootstraps: failed = 0, degenerate = 1
glance(hks_mod)
#> # A tibble: 1 × 25
#>    nobs sum_weights  npar family  alpha parameter    aic    bic logLik converged
#>   <int>       <dbl> <int> <chr>   <dbl>     <dbl>  <dbl>  <dbl>  <dbl> <lgl>    
#> 1  3746        3746    31 nbinom…  2.67       209 11221. 11414. -5580. TRUE     
#> # ℹ 15 more variables: c_converged <lgl>, n_above_c <int>, n_em_steps <int>,
#> #   min_alpha_pl <dbl>, n_bootstraps <int>, n_failed_bootstraps <int>,
#> #   n_degenerate_bootstraps <int>, n_c_on_boundary <int>,
#> #   oob_fraction_mean <dbl>, oob_fraction_min <dbl>, oob_fraction_max <dbl>,
#> #   n_starts <int>, n_starts_at_best <int>, median_boot_em_steps <int>,
#> #   n_boot_c_capped <int>
predict(hks_mod, type = "harmonic")
#>    [1] 1.554030e+02 1.952972e+02 2.114395e+02 1.840180e+02 2.041893e+02
#>    [6] 3.846727e+01 3.842307e+01 3.837882e+01 3.833451e+01 3.829015e+01
#>   [11] 3.824573e+01 4.020209e+01 4.015423e+01 4.010630e+01 4.005832e+01
#>   [16] 4.001028e+01 3.996217e+01 3.991400e+01 3.986578e+01 3.981749e+01
#>   [21] 3.976914e+01 3.972073e+01 3.967226e+01 3.749866e+01 3.745374e+01
#>   [26] 3.740877e+01 3.736374e+01 3.731865e+01 3.727351e+01 3.722831e+01
#>   [31] 3.718305e+01 3.713774e+01 1.391059e+02 3.704693e+01 3.700145e+01
#>   [36] 1.368354e+01 2.557662e+01 2.700476e+01 1.362299e+01 2.653452e+01
#>   [41] 1.112325e+00 1.114058e+00 1.115796e+00 1.117539e+00 1.119287e+00
#>   [46] 1.121041e+00 1.157917e+00 1.159668e+00 1.161425e+00 1.163186e+00
#>   [51] 1.164953e+00 1.166725e+00 1.168502e+00 1.170284e+00 1.172071e+00
#>   [56] 1.173864e+00 1.175661e+00 1.177464e+00 1.141372e+00 1.143199e+00
#>   [61] 1.145031e+00 1.146868e+00 1.148710e+00 1.150557e+00 1.152409e+00
#>   [66] 1.154267e+00 1.156129e+00 1.301312e+01 1.159869e+00 1.161747e+00
#>   [71] 1.387652e+02 1.383017e+02 1.378388e+02 1.373766e+02 1.369149e+02
#>   [76] 1.364538e+02 1.359933e+02 1.355335e+02 1.350743e+02 1.346158e+02
#>   [81] 1.341579e+02 1.826180e+02 1.911390e+02 1.856825e+02 1.713547e+02
#>   [86] 1.617796e+02 1.714213e+02 1.339839e+02 1.335222e+02 1.330612e+02
#>   [91] 3.195469e+01 1.321412e+02 3.188455e+01 1.338551e+02 1.333896e+02
#>   [96] 1.329248e+02 1.324607e+02 1.319974e+02 1.315349e+02 1.310730e+02
#>  [101] 1.306120e+02 1.714245e+02 1.659612e+02 1.292335e+02 1.731289e+02
#>  [106] 1.675138e+02 1.652026e+02 1.301319e+02 1.702778e+02 1.554227e+02
#>  [111] 1.482625e+02 1.496392e+02 1.594530e+02 1.616176e+02 1.371742e+02
#>  [116] 1.373354e+02 1.389461e+02 1.482139e+02 1.273095e+02 1.411031e+02
#>  [121] 1.364987e+02 7.471393e+01 1.418707e+02 1.601854e+02 1.638107e+02
#>  [126] 1.897787e+02 1.510926e+02 1.644469e+02 7.258885e+01 1.558283e+02
#>  [131] 1.290745e+02 5.208489e+01 1.596858e+02 5.767271e+01 2.083220e+02
#>  [136] 2.146828e+02 2.335475e+02 8.615309e+01 8.619965e+01 1.086119e+02
#>  [141] 2.777508e+02 2.724060e+02 2.587171e+02 2.540182e+02 2.601313e+02
#>  [146] 2.586085e+02 2.602907e+02 2.369701e+02 7.948005e+01 2.273179e+02
#>  [151] 2.275590e+02 2.084598e+02 2.117994e+02 2.544575e+02 2.580806e+02
#>  [156] 2.456230e+02 2.862529e+02 2.952640e+02 3.073525e+02 2.549320e+02
#>  [161] 2.577740e+02 2.256960e+02 7.176223e+01 1.979498e+02 6.690589e+01
#>  [166] 1.833238e+01 1.422933e+01 1.462344e+01 1.682660e+01 1.644004e+01
#>  [171] 1.342842e+01 1.358336e+01 1.370113e+01 1.457735e+01 1.514551e+01
#>  [176] 1.860724e+01 1.787384e+01 1.376803e+01 1.512309e+01 1.586856e+01
#>  [181] 1.606310e+01 1.576549e+01 1.573600e+01 1.738933e+01 2.552842e+01
#>  [186] 2.271316e+01 2.474117e+01 2.160152e+01 2.332181e+00 2.337780e+00
#>  [191] 2.343396e+00 2.349029e+00 2.354678e+00 2.360343e+00 2.366025e+00
#>  [196] 2.371724e+00 2.377439e+00 2.383171e+00 1.093882e+01 2.443381e+00
#>  [201] 2.449190e+00 2.455016e+00 2.460858e+00 1.103522e+01 2.472593e+00
#>  [206] 2.478485e+00 2.484394e+00 2.490320e+00 2.496263e+00 2.502222e+00
#>  [211] 2.508198e+00 2.509155e+00 2.515162e+00 2.521185e+00 2.527225e+00
#>  [216] 2.533282e+00 2.539356e+00 2.545446e+00 2.551553e+00 2.557677e+00
#>  [221] 2.563818e+00 2.569975e+00 2.576149e+00 2.613998e+00 2.620223e+00
#>  [226] 2.626464e+00 2.632723e+00 2.638998e+00 2.645290e+00 2.651599e+00
#>  [231] 2.657924e+00 2.664267e+00 6.722908e+00 2.677002e+00 7.616376e+00
#>  [236] 2.722081e+00 3.556980e+00 5.358421e+00 4.841396e+01 4.053370e+01
#>  [241] 7.931525e+00 4.183711e+01 4.284273e+01 3.467510e+01 4.025624e+01
#>  [246] 4.888679e+01 4.647405e+01 5.007043e+01 3.465122e+01 5.097986e+01
#>  [251] 4.869363e+01 4.579549e+01 4.955391e+01 5.004318e+01 5.125729e+01
#>  [256] 3.441334e+01 5.137425e+01 4.783263e+01 4.879462e+00 4.522516e+01
#>  [261] 4.660993e+01 2.960065e+00 5.108725e+01 5.139535e+01 5.346856e+01
#>  [266] 4.505281e+01 3.414617e+01 4.876653e+01 5.051540e+01 4.621470e+01
#>  [271] 5.329585e+01 4.559587e+01 5.454814e+01 4.028757e+01 6.031272e+00
#>  [276] 5.112527e+01 8.263325e+00 1.300966e+01 4.708955e+01 4.710366e+01
#>  [281] 1.243104e+01 3.384762e+01 4.055425e+01 4.572063e+01 3.383328e+01
#>  [286] 3.245926e+00 7.407440e+00 7.658650e+00 4.089441e+01 3.362225e+01
#>  [291] 4.289171e+01 4.749673e+01 4.976325e+01 4.576368e+01 4.777378e+01
#>  [296] 1.106360e+01 4.908473e+01 4.766772e+01 4.599199e+01 5.156079e+01
#>  [301] 4.719753e+01 4.650780e+01 4.688302e+01 4.824688e+01 4.979848e+01
#>  [306] 4.396757e+01 1.337082e+01 4.670285e+01 1.398014e+01 1.276885e+01
#>  [311] 1.366606e+01 1.410867e+01 4.579714e+01 4.303170e+01 4.614024e+01
#>  [316] 4.570570e+01 4.559861e+01 4.516870e+01 4.696243e+01 4.519970e+01
#>  [321] 4.306565e+01 4.498690e+01 4.342727e+01 1.507084e+01 4.277837e+01
#>  [326] 1.113543e+01 4.763206e+01 4.532082e+01 4.172563e+01 4.106374e+01
#>  [331] 1.409539e+01 4.177860e+01 4.327135e+01 4.209446e+01 4.546633e+01
#>  [336] 4.296702e+01 3.196277e+01 4.233182e+01 4.216421e+01 9.529971e+00
#>  [341] 3.925216e+01 3.937756e+01 4.237256e+01 4.245200e+01 4.299137e+01
#>  [346] 3.948184e+01 1.148425e+01 1.448761e+01 4.026823e+01 4.022137e+00
#>  [351] 8.121634e+00 5.634636e+01 2.126605e+01 1.074025e+01 1.195396e+01
#>  [356] 7.952234e+00 1.205367e+01 4.356423e+00 3.519585e+01 9.857721e+00
#>  [361] 4.530642e+01 4.191584e+01 4.061901e+01 1.085760e+01 3.733568e+01
#>  [366] 4.427182e+01 4.517218e+01 4.448824e+01 1.500765e+01 1.271494e+01
#>  [371] 1.017896e+01 1.043600e+01 1.228485e+01 6.199280e+00 1.337694e+01
#>  [376] 6.215301e+00 5.008391e+00 4.806823e+00 5.245640e+00 4.417853e+00
#>  [381] 4.284056e+00 4.292204e+00 4.393284e+00 4.433281e+00 4.441734e+00
#>  [386] 4.450192e+00 4.458653e+00 4.518511e+00 4.527107e+00 4.601733e+00
#>  [391] 4.610484e+00 2.390241e+00 2.395873e+00 2.401521e+00 2.407186e+00
#>  [396] 2.412868e+00 2.418566e+00 2.424280e+00 2.430011e+00 2.435759e+00
#>  [401] 4.550253e+01 2.447304e+00 2.495096e+00 6.454690e+01 3.582462e+01
#>  [406] 2.512653e+00 2.518539e+00 2.524441e+00 9.095828e+00 6.833146e+00
#>  [411] 5.210159e+01 2.548218e+00 2.554204e+00 2.560207e+00 2.607227e+00
#>  [416] 3.576886e+01 2.619355e+00 3.571319e+01 3.568509e+01 3.445709e+00
#>  [421] 2.643814e+00 2.649970e+00 2.656143e+00 2.662334e+00 2.668541e+00
#>  [426] 2.674765e+00 2.716087e+00 6.864733e+00 2.728649e+00 2.734955e+00
#>  [431] 2.741279e+00 2.747618e+00 2.753975e+00 2.760349e+00 3.915898e+01
#>  [436] 2.773146e+00 2.779570e+00 2.786011e+00 2.829412e+00 2.835901e+00
#>  [441] 2.842406e+00 2.848929e+00 3.528561e+01 3.525122e+01 2.868596e+00
#>  [446] 2.875186e+00 3.514716e+01 3.511218e+01 2.895054e+00 2.901710e+00
#>  [451] 2.963879e+00 2.970588e+00 2.977313e+00 2.984054e+00 2.990812e+00
#>  [456] 2.997587e+00 3.004378e+00 3.011186e+00 3.018010e+00 3.024851e+00
#>  [461] 3.031708e+00 3.038581e+00 3.103701e+00 3.110624e+00 3.117562e+00
#>  [466] 3.489177e+01 3.131489e+00 3.138476e+00 3.476945e+01 3.152500e+00
#>  [471] 3.159536e+00 3.166588e+00 4.967902e+01 3.180740e+00 4.005800e+01
#>  [476] 3.374516e+00 3.506581e+01 4.865817e+01 3.396041e+00 3.403248e+00
#>  [481] 3.410469e+00 3.417707e+00 3.424960e+00 3.432228e+00 3.439511e+00
#>  [486] 3.446810e+00 3.514330e+00 3.521661e+00 3.529007e+00 1.119106e+01
#>  [491] 3.543744e+00 3.551134e+00 3.558540e+00 3.565961e+00 3.573396e+00
#>  [496] 3.580845e+00 3.588310e+00 1.289234e+01 1.007764e+01 7.650234e+00
#>  [501] 8.242740e+00 3.688004e+00 6.562015e+00 3.703113e+00 3.710688e+00
#>  [506] 9.236800e+00 3.725880e+00 3.733497e+00 5.555676e+00 3.748772e+00
#>  [511] 1.224891e+01 3.982682e+00 3.990324e+00 1.020411e+01 4.005646e+00
#>  [516] 4.013326e+00 4.021018e+00 1.043105e+01 4.036439e+00 4.044168e+00
#>  [521] 4.051910e+00 4.059663e+00 4.302726e+00 4.310409e+00 4.318103e+00
#>  [526] 4.325808e+00 4.333524e+00 8.019370e+00 4.348987e+00 4.356734e+00
#>  [531] 4.364492e+00 4.372260e+00 2.337732e+01 4.387825e+00 4.507391e+00
#>  [536] 4.515136e+00 4.522890e+00 4.530653e+00 4.538425e+00 4.546206e+00
#>  [541] 4.553996e+00 4.561794e+00 4.569601e+00 4.577416e+00 4.585240e+00
#>  [546] 4.593072e+00 4.711910e+00 4.719681e+00 3.318433e+01 3.312009e+01
#>  [551] 3.305590e+01 3.299176e+01 3.292767e+01 4.766458e+00 4.774278e+00
#>  [556] 4.782104e+00 4.789936e+00 4.797775e+00 4.913592e+00 4.921354e+00
#>  [561] 4.929121e+00 4.936893e+00 4.944670e+00 4.952452e+00 4.960238e+00
#>  [566] 4.968029e+00 4.975824e+00 3.217387e+01 4.991427e+00 4.999234e+00
#>  [571] 3.217041e+01 5.117463e+00 3.203800e+01 4.437979e+01 4.414731e+01
#>  [576] 3.183996e+01 3.882258e+01 5.163817e+00 5.171551e+00 4.312350e+01
#>  [581] 4.162108e+01 3.884092e+01 4.246401e+01 8.863273e+00 1.470361e+01
#>  [586] 3.135146e+01 5.331920e+00 5.339553e+00 5.347186e+00 5.354818e+00
#>  [591] 5.463232e+00 5.471085e+00 5.478938e+00 1.331298e+01 2.882411e+00
#>  [596] 2.888183e+00 2.893972e+00 2.899778e+00 2.905601e+00 2.911440e+00
#>  [601] 2.917297e+00 2.923170e+00 7.648738e+00 5.793277e+00 2.940892e+00
#>  [606] 1.350291e+01 1.414316e+01 2.145369e+01 2.964440e+00 2.970448e+00
#>  [611] 2.976474e+00 2.982517e+00 2.988576e+00 2.994653e+00 3.000746e+00
#>  [616] 3.006857e+00 3.012984e+00 3.212184e+00 3.218347e+00 3.224527e+00
#>  [621] 3.230724e+00 3.236938e+00 3.243169e+00 3.249417e+00 3.255681e+00
#>  [626] 3.261962e+00 3.268260e+00 3.274575e+00 3.280906e+00 3.486091e+00
#>  [631] 3.492432e+00 3.498789e+00 3.505162e+00 3.511552e+00 3.517959e+00
#>  [636] 3.524382e+00 3.530821e+00 3.537277e+00 3.543749e+00 3.550238e+00
#>  [641] 3.556743e+00 2.516924e+01 1.282141e+01 5.821026e+01 1.281321e+01
#>  [646] 2.717137e+01 1.280503e+01 1.280094e+01 1.279686e+01 8.219888e+01
#>  [651] 8.193859e+01 8.167881e+01 1.335261e+01 1.143314e+02 8.291385e+01
#>  [656] 8.264709e+01 8.238088e+01 8.211523e+01 1.103202e+02 1.068974e+02
#>  [661] 1.331360e+01 1.330873e+01 8.079529e+01 1.329898e+01 8.223347e+01
#>  [666] 8.196336e+01 8.709167e+01 8.142493e+01 8.115661e+01 8.088888e+01
#>  [671] 8.062175e+01 8.035523e+01 3.605270e+01 1.382625e+01 1.382049e+01
#>  [676] 7.929520e+01 2.748399e+01 8.065326e+01 8.547786e+01 1.438303e+01
#>  [681] 7.984224e+01 1.436963e+01 1.436292e+01 7.903701e+01 1.434949e+01
#>  [686] 1.434276e+01 7.823760e+01 1.432930e+01 1.275466e+01 9.809618e+01
#>  [691] 7.235105e+01 7.211173e+01 7.187300e+01 1.273339e+01 7.139734e+01
#>  [696] 7.116040e+01 9.313619e+01 7.068834e+01 1.271205e+01 7.021868e+01
#>  [701] 2.961113e+01 7.770524e+01 9.540587e+01 7.046047e+01 8.795566e+01
#>  [706] 8.935735e+01 6.975032e+01 2.345870e+01 6.928003e+01 3.405266e+01
#>  [711] 8.401430e+01 8.205418e+01 4.208908e+01 9.993397e+01 9.060720e+01
#>  [716] 9.964897e+01 1.933023e+01 6.665427e+01 1.929945e+01 4.505670e+01
#>  [721] 9.556458e+01 4.848556e+01 3.837022e+01 2.667487e+01 9.135683e+01
#>  [726] 9.783850e+01 4.864078e+01 1.361634e+01 8.677289e+01 9.211324e+01
#>  [731] 9.358476e+01 6.150411e+01 8.347049e+01 6.398664e+01 5.045311e+01
#>  [736] 5.093006e+01 4.646614e+01 8.967177e+01 8.277619e+01 5.133419e+01
#>  [741] 5.631980e+01 5.147455e+01 8.080578e+01 5.450064e+01 2.422107e+01
#>  [746] 1.973676e+01 6.341962e+01 1.386860e+01 3.417738e+01 6.373399e+01
#>  [751] 6.350727e+01 6.850498e+01 5.443450e+01 8.599150e+01 6.463217e+01
#>  [756] 5.567095e+01 1.419153e+01 5.285497e+01 2.002020e+01 5.831668e+01
#>  [761] 7.681567e+01 9.005451e+01 4.098128e+01 1.552422e+01 8.538891e+01
#>  [766] 8.864182e+01 5.982859e+01 7.897141e+01 8.037566e+01 5.381328e+01
#>  [771] 5.925502e+01 8.000469e+01 5.635317e+01 8.509509e+01 6.733697e+01
#>  [776] 8.363998e+01 6.211079e+01 6.188025e+01 5.351541e+01 1.588743e+01
#>  [781] 1.587465e+01 6.421512e+01 6.073931e+01 1.583620e+01 6.149396e+01
#>  [786] 7.551719e+01 1.635908e+01 6.080418e+01 8.080942e+01 6.034834e+01
#>  [791] 6.012163e+01 7.363352e+01 7.820677e+01 5.944628e+01 7.374327e+01
#>  [796] 7.180520e+01 8.440504e+01 8.460988e+01 8.307248e+01 8.070597e+01
#>  [801] 5.941623e+01 5.919000e+01 7.200025e+01 7.147339e+01 6.264130e+01
#>  [806] 6.922471e+01 6.393038e+01 7.149332e+01 6.382306e+01 5.850553e+01
#>  [811] 1.730642e+01 5.700949e+01 6.321262e+01 5.905327e+01 3.357726e+01
#>  [816] 6.098753e+01 6.304659e+01 7.676133e+01 5.790556e+01 7.206852e+01
#>  [821] 9.579137e+01 1.325100e+02 9.041301e+01 7.037459e+01 6.738115e+01
#>  [826] 5.630454e+01 5.518971e+01 6.387847e+01 5.779560e+01 6.726274e+01
#>  [831] 5.836561e+01 3.825636e+01 3.722944e+01 3.153872e+01 3.297944e+01
#>  [836] 3.489483e+01 3.941093e+01 3.646609e+01 3.368133e+01 3.707692e+01
#>  [841] 2.662534e+01 3.669954e+01 2.854727e+01 7.144859e+00 1.477683e+01
#>  [846] 6.005342e+01 9.554914e+01 5.975301e+01 9.051631e+01 3.558096e+01
#>  [851] 6.911040e+01 9.054063e+00 7.000988e+00 7.004185e+00 7.007395e+00
#>  [856] 7.329337e+00 2.990138e+01 7.335133e+00 7.338049e+00 7.340977e+00
#>  [861] 7.343917e+00 7.346868e+00 7.349832e+00 1.097510e+01 7.355794e+00
#>  [866] 7.358793e+00 7.361803e+00 7.707468e+00 7.710108e+00 7.712759e+00
#>  [871] 7.715420e+00 7.718092e+00 7.720775e+00 7.723468e+00 7.726171e+00
#>  [876] 7.728885e+00 7.731608e+00 7.734342e+00 7.737086e+00 8.083891e+00
#>  [881] 8.086235e+00 2.244799e+01 1.323094e+01 8.093321e+00 8.095701e+00
#>  [886] 8.098090e+00 3.379739e+01 1.040200e+01 8.105307e+00 8.107729e+00
#>  [891] 8.110160e+00 8.210985e+00 2.078261e+01 8.215634e+00 8.217970e+00
#>  [896] 7.756278e+01 6.945300e+01 8.225024e+00 6.400504e+01 3.536798e+01
#>  [901] 8.232142e+00 2.983906e+01 8.236923e+00 8.553172e+00 1.579591e+01
#>  [906] 7.816832e+01 7.294101e+01 3.246099e+01 3.491781e+01 7.099668e+01
#>  [911] 7.680983e+01 7.572611e+01 5.508983e+01 7.129444e+01 4.564261e+01
#>  [916] 3.245805e+01 7.689925e+01 6.101065e+01 6.323207e+01 7.798371e+01
#>  [921] 6.831078e+01 7.068451e+01 3.285097e+01 7.259775e+01 7.414448e+01
#>  [926] 7.374764e+01 6.409439e+01 7.113979e+01 6.592926e+01 6.822532e+01
#>  [931] 7.329244e+01 7.086435e+01 7.460217e+01 6.176398e+01 7.130831e+01
#>  [936] 6.561027e+01 7.006540e+01 2.749331e+01 6.731304e+01 6.901371e+01
#>  [941] 6.448552e+01 6.047461e+01 6.814012e+01 6.527428e+01 6.841529e+01
#>  [946] 9.314032e+00 5.177607e+01 2.394511e+01 6.166011e+01 5.130209e+01
#>  [951] 6.432075e+01 7.019272e+01 6.979145e+01 6.215676e+01 6.318134e+01
#>  [956] 6.432799e+01 6.542974e+01 6.533403e+01 5.951464e+01 5.840353e+01
#>  [961] 6.149246e+01 6.042599e+01 6.021323e+01 3.108368e+01 2.912938e+01
#>  [966] 1.743852e+01 6.488967e+01 5.953261e+01 6.258716e+01 2.605246e+01
#>  [971] 2.301616e+01 6.271290e+01 1.827708e+01 5.827544e+01 5.644321e+01
#>  [976] 1.088300e+01 1.088104e+01 6.223437e+01 6.693904e+01 1.087513e+01
#>  [981] 2.592845e+01 6.522711e+01 6.789938e+01 6.571771e+01 6.362275e+01
#>  [986] 6.297338e+01 6.607908e+01 6.264127e+01 6.340301e+01 6.038544e+01
#>  [991] 6.282939e+01 6.191179e+01 6.294635e+01 6.250037e+01 6.184439e+01
#>  [996] 6.180800e+01 6.028564e+01 6.207714e+01 5.991156e+01 2.593312e+01
#> [1001] 6.433452e+01 6.476920e+01 6.452322e+01 6.249592e+01 6.355227e+01
#> [1006] 6.340826e+01 6.290247e+01 6.203932e+01 3.491161e+01 6.174058e+01
#> [1011] 5.663150e+01 6.051251e+01 6.107185e+01 5.922795e+01 5.932485e+01
#> [1016] 6.031649e+01 5.866926e+01 5.865420e+01 5.674684e+01 5.911336e+01
#> [1021] 5.855976e+01 5.786862e+01 5.406105e+01 5.871001e+01 5.768386e+01
#> [1026] 5.497184e+01 5.607644e+01 5.567198e+01 2.905979e+01 5.594020e+01
#> [1031] 2.301741e+01 5.528356e+01 1.240462e+01 1.519675e+01 1.239090e+01
#> [1036] 1.284946e+01 1.284136e+01 1.283323e+01 5.669409e+01 5.418090e+01
#> [1041] 1.280870e+01 1.280047e+01 1.279222e+01 1.278394e+01 1.277563e+01
#> [1046] 1.557234e+01 1.275894e+01 4.025367e+00 1.717103e+00 5.736071e+00
#> [1051] 2.270824e+00 1.732099e+00 1.737127e+00 1.742171e+00 1.747230e+00
#> [1056] 1.752304e+00 1.757393e+00 1.762497e+00 1.814802e+00 1.820007e+00
#> [1061] 1.825228e+00 1.830464e+00 1.835716e+00 1.840983e+00 1.846266e+00
#> [1066] 1.851564e+00 1.856878e+00 1.862207e+00 1.867552e+00 1.872913e+00
#> [1071] 1.918775e+00 1.924224e+00 1.929689e+00 6.619674e+00 1.940668e+00
#> [1076] 1.946180e+00 1.951709e+00 1.957253e+00 1.962814e+00 1.968390e+00
#> [1081] 1.973982e+00 1.979590e+00 2.023116e+00 2.028808e+00 2.034515e+00
#> [1086] 2.040239e+00 2.045979e+00 2.051736e+00 2.057508e+00 2.063296e+00
#> [1091] 2.069101e+00 2.074922e+00 2.080758e+00 2.086611e+00 2.093067e+00
#> [1096] 2.098954e+00 2.104856e+00 2.110775e+00 2.116710e+00 2.122661e+00
#> [1101] 2.128628e+00 2.134612e+00 2.140612e+00 2.146628e+00 2.152660e+00
#> [1106] 2.158709e+00 1.912364e+00 3.847917e+00 1.922855e+00 1.928124e+00
#> [1111] 1.933408e+00 1.938709e+00 8.120154e+00 4.549857e+00 1.954704e+00
#> [1116] 1.960067e+00 1.965446e+00 2.019396e+00 2.024863e+00 2.030345e+00
#> [1121] 2.035844e+00 2.041358e+00 2.046888e+00 2.052435e+00 2.057997e+00
#> [1126] 2.063576e+00 2.069170e+00 2.074781e+00 2.080408e+00 2.133367e+00
#> [1131] 2.139077e+00 2.144804e+00 2.150547e+00 2.156306e+00 2.162081e+00
#> [1136] 2.167873e+00 2.173681e+00 2.179506e+00 2.185346e+00 2.191204e+00
#> [1141] 2.197077e+00 9.002883e+01 4.190811e+01 4.864343e+01 3.105767e+01
#> [1146] 1.496281e+01 1.230772e+02 2.844309e+01 2.635440e+01 5.446920e+01
#> [1151] 7.144062e+01 7.902907e+01 3.169925e+01 7.579065e+01 2.002602e+01
#> [1156] 7.037259e+01 7.438948e+01 1.000951e+01 1.358075e+01 6.649714e+00
#> [1161] 4.930853e+01 5.354434e+01 2.683916e+01 3.255063e+01 6.892048e+01
#> [1166] 2.402264e+01 1.145860e+01 2.793994e+01 6.037053e+01 9.546292e+00
#> [1171] 2.157223e+01 1.316307e+02 1.792358e+02 4.864115e+01 2.947744e+01
#> [1176] 2.226507e+01 2.181341e+01 1.682727e+01 1.580046e+01 6.560530e+00
#> [1181] 1.083774e+01 7.040495e+01 7.162979e+00 7.345930e+00 6.093294e+00
#> [1186] 5.892877e+00 5.789739e+00 5.973506e+00 9.923102e+00 5.704424e+00
#> [1191] 5.957486e+00 6.296806e+00 5.869820e+00 7.122280e+00 8.521461e+00
#> [1196] 6.734163e+00 9.532865e+00 8.006599e+00 6.914629e+00 6.908659e+00
#> [1201] 6.349939e+00 5.679468e+00 6.031634e+00 4.974164e+00 4.977634e+00
#> [1206] 4.137464e+00 3.087623e+00 4.712900e+00 6.812291e+00 2.832173e+00
#> [1211] 5.527126e+00 2.352083e+00 1.895034e+01 2.275808e+01 2.481113e+01
#> [1216] 2.431293e+01 2.651225e+01 3.458636e+01 4.178883e+01 1.877209e+01
#> [1221] 5.983227e+01 5.657916e+01 5.361188e+01 5.579076e+01 5.763511e+01
#> [1226] 4.966060e+01 4.581349e+01 2.241999e+01 1.618725e+01 4.278142e+01
#> [1231] 4.765905e+01 5.329569e+01 5.230866e+01 3.951270e+01 3.942395e+01
#> [1236] 4.789629e+01 4.892502e+01 5.098724e+01 4.838828e+01 3.909949e+01
#> [1241] 5.251379e+01 3.892856e+01 5.073696e+01 5.192575e+01 4.832302e+01
#> [1246] 1.163962e+01 5.028629e+01 5.084534e+01 5.158884e+01 4.307299e+01
#> [1251] 5.029379e+01 5.233188e+01 4.881485e+01 4.027663e+01 5.643148e+01
#> [1256] 5.459251e+01 4.541163e+01 3.986901e+01 6.620243e+00 6.625643e+00
#> [1261] 6.631044e+00 6.636448e+00 6.641853e+00 6.647259e+00 6.652667e+00
#> [1266] 6.857238e+00 6.862357e+00 6.867476e+00 3.924739e+01 3.914526e+01
#> [1271] 3.904336e+01 3.894170e+01 3.884027e+01 3.873908e+01 6.903298e+00
#> [1276] 3.853741e+01 6.913526e+00 3.880221e+01 7.124207e+00 7.128996e+00
#> [1281] 7.133783e+00 7.138567e+00 7.143348e+00 7.148126e+00 7.152900e+00
#> [1286] 7.157671e+00 7.162437e+00 7.167200e+00 7.171959e+00 1.135324e+00
#> [1291] 1.137020e+00 1.138721e+00 1.140427e+00 1.142139e+00 1.143855e+00
#> [1296] 1.145577e+00 1.147304e+00 1.149036e+00 1.150773e+00 1.152515e+00
#> [1301] 1.190045e+00 1.191783e+00 1.193526e+00 1.195275e+00 1.197028e+00
#> [1306] 1.198787e+00 1.200551e+00 1.202320e+00 1.204094e+00 1.205874e+00
#> [1311] 1.207658e+00 1.209447e+00 1.249354e+00 3.578138e+00 1.252923e+00
#> [1316] 1.254715e+00 1.256513e+00 1.258315e+00 1.260122e+00 1.908959e+00
#> [1321] 5.017014e+00 3.305749e+00 1.267401e+00 1.269233e+00 3.947604e+00
#> [1326] 1.307480e+00 1.309309e+00 1.311143e+00 1.312982e+00 1.314826e+00
#> [1331] 1.316674e+00 1.318528e+00 1.320386e+00 1.322249e+00 1.324117e+00
#> [1336] 1.325989e+00 1.317368e+00 1.319256e+00 1.321149e+00 1.323047e+00
#> [1341] 2.211619e+00 1.326856e+00 1.328768e+00 1.330685e+00 1.332606e+00
#> [1346] 1.334532e+00 1.336463e+00 1.338398e+00 1.370685e+00 1.372612e+00
#> [1351] 1.374543e+00 1.376479e+00 1.378420e+00 1.380365e+00 1.382315e+00
#> [1356] 1.384269e+00 1.386228e+00 1.388191e+00 1.390158e+00 1.392130e+00
#> [1361] 1.425117e+00 1.427078e+00 1.429043e+00 1.431013e+00 1.432986e+00
#> [1366] 1.434964e+00 1.436946e+00 1.438933e+00 1.440923e+00 1.442918e+00
#> [1371] 1.444917e+00 1.446920e+00 1.508748e+00 1.510716e+00 1.512687e+00
#> [1376] 1.514663e+00 1.516642e+00 1.518625e+00 4.443495e+00 1.522603e+00
#> [1381] 1.524598e+00 1.526596e+00 1.528598e+00 1.530604e+00 1.572289e+00
#> [1386] 1.574269e+00 1.268619e+01 1.578241e+00 1.580232e+00 1.260273e+01
#> [1391] 1.584225e+00 1.586227e+00 1.251962e+01 1.590240e+00 1.246440e+01
#> [1396] 1.243684e+01 1.254874e+01 1.639414e+00 1.641401e+00 1.643392e+00
#> [1401] 1.645385e+00 1.647381e+00 1.649380e+00 1.651382e+00 1.653387e+00
#> [1406] 1.229468e+01 1.226667e+01 1.659419e+00 1.234337e+01 1.231460e+01
#> [1411] 1.228588e+01 1.225721e+01 1.711304e+00 1.713291e+00 1.217150e+01
#> [1416] 1.214303e+01 1.211461e+01 1.208625e+01 1.205793e+01 1.725262e+00
#> [1421] 1.770107e+00 1.772065e+00 1.774025e+00 1.775987e+00 4.795981e+00
#> [1426] 1.779917e+00 1.781885e+00 1.783854e+00 1.785825e+00 1.787798e+00
#> [1431] 1.789772e+00 2.860152e+00 1.837528e+00 1.187470e+01 1.184557e+01
#> [1436] 1.181649e+01 1.845242e+00 1.175852e+01 1.172962e+01 1.170079e+01
#> [1441] 1.852975e+00 1.854910e+00 1.856847e+00 1.858785e+00 1.904438e+00
#> [1446] 1.478320e+01 1.161179e+01 1.158266e+01 1.447588e+01 1.152460e+01
#> [1451] 1.465085e+01 1.486330e+01 1.143800e+01 1.140927e+01 3.239142e+00
#> [1456] 1.135199e+01 9.123365e+01 7.602488e+01 8.251351e+01 5.328530e+01
#> [1461] 6.877448e+01 5.305228e+01 5.631536e+00 7.712782e+01 5.748312e+01
#> [1466] 5.258739e+01 7.417779e+01 6.827860e+01 6.461559e+01 6.533059e+01
#> [1471] 5.280382e+01 5.843503e+00 5.256453e+01 5.244507e+01 5.232574e+01
#> [1476] 5.220653e+01 1.818656e+01 2.209425e+01 5.875059e+00 6.253145e+00
#> [1481] 6.740603e+00 1.384156e+01 1.209039e+01 1.225003e+01 1.306222e+01
#> [1486] 1.315138e+01 1.448243e+01 1.487905e+01 1.630443e+01 1.853850e+01
#> [1491] 1.519228e+01 1.836949e+01 1.931994e+01 1.087596e+01 8.843954e+00
#> [1496] 8.594499e+00 4.451671e+00 2.916467e+00 1.623722e+00 1.660016e+00
#> [1501] 1.173446e+00 1.074355e+01 6.478390e-01 2.309117e+01 3.326376e+00
#> [1506] 3.332137e+00 1.100668e+01 3.343710e+00 3.349521e+00 3.355350e+00
#> [1511] 3.361195e+00 3.367058e+00 3.372937e+00 3.378833e+00 3.429862e+00
#> [1516] 4.108669e+01 3.441720e+00 3.447674e+00 3.453645e+00 3.459633e+00
#> [1521] 4.903426e+00 4.606187e+00 4.475216e+00 4.476785e+00 4.614702e+00
#> [1526] 4.603025e+00 4.498133e+00 3.499147e+00 3.503077e+00 1.503312e+01
#> [1531] 1.764559e+01 2.048591e+01 2.546665e+01 3.086794e+01 3.730644e+01
#> [1536] 3.783491e+01 3.135759e+01 2.888438e+01 2.454756e+01 1.787367e+01
#> [1541] 1.824961e+01 1.524427e+01 1.555975e+01 1.466923e+01 1.582094e+01
#> [1546] 1.551457e+01 1.203454e+01 1.204798e+01 7.416939e+00 6.330910e+00
#> [1551] 3.877606e+00 3.884091e+00 3.890591e+00 3.897108e+00 3.903641e+00
#> [1556] 3.910189e+00 5.885284e+00 3.923332e+00 3.929928e+00 3.936539e+00
#> [1561] 4.060415e+00 4.067009e+00 4.073619e+00 4.080243e+00 4.086883e+00
#> [1566] 4.093538e+00 4.100208e+00 5.309734e+00 4.113594e+00 4.120309e+00
#> [1571] 4.127039e+00 4.133784e+00 6.396341e+00 4.288754e+00 4.295474e+00
#> [1576] 4.302207e+00 4.308955e+00 4.315717e+00 4.322493e+00 4.329284e+00
#> [1581] 4.336088e+00 4.342907e+00 4.349739e+00 4.356585e+00 4.018796e+00
#> [1586] 4.025832e+00 4.032882e+00 4.039946e+00 4.047025e+00 4.054118e+00
#> [1591] 4.061225e+00 4.068346e+00 4.075482e+00 4.082631e+00 4.089795e+00
#> [1596] 4.096972e+00 8.269974e+00 3.806526e+00 3.812140e+00 3.817771e+00
#> [1601] 9.580879e+00 5.811773e+00 1.231705e+01 7.513716e+00 3.846174e+00
#> [1606] 7.530460e+00 3.857653e+00 3.061257e+00 3.067220e+00 3.073200e+00
#> [1611] 3.079196e+00 6.736849e+00 8.621036e+00 3.097288e+00 3.103352e+00
#> [1616] 3.109433e+00 3.115531e+00 3.121646e+00 3.127778e+00 3.196060e+00
#> [1621] 3.202224e+00 3.208405e+00 3.214603e+00 3.220818e+00 3.227049e+00
#> [1626] 3.233298e+00 3.239563e+00 3.245845e+00 3.252144e+00 3.258459e+00
#> [1631] 3.264791e+00 3.333508e+00 3.339870e+00 3.346248e+00 3.352642e+00
#> [1636] 3.359054e+00 3.365481e+00 3.371926e+00 3.378387e+00 3.384865e+00
#> [1641] 3.391359e+00 3.397870e+00 3.404397e+00 3.249186e+00 3.254942e+00
#> [1646] 3.260715e+00 3.266506e+00 3.272313e+00 3.278137e+00 3.283978e+00
#> [1651] 3.289836e+00 3.295711e+00 1.377387e+01 1.139560e+01 1.769267e+01
#> [1656] 5.736294e+01 7.030507e+00 1.421874e+01 5.815581e+01 5.198998e+01
#> [1661] 5.696738e+01 9.036886e+00 1.181635e+01 1.056052e+01 5.668229e+01
#> [1666] 9.296034e+00 2.924940e+01 2.926313e+01 1.471414e+00 1.475968e+00
#> [1671] 2.930303e+01 2.931588e+01 2.932852e+01 1.494323e+00 1.498947e+00
#> [1676] 1.503585e+00 1.508237e+00 1.530543e+00 1.535259e+00 1.539990e+00
#> [1681] 1.544735e+00 1.549494e+00 1.554268e+00 1.559056e+00 1.563858e+00
#> [1686] 1.568675e+00 2.961467e+01 2.962329e+01 2.963168e+01 1.606061e+00
#> [1691] 2.977565e+01 2.978270e+01 2.978953e+01 1.625851e+00 2.980250e+01
#> [1696] 2.980865e+01 1.640848e+00 1.645877e+00 2.167066e+00 4.779309e+00
#> [1701] 8.130482e+01 1.694299e+02 2.170973e+02 2.310037e+02 2.705099e+01
#> [1706] 2.601612e+01 2.244598e+02 2.213869e+02 2.183560e+02 1.791721e+02
#> [1711] 1.802946e+02 5.836135e+01 3.410392e+00 5.036203e+01 3.402602e+00
#> [1716] 4.966547e+01 3.515249e+00 5.063113e+01 4.574618e+01 3.101878e+00
#> [1721] 3.112317e+00 2.682294e+00 2.813083e+00 3.226722e+00 3.323957e+00
#> [1726] 3.364520e+00 3.861209e+00 3.873018e+00 3.955369e+00 2.347378e+00
#> [1731] 2.083016e+00 2.071160e+00 2.042578e+00 2.048753e+00 3.202174e+01
#> [1736] 2.191232e+00 3.575142e+00 4.843459e+01 3.628000e+00 4.052205e+00
#> [1741] 4.101109e+00 4.150597e+00 4.163132e+00 4.175695e+00 4.188288e+00
#> [1746] 2.145121e+00 1.989713e+00 1.995657e+00 2.001617e+00 1.641413e+00
#> [1751] 1.646269e+00 1.651139e+00 1.656025e+00 1.660925e+00 1.665839e+00
#> [1756] 7.851523e+00 2.582438e+00 1.680673e+00 1.685647e+00 1.690637e+00
#> [1761] 1.695924e+00 1.700943e+00 2.242471e+00 1.711028e+00 5.512106e+00
#> [1766] 6.775559e+00 1.726268e+00 1.731378e+00 1.736503e+00 1.741643e+00
#> [1771] 1.746799e+00 1.751970e+00 4.507394e+00 1.733908e+00 6.752502e+00
#> [1776] 2.682520e+00 3.721270e+00 1.341465e+01 1.170419e+01 9.870958e+00
#> [1781] 1.770404e+00 5.354519e+00 1.780969e+00 1.786274e+00 1.815567e+00
#> [1786] 2.795828e+00 6.261888e+00 4.885749e+00 4.091963e+00 1.842610e+00
#> [1791] 5.428287e+00 4.641626e+00 5.455035e+00 2.444937e+00 9.512051e-01
#> [1796] 9.527678e-01 3.202283e-01 1.735948e-01 1.272051e-01 1.151344e-01
#> [1801] 8.344166e-02 4.421497e-02 1.862606e-02 1.573445e-02 1.487034e-02
#> [1806] 1.500242e-02 1.447765e-02 1.414266e-02 1.445428e-02 1.524760e-02
#> [1811] 1.626819e-02 1.575458e-02 1.702134e-02 1.700821e-02 1.692459e-02
#> [1816] 1.550567e-02 1.569644e-02 1.547077e-02 1.469785e-02 1.515167e-02
#> [1821] 1.773735e+00 1.778794e+00 1.783868e+00 1.788958e+00 1.794062e+00
#> [1826] 1.799182e+00 1.804317e+00 1.809468e+00 1.814634e+00 1.819815e+00
#> [1831] 1.825012e+00 1.861326e+00 1.866597e+00 1.871883e+00 1.877185e+00
#> [1836] 1.882503e+00 1.887836e+00 1.893186e+00 1.898550e+00 1.903931e+00
#> [1841] 1.909327e+00 1.914738e+00 1.920166e+00 1.957767e+00 1.963270e+00
#> [1846] 1.968788e+00 1.974322e+00 1.979873e+00 1.985439e+00 1.991021e+00
#> [1851] 1.996619e+00 2.002233e+00 2.007863e+00 2.013509e+00 2.019171e+00
#> [1856] 1.368354e+01 1.366338e+01 1.364320e+01 1.362299e+01 1.360277e+01
#> [1861] 1.112325e+00 1.114058e+00 1.115796e+00 1.117539e+00 1.119287e+00
#> [1866] 1.121041e+00 1.157917e+00 1.159668e+00 1.161425e+00 1.163186e+00
#> [1871] 1.164953e+00 1.166725e+00 1.168502e+00 1.170284e+00 1.172071e+00
#> [1876] 1.173864e+00 1.175661e+00 1.177464e+00 1.141372e+00 1.143199e+00
#> [1881] 1.145031e+00 1.146868e+00 1.148710e+00 1.150557e+00 1.152409e+00
#> [1886] 1.154267e+00 1.156129e+00 1.301312e+01 1.159869e+00 1.161747e+00
#> [1891] 1.209771e+00 1.211435e+00 1.213105e+00 1.214780e+00 1.216460e+00
#> [1896] 1.218145e+00 1.219836e+00 1.221531e+00 1.223232e+00 1.224938e+00
#> [1901] 1.226649e+00 1.263371e+00 1.265075e+00 1.266784e+00 1.268498e+00
#> [1906] 1.270218e+00 1.271942e+00 1.273671e+00 1.275406e+00 1.277145e+00
#> [1911] 1.278889e+00 1.280639e+00 1.282393e+00 1.273519e+00 1.275289e+00
#> [1916] 1.277064e+00 1.278845e+00 1.280630e+00 1.282420e+00 1.284215e+00
#> [1921] 1.286015e+00 1.287820e+00 1.289630e+00 1.291445e+00 1.293265e+00
#> [1926] 3.824356e-01 3.836697e-01 3.849078e-01 8.449409e+00 9.896995e+00
#> [1931] 9.900929e+00 8.469134e+00 3.911570e-01 3.924186e-01 1.505915e+01
#> [1936] 1.550677e+01 1.558454e+01 4.033079e-01 4.046034e-01 4.059029e-01
#> [1941] 4.072065e-01 4.085141e-01 4.098259e-01 4.111416e-01 4.124615e-01
#> [1946] 4.137855e-01 4.151135e-01 4.164457e-01 4.237828e-01 4.251374e-01
#> [1951] 4.264962e-01 4.278592e-01 4.292263e-01 4.305977e-01 4.319732e-01
#> [1956] 4.333530e-01 4.347369e-01 4.361251e-01 4.375175e-01 4.389142e-01
#> [1961] 4.141952e-01 4.155027e-01 4.168143e-01 4.181299e-01 1.076507e+00
#> [1966] 4.207736e-01 4.221016e-01 8.004906e-01 1.193904e+00 4.261106e-01
#> [1971] 4.274552e-01 7.680458e-01 7.252199e-01 3.564701e-01 3.576406e-01
#> [1976] 1.251336e+00 3.599929e-01 3.611747e-01 1.054934e+00 4.798949e-01
#> [1981] 3.647425e-01 3.659394e-01 3.727731e-01 3.739921e-01 3.752150e-01
#> [1986] 3.764417e-01 3.776722e-01 3.789067e-01 3.801450e-01 3.813872e-01
#> [1991] 8.322346e+00 8.328854e+00 3.851373e-01 3.863951e-01 3.935592e-01
#> [1996] 3.948399e-01 3.961247e-01 3.974135e-01 3.987063e-01 4.000031e-01
#> [2001] 4.013039e-01 4.026087e-01 4.039176e-01 4.052305e-01 4.065475e-01
#> [2006] 4.078685e-01 4.153697e-01 4.167144e-01 4.180632e-01 4.194162e-01
#> [2011] 4.207733e-01 4.221346e-01 4.235000e-01 4.248696e-01 4.262433e-01
#> [2016] 4.276213e-01 4.290034e-01 4.303897e-01 6.561079e+00 2.650734e+00
#> [2021] 2.656483e+00 2.662248e+00 4.546162e+00 2.673828e+00 2.679643e+00
#> [2026] 1.362480e+01 2.691324e+00 2.697190e+00 1.124672e+01 2.566388e+00
#> [2031] 3.626602e+01 3.623941e+01 4.455183e+01 9.765792e+00 2.595894e+00
#> [2036] 4.426400e+00 3.610363e+01 3.607594e+01 6.006527e+00 2.625819e+00
#> [2041] 2.631854e+00 2.482254e+00 3.508660e+01 5.182556e+01 3.503912e+01
#> [2046] 2.506319e+00 2.512377e+00 2.573680e+00 2.579952e+00 3.166791e+00
#> [2051] 4.782232e+00 4.492336e+00 7.649333e+01 1.443098e+01 2.015573e+01
#> [2056] 2.971223e+01 3.727492e+01 3.126052e+02 2.847909e+02 2.879886e+02
#> [2061] 1.210702e+02 1.688529e+02 1.189363e+01 9.141055e+01 6.904293e+01
#> [2066] 6.604834e+01 5.573879e+00 5.588777e+00 5.869803e+00 8.067765e+01
#> [2071] 6.455936e+00 6.171314e+00 1.016259e+01 1.480536e+01 1.136415e+02
#> [2076] 1.392496e+01 9.545082e+01 1.047401e+01 6.434363e+00 6.337386e+00
#> [2081] 2.690577e+00 3.367759e+01 3.364963e+01 3.362150e+01 3.359321e+01
#> [2086] 2.676866e+00 2.683539e+00 3.350734e+01 2.696936e+00 3.402912e+01
#> [2091] 3.399572e+01 3.396217e+01 3.392847e+01 4.132876e+01 5.018741e+01
#> [2096] 3.711824e+01 5.267258e+01 4.426182e+01 2.913121e+00 5.141526e+01
#> [2101] 3.868033e+01 5.207602e+01 4.981451e+01 5.516246e+01 4.977050e+01
#> [2106] 5.186835e+01 3.420721e+01 3.416575e+01 5.097549e+01 5.137756e+01
#> [2111] 5.073577e+01 4.370385e+01 1.882588e+01 3.442953e+00 3.450282e+00
#> [2116] 3.457626e+00 3.464985e+00 3.472359e+00 3.479748e+00 3.487152e+00
#> [2121] 3.494571e+00 3.502004e+00 3.509453e+00 3.516916e+00 3.524395e+00
#> [2126] 3.446803e+01 1.808005e+01 3.706987e+00 3.714518e+00 3.722063e+00
#> [2131] 3.729621e+00 3.737194e+00 3.744781e+00 3.752381e+00 3.759995e+00
#> [2136] 3.767623e+00 3.775265e+00 3.785184e+00 3.792852e+00 3.800533e+00
#> [2141] 3.808228e+00 3.815936e+00 1.492203e+01 2.000585e+01 1.000486e+01
#> [2146] 3.846898e+00 3.854670e+00 3.862456e+00 3.870254e+00 4.138816e+00
#> [2151] 4.146569e+00 4.154333e+00 1.507218e+01 4.169896e+00 4.177694e+00
#> [2156] 4.185503e+00 4.193323e+00 4.201154e+00 4.208997e+00 4.216849e+00
#> [2161] 4.224713e+00 4.280639e+00 4.288502e+00 4.296376e+00 4.304260e+00
#> [2166] 4.312154e+00 4.320057e+00 4.327971e+00 4.335894e+00 4.343826e+00
#> [2171] 4.351768e+00 4.359719e+00 4.367679e+00 4.419988e+00 4.427943e+00
#> [2176] 4.435907e+00 4.443879e+00 4.451860e+00 4.459849e+00 4.467846e+00
#> [2181] 4.475851e+00 4.483864e+00 4.491884e+00 4.499913e+00 4.507948e+00
#> [2186] 3.290010e-01 3.300967e-01 3.311959e-01 3.322987e-01 3.334050e-01
#> [2191] 3.345148e-01 3.356282e-01 3.367452e-01 3.378658e-01 3.389899e-01
#> [2196] 3.401176e-01 3.453094e-01 3.464555e-01 3.476053e-01 3.487587e-01
#> [2201] 3.499158e-01 3.510766e-01 3.522411e-01 3.534093e-01 3.545812e-01
#> [2206] 1.699504e+01 8.051973e+00 8.059370e+00 1.036458e+01 1.216488e+01
#> [2211] 8.124430e+00 1.685658e+01 1.747638e+01 8.145282e+00 1.222701e+00
#> [2216] 4.907755e-01 3.732560e-01 3.744886e-01 3.757250e-01 4.971849e-01
#> [2221] 3.826281e-01 3.838881e-01 3.851519e-01 3.864197e-01 3.876914e-01
#> [2226] 3.889671e-01 3.902467e-01 3.915303e-01 3.928178e-01 3.941094e-01
#> [2231] 3.954049e-01 3.967044e-01 4.026429e-01 6.991775e-01 8.354366e+00
#> [2236] 1.183523e+01 1.298928e+01 4.092826e-01 1.182303e+01 1.503382e+01
#> [2241] 1.181447e+01 8.906456e-01 2.833773e+00 1.955116e+00 4.235843e-01
#> [2246] 4.249659e-01 4.263517e-01 4.277417e-01 4.291359e-01 4.305342e-01
#> [2251] 4.319368e-01 4.333436e-01 4.347546e-01 4.361698e-01 4.375893e-01
#> [2256] 4.390129e-01 4.454783e-01 6.888657e-01 4.483730e-01 1.710379e+00
#> [2261] 4.512849e-01 4.527474e-01 4.542141e-01 4.556852e-01 2.128424e+00
#> [2266] 1.539363e+01 2.751866e+00 1.409629e+01 8.587238e+00 9.805296e+00
#> [2271] 8.591922e+00 8.594161e+00 1.282128e+01 1.501737e+01 1.202005e+00
#> [2276] 1.049523e+01 2.101126e+00 4.813309e-01 1.513998e+01 1.383695e+01
#> [2281] 4.913836e-01 4.929589e-01 4.945387e-01 4.961230e-01 4.977118e-01
#> [2286] 4.993051e-01 5.009029e-01 5.025053e-01 5.041122e-01 5.057236e-01
#> [2291] 5.073396e-01 5.089601e-01 5.161675e-01 5.178109e-01 9.671445e-01
#> [2296] 7.985408e-01 1.728990e+00 8.033899e-01 5.260964e-01 6.914729e-01
#> [2301] 5.294429e-01 5.311231e-01 5.328078e-01 5.344972e-01 7.078334e-01
#> [2306] 5.422887e-01 1.009494e+00 5.457211e-01 5.474443e-01 1.486633e+00
#> [2311] 9.390732e-01 5.526419e-01 5.543838e-01 5.561304e-01 1.033145e+00
#> [2316] 5.596376e-01 8.659506e-01 5.690849e-01 8.706218e+00 8.704501e+00
#> [2321] 8.702719e+00 8.700872e+00 5.780517e-01 5.798592e-01 5.816715e-01
#> [2326] 5.834884e-01 5.853101e-01 5.871366e-01 5.951142e-01 2.245134e+00
#> [2331] 5.988189e-01 6.006784e-01 6.025426e-01 7.887938e-01 6.062854e-01
#> [2336] 6.081639e-01 6.100471e-01 6.119351e-01 8.006461e-01 6.157254e-01
#> [2341] 6.274692e-01 6.293981e-01 6.313317e-01 6.332701e-01 6.352133e-01
#> [2346] 6.371613e-01 6.391140e-01 6.410714e-01 6.430337e-01 6.450007e-01
#> [2351] 6.469724e-01 6.489489e-01 6.580040e-01 6.600051e-01 6.620110e-01
#> [2356] 6.640217e-01 6.660371e-01 6.680572e-01 6.700821e-01 6.721118e-01
#> [2361] 6.741461e-01 6.761852e-01 6.782291e-01 6.802776e-01 7.416035e-01
#> [2366] 7.441768e-01 7.467587e-01 7.493492e-01 7.519483e-01 7.545562e-01
#> [2371] 7.571727e-01 7.597980e-01 7.624320e-01 7.650748e-01 2.950889e+00
#> [2376] 2.430549e+00 2.599452e+00 2.388746e+00 7.884938e-01 7.912223e-01
#> [2381] 7.939599e-01 7.967066e-01 7.994624e-01 1.403778e+00 1.857823e+00
#> [2386] 8.077846e-01 8.105770e-01 8.229744e-01 8.258161e-01 8.286673e-01
#> [2391] 8.315278e-01 8.343978e-01 8.372772e-01 8.401661e-01 1.313242e+00
#> [2396] 8.459724e-01 1.622032e+00 8.518169e-01 1.130933e+00 8.669052e-01
#> [2401] 8.698906e-01 1.358936e+00 1.158640e+00 8.789053e-01 8.819299e-01
#> [2406] 8.849644e-01 8.880087e-01 8.910630e-01 8.941272e-01 8.972014e-01
#> [2411] 9.002856e-01 9.134280e-01 9.165642e-01 9.197105e-01 9.228670e-01
#> [2416] 9.260338e-01 9.292108e-01 9.323981e-01 9.355956e-01 9.388036e-01
#> [2421] 9.420218e-01 9.452505e-01 9.484895e-01 9.620176e-01 9.653099e-01
#> [2426] 9.686128e-01 9.719263e-01 9.752505e-01 9.785853e-01 9.819307e-01
#> [2431] 9.852869e-01 9.886538e-01 9.920315e-01 9.954199e-01 9.988192e-01
#> [2436] 1.048564e+00 7.945831e-01 7.973275e-01 2.170479e+00 8.028436e-01
#> [2441] 8.056155e-01 1.547080e+00 8.111867e-01 8.139862e-01 8.167949e-01
#> [2446] 8.196129e-01 8.272975e-01 8.301497e-01 8.330113e-01 8.358823e-01
#> [2451] 8.387628e-01 8.416528e-01 8.445524e-01 8.474614e-01 8.503801e-01
#> [2456] 8.533083e-01 8.562462e-01 8.591937e-01 8.621509e-01 8.651178e-01
#> [2461] 8.680944e-01 8.710808e-01 8.740770e-01 8.770830e-01 8.800988e-01
#> [2466] 8.831245e-01 8.861601e-01 8.892056e-01 8.922610e-01 8.953264e-01
#> [2471] 1.886364e+00 1.891571e+00 5.102852e+01 5.494628e+01 6.103428e+01
#> [2476] 8.006616e+00 4.076549e+00 6.715305e+00 1.928454e+00 7.814450e+00
#> [2481] 1.939134e+00 4.358232e+01 5.076576e+00 6.641506e+00 6.358586e+00
#> [2486] 1.973429e+00 7.952259e+00 1.984347e+00 7.243513e+00 6.322649e+00
#> [2491] 2.000844e+00 4.006802e+00 2.011922e+00 6.901043e+00 6.186040e+00
#> [2496] 6.636520e+00 2.040832e+00 2.046467e+00 2.052118e+00 6.134942e+00
#> [2501] 4.968094e+00 5.323268e+00 9.407462e+00 2.080615e+00 5.015234e+00
#> [2506] 5.141094e+01 5.033890e+01 4.827669e+01 7.688437e+00 5.023909e+01
#> [2511] 5.239141e+01 5.634285e+00 2.806027e+00 4.989842e+01 7.611102e+00
#> [2516] 1.067852e+01 2.167714e+00 2.193590e+00 4.699995e+01 8.165258e+00
#> [2521] 9.776613e+00 5.104946e+01 4.375132e+01 1.224091e+01 1.405187e+01
#> [2526] 5.698011e+00 7.856540e+00 1.015038e+01 3.946488e+01 1.151428e+01
#> [2531] 7.530565e+00 3.801271e+01 3.229024e+01 2.322163e+00 3.225416e+01
#> [2536] 3.223583e+01 2.341145e+00 5.124509e+00 5.935782e+00 3.216057e+01
#> [2541] 2.366685e+00 2.412466e+00 2.418952e+00 6.427995e+00 2.431973e+00
#> [2546] 2.438509e+00 7.662627e+00 4.390689e+01 8.751960e+00 3.971770e+01
#> [2551] 6.182417e+00 3.731407e+01 5.132295e+00 3.838000e+00 9.523992e+00
#> [2556] 3.212784e+01 4.728561e+01 4.507930e+01 1.479406e+01 3.202849e+01
#> [2561] 6.286415e+00 3.491135e+00 3.589709e+00 2.249116e+01 1.775289e+01
#> [2566] 2.250520e+01 5.385172e+01 1.278977e+01 7.442043e+00 4.721871e+01
#> [2571] 4.039251e+00 3.309147e+00 4.051473e+00 4.895041e+00 7.727377e+00
#> [2576] 1.402472e+01 1.455551e+01 8.553852e+00 6.374593e+01 3.292374e+03
#> [2581] 3.865117e+00 3.132179e+00 1.038150e+01 3.058322e+00 2.632834e+00
#> [2586] 4.462259e+01 7.333071e+00 3.045107e+01 1.952050e+00 2.298595e+00
#> [2591] 2.337303e+00 2.516741e+00 2.038956e+00 2.105676e+00 2.080006e+00
#> [2596] 2.078678e+00 2.137721e+00 2.252325e+00 2.149843e+00 2.149954e+00
#> [2601] 2.237784e+00 2.350948e+00 2.381817e+00 2.207288e+00 2.037048e+00
#> [2606] 2.020493e+00 1.999550e+00 2.005427e+00 2.058570e+00 2.260796e+00
#> [2611] 2.378372e+00 2.553502e+00 2.487720e+00 3.660615e-01 3.672552e-01
#> [2616] 3.684527e-01 3.696540e-01 1.503113e+00 3.720680e-01 3.732808e-01
#> [2621] 3.744973e-01 3.757178e-01 3.769421e-01 3.781702e-01 3.987904e-01
#> [2626] 4.000741e-01 4.013619e-01 4.026537e-01 4.039495e-01 4.052493e-01
#> [2631] 4.065532e-01 4.078612e-01 4.091732e-01 4.104893e-01 4.118094e-01
#> [2636] 4.131337e-01 4.058708e-01 4.071823e-01 4.084979e-01 4.098176e-01
#> [2641] 4.111414e-01 4.124692e-01 4.138012e-01 4.151372e-01 4.164774e-01
#> [2646] 4.178217e-01 4.191701e-01 4.205226e-01 4.280832e-01 6.634917e-01
#> [2651] 4.308393e-01 6.676214e-01 4.336123e-01 4.350051e-01 4.364022e-01
#> [2656] 4.378035e-01 4.392090e-01 4.406188e-01 4.420329e-01 4.434512e-01
#> [2661] 4.515047e-01 4.529475e-01 4.543946e-01 4.558461e-01 4.573018e-01
#> [2666] 4.587620e-01 4.602265e-01 4.616953e-01 4.631686e-01 4.646462e-01
#> [2671] 4.661282e-01 4.676145e-01 4.762082e-01 4.777201e-01 4.792364e-01
#> [2676] 4.807572e-01 4.822825e-01 7.445687e-01 4.853465e-01 4.868852e-01
#> [2681] 4.884284e-01 4.899761e-01 4.915283e-01 4.930850e-01 5.021504e-01
#> [2686] 5.037334e-01 9.447260e-01 7.788127e-01 5.085099e-01 5.101112e-01
#> [2691] 5.117171e-01 5.133276e-01 5.149427e-01 5.165624e-01 5.181867e-01
#> [2696] 5.198156e-01 5.317956e-01 5.334570e-01 5.351230e-01 1.878645e+00
#> [2701] 5.384691e-01 5.401492e-01 5.418339e-01 5.435234e-01 5.452175e-01
#> [2706] 5.469164e-01 5.486200e-01 5.503282e-01 5.597893e-01 5.615239e-01
#> [2711] 5.632632e-01 5.650074e-01 5.667562e-01 5.685099e-01 5.702683e-01
#> [2716] 5.720314e-01 5.737994e-01 5.755721e-01 5.773496e-01 5.791319e-01
#> [2721] 5.890158e-01 5.908248e-01 5.926386e-01 5.944573e-01 5.962807e-01
#> [2726] 5.981090e-01 5.999421e-01 6.017800e-01 6.036228e-01 6.054704e-01
#> [2731] 6.073228e-01 6.091801e-01 4.582514e-01 4.596490e-01 4.610510e-01
#> [2736] 4.624573e-01 4.638679e-01 4.652829e-01 4.667023e-01 4.681260e-01
#> [2741] 4.695540e-01 4.709865e-01 4.724233e-01 4.814936e-01 4.829536e-01
#> [2746] 4.844181e-01 9.255341e+00 9.258033e+00 9.260656e+00 9.263209e+00
#> [2751] 9.265693e+00 9.268108e+00 4.947943e-01 9.272728e+00 4.977993e-01
#> [2756] 9.334933e+00 1.666955e+00 1.520041e+00 5.117469e-01 5.132889e-01
#> [2761] 5.148355e-01 5.163867e-01 5.179425e-01 5.195029e-01 5.210680e-01
#> [2766] 5.226376e-01 5.242119e-01 5.336790e-01 5.352767e-01 5.368792e-01
#> [2771] 5.384863e-01 5.400981e-01 5.417146e-01 5.433359e-01 5.449618e-01
#> [2776] 5.465925e-01 5.482278e-01 5.498679e-01 5.515128e-01 5.609794e-01
#> [2781] 5.626474e-01 5.643203e-01 5.659979e-01 5.676803e-01 5.693676e-01
#> [2786] 5.710596e-01 5.727564e-01 5.744580e-01 5.761645e-01 5.778757e-01
#> [2791] 5.795918e-01 5.891689e-01 5.909081e-01 5.926521e-01 5.944011e-01
#> [2796] 5.961549e-01 5.979136e-01 5.996771e-01 6.014456e-01 6.032189e-01
#> [2801] 6.049971e-01 6.067802e-01 6.085683e-01 1.171263e+00 1.172945e+00
#> [2806] 1.174631e+00 1.176323e+00 1.178020e+00 1.179723e+00 1.181430e+00
#> [2811] 1.183143e+00 1.184860e+00 1.186583e+00 1.188311e+00 1.228365e+00
#> [2816] 1.230086e+00 1.231812e+00 1.233543e+00 1.235280e+00 1.237021e+00
#> [2821] 1.238768e+00 1.240520e+00 1.242276e+00 1.244038e+00 1.245805e+00
#> [2826] 1.247577e+00 1.284153e+00 1.285918e+00 1.287687e+00 1.289462e+00
#> [2831] 1.291241e+00 1.293026e+00 1.294815e+00 1.296610e+00 1.298409e+00
#> [2836] 1.300214e+00 1.302023e+00 1.303837e+00 1.295089e+00 1.296919e+00
#> [2841] 1.298754e+00 1.300593e+00 1.302438e+00 1.304287e+00 1.306141e+00
#> [2846] 1.308000e+00 1.309864e+00 1.311733e+00 1.313607e+00 1.315485e+00
#> [2851] 1.347924e+00 1.349795e+00 1.351671e+00 1.353551e+00 1.355436e+00
#> [2856] 1.357326e+00 1.359220e+00 1.361119e+00 1.363023e+00 1.364932e+00
#> [2861] 1.366845e+00 1.368762e+00 1.401932e+00 1.403840e+00 1.405752e+00
#> [2866] 1.407668e+00 1.409589e+00 1.411515e+00 1.413445e+00 1.415379e+00
#> [2871] 1.417318e+00 1.419261e+00 1.421209e+00 1.423161e+00 1.485454e+00
#> [2876] 1.487373e+00 1.489295e+00 1.491222e+00 1.493153e+00 1.495088e+00
#> [2881] 1.497027e+00 1.498971e+00 1.500918e+00 1.502870e+00 1.504825e+00
#> [2886] 1.506785e+00 3.651816e-01 3.663730e-01 3.675683e-01 3.687673e-01
#> [2891] 3.699702e-01 3.711769e-01 3.723874e-01 3.736017e-01 3.748199e-01
#> [2896] 3.760419e-01 3.772678e-01 3.843120e-01 3.855603e-01 3.868124e-01
#> [2901] 3.880685e-01 3.893286e-01 3.905926e-01 3.918605e-01 3.931324e-01
#> [2906] 3.944083e-01 3.956882e-01 3.969720e-01 3.982599e-01 4.048688e-01
#> [2911] 4.061779e-01 4.074910e-01 4.088082e-01 4.101295e-01 4.114548e-01
#> [2916] 4.127843e-01 4.141178e-01 4.154555e-01 8.556008e+00 4.181431e-01
#> [2921] 4.194931e-01 4.268596e-01 4.282326e-01 4.296098e-01 4.309912e-01
#> [2926] 4.323768e-01 4.337666e-01 4.351607e-01 4.365590e-01 4.379615e-01
#> [2931] 4.393682e-01 4.407793e-01 4.421945e-01 1.190807e-01 1.195031e-01
#> [2936] 1.199271e-01 1.203525e-01 1.207794e-01 1.212077e-01 1.216376e-01
#> [2941] 1.220690e-01 5.546317e-01 1.229362e-01 1.233721e-01 1.243533e-01
#> [2946] 1.247941e-01 1.252364e-01 1.256803e-01 1.261257e-01 1.265726e-01
#> [2951] 1.270212e-01 1.274712e-01 1.279229e-01 1.283760e-01 1.288308e-01
#> [2956] 1.292872e-01 1.311866e-01 1.316511e-01 1.321173e-01 1.325850e-01
#> [2961] 1.330544e-01 1.335254e-01 1.339980e-01 1.344722e-01 1.349481e-01
#> [2966] 1.354257e-01 1.359049e-01 1.363857e-01 1.440366e+00 1.444846e+00
#> [2971] 1.449341e+00 1.453849e+00 1.458372e+00 1.462908e+00 1.467458e+00
#> [2976] 1.472022e+00 1.476600e+00 1.481192e+00 1.485799e+00 4.220684e+00
#> [2981] 3.065100e+00 1.521687e+00 1.526395e+00 1.531117e+00 1.535854e+00
#> [2986] 1.540604e+00 1.545370e+00 1.550149e+00 1.554943e+00 1.559752e+00
#> [2991] 1.564575e+00 1.591991e+00 1.596888e+00 1.601800e+00 1.606727e+00
#> [2996] 1.611669e+00 1.616625e+00 1.621596e+00 1.626582e+00 1.631582e+00
#> [3001] 1.636597e+00 1.641627e+00 1.646672e+00 1.674695e+00 1.679815e+00
#> [3006] 1.684950e+00 1.690100e+00 1.695266e+00 1.700446e+00 1.705641e+00
#> [3011] 1.710852e+00 1.716077e+00 1.721318e+00 1.726574e+00 1.731845e+00
#> [3016] 1.761096e+00 1.766444e+00 1.771807e+00 1.777185e+00 1.782579e+00
#> [3021] 1.787989e+00 3.007363e+01 1.798854e+00 3.007204e+01 3.007091e+01
#> [3026] 4.922555e+01 3.006797e+01 1.838955e+00 1.844512e+00 1.850085e+00
#> [3031] 1.855674e+00 3.979746e+01 1.866897e+00 1.872532e+00 1.878183e+00
#> [3036] 3.871594e+01 4.212775e+01 3.750174e+01 6.248327e+00 5.018581e+01
#> [3041] 4.594369e+01 2.971318e+00 4.492800e+00 3.878719e+00 4.391067e+01
#> [3046] 1.455553e+01 9.510492e+00 1.978903e+00 8.903038e+00 3.683421e+00
#> [3051] 4.586983e+00 4.099297e+01 2.034522e+00 2.040567e+00 2.046629e+00
#> [3056] 2.052706e+00 2.058800e+00 2.064909e+00 2.071035e+00 2.077176e+00
#> [3061] 2.083334e+00 2.089508e+00 2.095697e+00 2.112198e+00 2.118438e+00
#> [3066] 2.124694e+00 2.130966e+00 2.137254e+00 2.143558e+00 2.149878e+00
#> [3071] 2.156215e+00 2.162567e+00 2.168935e+00 2.175320e+00 2.181720e+00
#> [3076] 2.216460e+00 2.222941e+00 2.229439e+00 2.235952e+00 4.613193e+01
#> [3081] 3.618939e+01 1.023583e+01 2.262165e+00 3.484704e+01 3.968366e+01
#> [3086] 2.991242e+01 2.288636e+00 2.995973e+01 2.331014e+00 2.337752e+00
#> [3091] 2.344505e+00 2.351274e+00 2.358059e+00 2.364861e+00 2.371678e+00
#> [3096] 2.378511e+00 2.385360e+00 2.392225e+00 2.399106e+00 2.436351e+00
#> [3101] 2.443311e+00 2.450287e+00 2.457279e+00 2.464287e+00 2.471310e+00
#> [3106] 2.478349e+00 2.485404e+00 2.492475e+00 2.499561e+00 2.506663e+00
#> [3111] 2.513780e+00 1.025649e+00 1.029097e+00 1.032557e+00 1.036027e+00
#> [3116] 1.039509e+00 7.386207e+00 1.046506e+00 1.050021e+00 1.053548e+00
#> [3121] 6.474137e+00 1.060635e+00 1.074371e+00 1.077972e+00 1.081585e+00
#> [3126] 1.085210e+00 1.088845e+00 1.092493e+00 1.096152e+00 1.099823e+00
#> [3131] 1.103506e+00 1.107200e+00 1.110906e+00 1.114624e+00 1.128773e+00
#> [3136] 1.132544e+00 1.136327e+00 1.140122e+00 1.143929e+00 1.147747e+00
#> [3141] 1.151578e+00 1.155421e+00 1.159277e+00 1.163144e+00 1.167023e+00
#> [3146] 1.170915e+00 1.173993e+00 1.177907e+00 1.181833e+00 1.185771e+00
#> [3151] 1.189722e+00 1.193685e+00 1.197660e+00 1.201648e+00 1.205648e+00
#> [3156] 1.209660e+00 1.213685e+00 1.217723e+00 1.327133e+00 1.331366e+00
#> [3161] 1.335612e+00 1.339871e+00 1.344144e+00 1.348430e+00 1.352729e+00
#> [3166] 1.357042e+00 3.087248e+00 1.365707e+00 1.370060e+00 1.388705e+00
#> [3171] 1.393118e+00 1.397544e+00 1.401983e+00 1.406437e+00 1.410904e+00
#> [3176] 1.415385e+00 1.419879e+00 1.424387e+00 1.428910e+00 1.433445e+00
#> [3181] 1.437995e+00 1.456818e+00 1.461428e+00 1.466051e+00 1.470689e+00
#> [3186] 1.475341e+00 1.480007e+00 1.484687e+00 1.489381e+00 1.494090e+00
#> [3191] 1.498812e+00 1.503549e+00 1.508300e+00 1.368354e+01 1.366338e+01
#> [3196] 1.364320e+01 1.362299e+01 2.487230e+01 1.112325e+00 1.114058e+00
#> [3201] 1.115796e+00 1.117539e+00 1.119287e+00 1.121041e+00 1.157917e+00
#> [3206] 1.159668e+00 1.161425e+00 1.163186e+00 1.164953e+00 1.166725e+00
#> [3211] 6.662340e+00 1.170284e+00 1.172071e+00 1.173864e+00 1.175661e+00
#> [3216] 1.177464e+00 1.141372e+00 1.143199e+00 1.145031e+00 1.146868e+00
#> [3221] 1.148710e+00 1.150557e+00 1.152409e+00 1.154267e+00 1.156129e+00
#> [3226] 1.301312e+01 1.159869e+00 1.161747e+00 1.197741e+00 1.199616e+00
#> [3231] 2.467992e+00 1.203381e+00 1.205271e+00 1.207166e+00 1.209066e+00
#> [3236] 1.210971e+00 1.212881e+00 1.214796e+00 1.216716e+00 1.218641e+00
#> [3241] 1.255817e+00 1.632575e+00 1.259660e+00 1.261589e+00 1.263522e+00
#> [3246] 1.265461e+00 1.267404e+00 1.269353e+00 1.271306e+00 1.273264e+00
#> [3251] 1.275226e+00 1.277194e+00 1.316622e+00 1.318580e+00 1.320542e+00
#> [3256] 1.322509e+00 1.324481e+00 1.326457e+00 1.328438e+00 1.330424e+00
#> [3261] 1.332414e+00 1.334409e+00 1.336408e+00 1.338412e+00 1.374332e+00
#> [3266] 1.376325e+00 1.378323e+00 1.380325e+00 1.382331e+00 1.384342e+00
#> [3271] 1.386356e+00 1.388376e+00 1.390399e+00 1.392427e+00 2.302576e+00
#> [3276] 1.396496e+00 1.388257e+00 1.390308e+00 1.392364e+00 1.394424e+00
#> [3281] 1.396488e+00 1.398556e+00 1.400628e+00 1.402704e+00 1.404785e+00
#> [3286] 1.406869e+00 1.408958e+00 1.411050e+00 1.442818e+00 3.931075e+00
#> [3291] 1.446983e+00 2.375757e+00 1.451164e+00 8.043719e+00 2.570658e+00
#> [3296] 8.026533e+00 4.708245e+00 1.461682e+00 1.463796e+00 1.465915e+00
#> [3301] 1.498281e+00 1.500384e+00 1.502490e+00 1.504600e+00 1.506713e+00
#> [3306] 1.508830e+00 1.510950e+00 1.513073e+00 1.515200e+00 1.517330e+00
#> [3311] 1.519463e+00 1.521599e+00 1.581888e+00 1.583981e+00 1.586077e+00
#> [3316] 1.588177e+00 1.590279e+00 1.592384e+00 1.594491e+00 1.596602e+00
#> [3321] 1.598715e+00 1.600831e+00 1.602950e+00 1.605071e+00 1.645621e+00
#> [3326] 1.647710e+00 1.170793e+01 3.431365e+00 1.653993e+00 1.162879e+01
#> [3331] 1.658194e+00 1.660297e+00 1.155004e+01 1.664510e+00 1.149775e+01
#> [3336] 1.147166e+01 1.155447e+01 1.712635e+00 1.714711e+00 1.716788e+00
#> [3341] 1.718867e+00 1.720947e+00 1.723029e+00 1.725113e+00 1.727198e+00
#> [3346] 1.131532e+01 1.128899e+01 1.733463e+00 1.133963e+01 1.131271e+01
#> [3351] 1.128583e+01 1.125901e+01 1.784083e+00 1.786132e+00 1.117886e+01
#> [3356] 1.115225e+01 1.112569e+01 1.109918e+01 1.107273e+01 1.798444e+00
#> [3361] 1.841545e+00 1.843547e+00 1.845550e+00 1.847552e+00 1.849555e+00
#> [3366] 1.851559e+00 1.853562e+00 1.855566e+00 1.857569e+00 1.859573e+00
#> [3371] 1.861577e+00 1.863580e+00 1.907381e+00 1.086413e+01 1.083718e+01
#> [3376] 1.081030e+01 1.915155e+00 1.075672e+01 1.073002e+01 1.070338e+01
#> [3381] 1.922919e+00 1.924858e+00 1.926796e+00 1.928733e+00 1.972208e+00
#> [3386] 1.063109e+01 1.060423e+01 1.057744e+01 1.055072e+01 1.052406e+01
#> [3391] 1.049746e+01 1.047093e+01 1.044446e+01 1.041805e+01 1.990854e+00
#> [3396] 1.036543e+01 1.694145e+00 1.699085e+00 1.704041e+00 1.709012e+00
#> [3401] 1.713997e+00 1.718998e+00 1.724013e+00 1.729044e+00 1.734090e+00
#> [3406] 1.739151e+00 1.744227e+00 1.721389e+00 1.726450e+00 1.731527e+00
#> [3411] 1.736619e+00 1.741726e+00 1.746848e+00 1.751986e+00 1.757138e+00
#> [3416] 1.762306e+00 1.767489e+00 1.772688e+00 1.777902e+00 1.852003e+00
#> [3421] 1.857356e+00 1.862724e+00 1.868108e+00 1.443573e+00 1.407262e+00
#> [3426] 1.340235e+00 1.335575e+00 1.338894e+00 1.334053e+00 1.336864e+00
#> [3431] 1.340195e+00 1.369435e+00 1.385891e+00 1.378924e+00 1.397645e+00
#> [3436] 1.399241e+00 1.395641e+00 1.405218e+00 1.408771e+00 1.379374e+00
#> [3441] 1.382805e+00 1.425150e+00 1.506922e+00 1.825612e+00 1.830741e+00
#> [3446] 1.835887e+00 1.841047e+00 4.656677e+01 4.230485e+01 3.206462e+01
#> [3451] 1.861844e+00 3.206178e+01 1.872335e+00 4.728972e+00 5.068945e+01
#> [3456] 1.907814e+00 1.913156e+00 1.918513e+00 1.923886e+00 1.929275e+00
#> [3461] 1.934680e+00 1.940100e+00 1.945536e+00 1.950988e+00 5.063959e+00
#> [3466] 4.563900e+00 1.987060e+00 1.992601e+00 1.998157e+00 2.003730e+00
#> [3471] 2.009319e+00 2.014923e+00 2.020544e+00 2.026180e+00 2.031833e+00
#> [3476] 2.037502e+00 2.043186e+00 2.048887e+00 2.074812e+00 2.080570e+00
#> [3481] 2.086345e+00 2.092136e+00 2.097944e+00 2.103767e+00 2.109607e+00
#> [3486] 2.115463e+00 2.121335e+00 2.127223e+00 2.133128e+00 2.139049e+00
#> [3491] 3.276888e+01 3.276524e+01 3.276139e+01 3.275731e+01 3.275300e+01
#> [3496] 1.950962e+00 1.956293e+00 1.961639e+00 3.273359e+01 4.948790e+01
#> [3501] 5.700053e+00 3.286033e+01 3.285346e+01 3.284637e+01 2.021752e+00
#> [3506] 3.283154e+01 2.032765e+00 2.038296e+00 2.043843e+00 3.341582e+01
#> [3511] 2.104244e+00 2.110001e+00 2.115774e+00 6.086485e+01 7.122032e+00
#> [3516] 7.124988e+00 7.127958e+00 7.130940e+00 7.133936e+00 5.994503e+01
#> [3521] 7.139965e+00 8.913517e+01 8.635729e+01 7.040017e+01 6.543999e+01
#> [3526] 7.853773e+01 6.012278e+01 5.996538e+01 7.583806e+01 7.452909e+00
#> [3531] 7.970764e+00 8.000497e+00 6.208684e+01 8.167127e+00 6.287851e+01
#> [3536] 9.207222e+00 9.867016e+00 7.004143e+01 1.002219e+01 6.983547e+01
#> [3541] 1.194000e+01 1.011182e+01 2.723447e+01 8.615488e+00 7.976417e+01
#> [3546] 6.721962e+00 5.738321e+01 2.016991e+01 6.960208e+00 6.951731e+00
#> [3551] 6.871212e+00 6.915006e+00 6.906173e+00 6.933338e+00 6.937686e+00
#> [3556] 7.202549e+00 7.263027e+00 1.662284e+00 1.234564e+00 1.179321e+00
#> [3561] 1.128444e+00 1.124395e+00 1.112007e+00 1.122978e+00 1.074901e+00
#> [3566] 1.073186e+00 9.868220e-01 9.780271e-01 5.276632e-01 4.186819e-01
#> [3571] 3.486501e-01 3.561979e-01 7.213050e+00 3.286978e+00 2.553053e+01
#> [3576] 3.286160e+00 2.536034e+01 3.285352e+00 3.284951e+00 3.284553e+00
#> [3581] 1.114689e+01 5.433331e+00 2.485518e+01 3.455944e+00 3.455328e+00
#> [3586] 3.454714e+00 3.454101e+00 3.453491e+00 3.452882e+00 3.452275e+00
#> [3591] 3.451669e+00 3.451066e+00 3.450464e+00 3.449864e+00 3.449265e+00
#> [3596] 3.623518e+00 3.622678e+00 3.621840e+00 3.621002e+00 3.620166e+00
#> [3601] 3.619331e+00 3.618497e+00 3.617664e+00 3.616832e+00 3.616001e+00
#> [3606] 3.615171e+00 3.614342e+00 3.287391e+00 3.286978e+00 2.553053e+01
#> [3611] 3.286160e+00 2.536034e+01 3.285352e+00 4.245584e+00 3.284553e+00
#> [3616] 6.254891e+00 3.283764e+00 2.485518e+01 3.455944e+00 3.455328e+00
#> [3621] 3.454714e+00 3.454101e+00 3.453491e+00 3.452882e+00 3.452275e+00
#> [3626] 3.451669e+00 2.488910e+01 3.450464e+00 3.449864e+00 3.449265e+00
#> [3631] 3.623518e+00 3.622678e+00 3.621840e+00 3.621002e+00 3.620166e+00
#> [3636] 3.619331e+00 3.618497e+00 3.617664e+00 3.616832e+00 3.616001e+00
#> [3641] 3.615171e+00 3.614342e+00 2.406017e+01 1.485786e+01 1.483107e+01
#> [3646] 1.480427e+01 1.477747e+01 1.475068e+01 1.472388e+01 1.469709e+01
#> [3651] 1.467030e+01 1.464351e+01 1.461672e+01 3.064215e+00 3.069989e+00
#> [3656] 3.075781e+00 3.081590e+00 3.087416e+00 3.093259e+00 3.099118e+00
#> [3661] 3.104995e+00 3.110889e+00 1.016085e+01 3.122726e+00 7.609882e+00
#> [3666] 5.930107e+00 5.469726e+00 5.948417e+00 3.243103e+00 3.249125e+00
#> [3671] 3.255165e+00 4.954338e+00 4.260359e+00 4.970821e+00 8.456573e+00
#> [3676] 3.285613e+00 3.324198e+00 3.330351e+00 3.336521e+00 3.342708e+00
#> [3681] 3.348912e+00 3.355132e+00 3.361370e+00 3.367624e+00 3.373894e+00
#> [3686] 3.380182e+00 3.386486e+00 3.392807e+00 3.563531e+00 3.569855e+00
#> [3691] 3.576196e+00 3.582553e+00 3.588927e+00 3.595318e+00 3.601724e+00
#> [3696] 3.608147e+00 3.614587e+00 3.621042e+00 3.627514e+00 3.634003e+00
#> [3701] 3.687249e+00 3.692890e+00 3.698547e+00 3.704222e+00 3.709913e+00
#> [3706] 3.715622e+00 3.721346e+00 3.727088e+00 7.308165e+00 1.546427e+01
#> [3711] 4.879272e+00 3.856401e+00 3.862185e+00 3.867985e+00 3.873802e+00
#> [3716] 3.879636e+00 3.885486e+00 3.891352e+00 3.897235e+00 3.903135e+00
#> [3721] 4.280930e+01 3.914983e+00 3.920932e+00 4.046135e+00 4.052063e+00
#> [3726] 4.058008e+00 4.063968e+00 4.069945e+00 4.075938e+00 4.081947e+00
#> [3731] 4.087973e+00 4.094014e+00 4.100071e+00 4.106145e+00 4.112234e+00
#> [3736] 5.175898e+00 1.516274e+01 1.005473e+01 6.743658e+00 5.194990e+00
#> [3741] 1.348393e+01 1.348727e+01 9.449236e+00 5.214330e+00 1.316694e+01
#> [3746] 8.738897e+00
# }
```
