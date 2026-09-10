# EM / ECME driver for the EVZINB and EVINB models

Estimates an EVZINB or EVINB model by the generalised-EM / ECME
algorithm of Appendix A1 in Randahl and Vegelius (2024). It alternates
between EM iterations at a fixed extreme-value threshold
([`em_fit_fixed_c`](em_fit_fixed_c.md)) and a grid update of that
threshold ([`em_profile_c`](em_profile_c.md)), first in a short
"warm-up" phase (`control$max.no.em.steps.warmup` EM steps) and then in
a "convergence" phase that runs to `control$max.diff.par`.

## Usage

``` r
em_fit(y, x.obj, ini.val, control, model = c("evzinb", "evinb"))
```

## Arguments

- y:

  Numeric response vector.

- x.obj:

  List of raw component design matrices (elements `X.multinom.ZC`,
  `X.multinom.PL`, `X.NB`, `X.PL`, each a numeric matrix \\n \times p\\
  without an intercept column or `NULL`; optional `offset.nb`).

- ini.val:

  List of starting values (`Beta.multinom.ZC`, `Beta.multinom.PL`,
  `Beta.NB`, `Alpha.NB`, `Beta.PL`, `C`).

- control:

  An [`evinf_control`](evinf_control.md) object.

- model:

  `"evzinb"` (all three latent states free) or `"evinb"` (the
  zero-inflation component is switched off: its coefficients are held at
  `ini.val$Beta.multinom.ZC` in the convergence phase and in every
  \\C\_{EV}\\ profile).

## Value

A list with, among others, `par.mat` (estimated parameters), `log.lik`,
`AIC`, `BIC`, `resp` (posterior state probabilities), `converge`,
`c_profile`, `c_trace`, `log.lik.vec.all`, `loglik_recomputed`, the
fitted-value vectors (`mu.nb.vec`, `alpha.pl.vec`, `y.hat.pl*`, ...) and
the design matrices. The exact set and names are consumed by
`run_evzinb()` / `run_evinb()`.

## Details

For `model = "evinb"` the zero-inflation multinomial block is held at
its initial value (`ini.val$Beta.multinom.ZC`) throughout both phases
and in every \\C\_{EV}\\ profile.

## See also

[`evzinb()`](evzinb.md), [`evinb()`](evinb.md)
