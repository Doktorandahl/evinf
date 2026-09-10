# EM iterations for the EVZINB / EVINB mixture at a fixed \\C\_{EV}\\

Runs [`em_step`](em_step.md) repeatedly, starting from `ini.val`, until
the largest absolute parameter change falls below `control$max.diff.par`
or `control$max.no.em.steps` iterations have been taken. The
extreme-value threshold is held at `ini.val$C` throughout; the driver
([`em_fit`](em_fit.md)) alternates calls to this routine with the ECME
update of \\C\_{EV}\\ ([`em_profile_c`](em_profile_c.md)).

## Usage

``` r
em_fit_fixed_c(y, x_obj, ini.val, control, fixed_zc = FALSE)
```

## Arguments

- y:

  Numeric response vector.

- x_obj:

  List of raw component design matrices (see
  [`em_extend_design`](em_extend_design.md)).

- ini.val:

  List of starting values: `Beta.multinom.ZC`, `Beta.multinom.PL`,
  `Beta.NB`, `Alpha.NB`, `Beta.PL`, `C`.

- control:

  An [`evinf_control`](evinf_control.md) object.

- fixed_zc:

  When `TRUE` (EVINB) the zero-inflation multinomial block is held at
  `ini.val$Beta.multinom.ZC`.

## Value

A list with

- par.mat:

  the estimated parameters plus `Props` (the E-step prior state
  probabilities) and `C`;

- par.all:

  the free parameters as a flat numeric vector;

- resp:

  the E-step responsibilities (posterior state probabilities);

- log.lik.vec:

  the log-likelihood after each EM step;

- log.lik:

  the final log-likelihood;

- converge:

  `TRUE` if the loop stopped on the tolerance, not the iteration cap;

- n_em_steps:

  the number of EM steps taken.

## See also

[`em_fit`](em_fit.md), [`evzinb()`](evzinb.md), [`evinb()`](evinb.md)
