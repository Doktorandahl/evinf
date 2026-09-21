# One EM step for the EVZINB / EVINB mixture at fixed \\C\_{EV}\\

Performs a single generalised-EM iteration: the C++ routine
`update_bfgs_fun()` takes a bounded quasi-Newton step in all free
parameters, then each component's step is accepted only if a
one-dimensional line search over the damping factor \\\eta\\ does not
decrease the log-likelihood (Appendix A1 of Randahl and Vegelius, 2024).
The extreme-value threshold \\C\_{EV}\\ is held fixed.

## Usage

``` r
em_step(y, ext, par, control, fixed_zc = FALSE)
```

## Arguments

- y:

  Numeric response vector.

- ext:

  Extended design matrices from
  [`em_extend_design`](em_extend_design.md).

- par:

  List of current parameter values: `Beta.multinom.ZC`,
  `Beta.multinom.PL`, `Beta.NB`, `Alpha.NB`, `Beta.PL`, `C`.

- control:

  An [`evinf_control`](evinf_control.md) object (for `max.upd.par.nb`,
  `no.m.bfgs.steps.nb`, `eta.int`).

- fixed_zc:

  When `TRUE` (the EVINB case) the zero-inflation multinomial block is
  not updated and its line search is skipped.

## Value

A list with

- par:

  the updated parameter list (same names as the input);

- par_end:

  the updated free parameters as a flat numeric vector;

- max_abs_par_diff:

  the largest absolute parameter change this step;

- log_lik:

  the log-likelihood at the updated parameters;

- upd_obj:

  the raw return of `update_bfgs_fun()` (its `prop` and `resp` elements
  are the E-step quantities).

## See also

[`em_fit_fixed_c`](em_fit_fixed_c.md), [`evzinb()`](evzinb.md),
[`evinb()`](evinb.md)
