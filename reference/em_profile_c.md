# Profile the log-likelihood over the extreme-value threshold

Evaluates the full-data log-likelihood at each candidate value of the
extreme-value threshold \\C\_{EV}\\, holding all other parameters fixed,
and returns the profile together with its maximiser. This is the ECME
update of \\C\_{EV}\\ in Appendix A1 of Randahl and Vegelius (2024); it
is run once per EM iteration, in both the warm-up and the convergence
phase.

## Usage

``` r
em_profile_c(y, x_obj, par, c_candidates)
```

## Arguments

- y:

  Numeric response vector.

- x_obj:

  List of component design matrices (see
  [`em_extend_design`](em_extend_design.md)).

- par:

  List of current parameter values with elements `Beta.multinom.ZC`,
  `Beta.multinom.PL`, `Beta.NB`, `Alpha.NB`, `Beta.PL`.

- c_candidates:

  Numeric vector of candidate \\C\_{EV}\\ values, from
  [`em_c_candidates`](em_c_candidates.md).

## Value

A list with

- profile:

  a data frame with columns `c` and `loglik`;

- c_hat:

  the candidate(s) attaining the maximum log-likelihood;

- loglik_max:

  that maximum.

## See also

[`evzinb()`](evzinb.md), [`evinb()`](evinb.md)
