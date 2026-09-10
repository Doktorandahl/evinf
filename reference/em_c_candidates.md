# Candidate set for the extreme-value threshold \\C\_{EV}\\

Builds the set of values over which the EM driver profiles the
log-likelihood for the extreme-value threshold \\C\_{EV}\\: the unique
observed response values inside `c.lim`, optionally thinned.

## Usage

``` r
em_c_candidates(y, c.lim, prune.c.range)
```

## Arguments

- y:

  Numeric response vector.

- c.lim:

  Length-2 numeric giving the lower and upper bound of the search.

- prune.c.range:

  `FALSE` for no thinning, or a number in `[0, 1]` giving the proportion
  of interior candidates to drop (sampled with probability proportional
  to the gaps between consecutive candidates, so the endpoints are
  always kept).

## Value

A sorted numeric vector of candidate \\C\_{EV}\\ values.

## Details

This is the ECME step of the algorithm in Appendix A1 of Randahl and
Vegelius (2024): \\C\_{EV}\\ is updated by a grid search over the
observed support rather than by a smooth optimiser. A warning is emitted
when the set has more than 100 elements and no pruning was requested,
because each extra candidate is one extra full-data log-likelihood
evaluation per EM iteration.

## See also

[`evzinb()`](evzinb.md), [`evinb()`](evinb.md)
