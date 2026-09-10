# A goodness-of-fit gof tibble for GOF metrics when using modelsummary

A goodness-of-fit gof tibble for GOF metrics when using modelsummary.
The GM tibble can be used to obtain correct table output when making
regression tables with modelsummary

## Usage

``` r
gm_evzinb
```

## Format

\## \`gm_evzinb\` A tibble with 11 rows and 3 columns:

- raw:

  The modelsummary/broom internal name for the statistic

- clean:

  The table output for the statistic

- fmt:

  The number of decimals reported for each statistic by default (can be
  adapted)

## See also

[`gof_map_evinf()`](gof_map_evinf.md), which builds this object (and
lets you append rows); the regeneration script lives in
\`data-raw/gm_evzinb.R\`.
