# Replication data for Hultman, Kathman, and Shannon (2013) United Nations Peacekeeping and Civilian Protection in Civil War

A reduced replication data set from Hultman et al. (2013) United Nations
Peacekeeping and Civilian Protection in Civil War. Used to reproduce the
results of Randahl and Vegelius (2024) . To reproduce any other results
from Hultman et al. (2013) please download the original replication
dataset using the link under source.

## Usage

``` r
hks
```

## Format

A tibble with 3746 rows and 13 columns:

- osvAll:

  The number of observed fatalities from one-sided violence against
  civilians in the specified conflict-month

- troopLag:

  The number of UN military troops, in thousands (lagged)

- policeLag:

  The number of UN police, in thousands (lagged)

- militaryobserversLag:

  The number of UN military observers, in thousands (lagged)

- brv_AllLag_log:

  The natural logarithm of the total number of battle-related deaths in
  the conflict in the previous month, log-transformed as in Hultman et
  al. (2013)

- osvAllLagDum:

  A dummy variable taking the value 1 if any one-sided violence against
  civilians took place in the previous conflict month

- incomp:

  UCDP/PRIO incompatibility: 1 = territory, 2 = government

- epduration:

  The number of months the current conflict-episode has been ongoing

- lntpop:

  The natural logarithm of the population of the country in which the
  conflict takes place

- troopLag_log:

  \`log1p(troopLag)\`

- epdur_log:

  \`log(epduration)\`

- policeLag_log:

  \`log1p(policeLag)\`

- militaryobserversLag_log:

  \`log1p(militaryobserversLag)\`

## Source

https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6EBCGA

## Details

The UN personnel counts (\`troopLag\`, \`policeLag\`,
\`militaryobserversLag\`) are expressed in \*\*thousands\*\* of
personnel. This is a rescaling of the original Hultman, Kathman, and
Shannon replication data introduced for Randahl and Vegelius (2024) .

The columns ending in \`\_log\` (\`troopLag_log\`, \`policeLag_log\`,
\`militaryobserversLag_log\`, \`epdur_log\`, \`brv_AllLag_log\`) are
transforms that were pre-computed for the Pareto component in Randahl
and Vegelius (2024) . The transform follows the standard rule in the
conflict literature: the personnel counts contain structural zeros and
use \`log1p()\`, while an episode duration is always at least one month
(at least two in these data) and uses \`log()\`, the conventional
transform for durations. Since evinf 0.9.4 the component formulas accept
in-formula transformations, so \`log1p(troopLag)\`, \`log(epduration)\`
etc. can be written directly in the formula and these pre-computed
columns are no longer needed.

The bundled data carry \*\*no conflict identifier\*\*. The
conflict-level cluster bootstrap reported in Randahl and Vegelius (2024)
therefore cannot be reproduced from the bundled data alone; download the
original replication data (link under \`source\`) and merge on your own
conflict key to use `block` in [`evzinb()`](evzinb.md) /
[`evinb()`](evinb.md).

## References

Hultman L, Kathman J, Shannon M (2013). “United Nations peacekeeping and
civilian protection in civil war.” *American Journal of Political
Science*, **57**(4), 875–891.  
  
Randahl D, Vegelius J (2024). “Inference with Extremes: Accounting for
Extreme Values in Count Regression Models.” *International Studies
Quarterly*, **68**(4), sqae137.
[doi:10.1093/isq/sqae137](https://doi.org/10.1093/isq/sqae137) .
