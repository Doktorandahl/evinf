## Regenerates data/hks.rda
##
## The only transformation applied here is the 0.10.0 rename of the
## log-transformed battle-related-deaths column `brv_AllLag` -> `brv_AllLag_log`,
## so that it follows the `_log` suffix convention of the other transformed
## columns (`troopLag_log`, `epdur_log`, `policeLag_log`,
## `militaryobserversLag_log`). The underlying values are unchanged.
##
## The bundled data are the reduced replication set from Randahl and Vegelius
## (2024); see ?hks. They carry no conflict identifier, so the conflict-level
## cluster bootstrap in that paper cannot be reproduced from the bundled data
## alone.

load("data/hks.rda")

names(hks)[names(hks) == "brv_AllLag"] <- "brv_AllLag_log"

usethis::use_data(hks, overwrite = TRUE, version = 3)
