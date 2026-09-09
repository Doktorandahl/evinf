#' Replication data for Hultman, Kathman, and Shannon (2013) United Nations Peacekeeping and Civilian Protection in Civil War
#'
#' A reduced replication data set from \insertCite{hultman2013united;textual}{evinf} United Nations Peacekeeping and Civilian Protection in Civil War. Used to reproduce the results of \insertCite{randahl2024inference;textual}{evinf}. To reproduce any other results from \insertCite{hultman2013united;textual}{evinf} please download the original replication dataset using the link under source.
#'
#' The UN personnel counts (`troopLag`, `policeLag`, `militaryobserversLag`) are
#' expressed in **thousands** of personnel. This is a rescaling of the original
#' Hultman, Kathman, and Shannon replication data introduced for
#' \insertCite{randahl2024inference;textual}{evinf}.
#'
#' The columns ending in `_log` (`troopLag_log`, `policeLag_log`,
#' `militaryobserversLag_log`, `epdur_log`) are transforms that were
#' pre-computed for the Pareto component in
#' \insertCite{randahl2024inference;textual}{evinf}: the three personnel columns
#' use `log1p()` and `epdur_log` uses `log()`. Since evinf 0.9.4 the component
#' formulas accept in-formula transformations, so `log1p(troopLag)` etc. can be
#' used directly in the formula and these pre-computed columns are no longer
#' needed.
#'
#' @format
#' A tibble with 3746 rows and 13 columns:
#' \describe{
#'   \item{osvAll}{The number of observed fatalities from one-sided violence against civilians in the specified conflict-month}
#'   \item{troopLag}{The number of UN military troops, in thousands (lagged)}
#'   \item{policeLag}{The number of UN police, in thousands (lagged)}
#'   \item{militaryobserversLag}{The number of UN military observers, in thousands (lagged)}
#'   \item{brv_AllLag}{The natural logarithm of the total number of battle-related deaths in the conflict in the previous month}
#'   \item{osvAllLagDum}{A dummy variable taking the value 1 if any one-sided violence against civilians took place in the previous conflict month}
#'   \item{incomp}{UCDP/PRIO incompatibility: 1 = territory, 2 = government}
#'   \item{epduration}{The number of months the current conflict-episode has been ongoing}
#'   \item{lntpop}{The natural logarithm of the population of the country in which the conflict takes place}
#'   \item{troopLag_log}{`log1p(troopLag)`}
#'   \item{epdur_log}{`log(epduration)`}
#'   \item{policeLag_log}{`log1p(policeLag)`}
#'   \item{militaryobserversLag_log}{`log1p(militaryobserversLag)`}
#' }
#' @source https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6EBCGA
#' @references
#' \insertAllCited{}
"hks"
