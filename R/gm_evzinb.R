#' A goodness-of-fit gof tibble for GOF metrics when using modelsummary
#'
#' A goodness-of-fit gof tibble for GOF metrics when using modelsummary. The GM tibble can be used to obtain correct table output when making regression tables with modelsummary
#'
#' @format ## `gm_evzinb`
#' A tibble with 11 rows and 3 columns:
#' \describe{
#'   \item{raw}{The modelsummary/broom internal name for the statistic}
#'   \item{clean}{The table output for the statistic}
#'   \item{fmt}{The number of decimals reported for each statistic by default (can be adapted)}}
#' @seealso \code{\link{gof_map_evinf}()}, which builds this object (and lets you
#'   append rows); the regeneration script lives in `data-raw/gm_evzinb.R`.
"gm_evzinb"
