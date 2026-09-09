# audit 4.13 - the modelsummary goodness-of-fit map as a function.

#' Goodness-of-fit map for \code{modelsummary}
#'
#' Returns the \code{gof_map} that \code{modelsummary::modelsummary()} needs to
#' label and format the goodness-of-fit rows produced by \code{glance.evzinb()} /
#' \code{glance.evinb()}. The bundled data object \code{\link{gm_evzinb}} is
#' simply \code{gof_map_evinf()}.
#'
#' @param extra Optional data frame / tibble of additional rows (columns
#'   \code{raw}, \code{clean}, \code{fmt}) to append, e.g. for a statistic added
#'   by a custom \code{glance()} method.
#'
#' @return A tibble with columns \code{raw} (the \code{glance()} column name),
#'   \code{clean} (the row label in the table) and \code{fmt} (default number of
#'   decimals). Rows: number of observations, number of parameters, alpha_NB,
#'   C_EV, observations at or above C_EV, log-likelihood, AIC, BIC, number of
#'   successful and failed bootstraps, and convergence.
#'
#' @details The coefficient table produced by \code{tidy()} with
#'   \code{component = "all"} has one row per coefficient \emph{per component}.
#'   \code{modelsummary::modelsummary()} cannot render that with its default
#'   arguments (the same limitation it has for \code{nnet::multinom()}); pass
#'   \code{shape = term + y.level ~ model} together with
#'   \code{gof_map = gof_map_evinf()}.
#' @export
#'
#' @examples
#' gof_map_evinf()
#' gof_map_evinf(extra = data.frame(raw = "n_em_steps", clean = "EM steps",
#'                                  fmt = 0))
#'
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' if (requireNamespace("modelsummary", quietly = TRUE)) {
#'   modelsummary::modelsummary(
#'     model,
#'     shape = term + y.level ~ model,
#'     gof_map = gof_map_evinf()
#'   )
#' }
#' }
gof_map_evinf <- function(extra = NULL) {
  base <- tibble::tribble(
    ~raw,                  ~clean,                ~fmt,
    "nobs",                "obs",                 0,
    "npar",                "par",                 0,
    "alpha",               "alpha_nb",            2,
    "parameter",           "C_EV",                0,
    "n_above_c",           "obs_above_c_ev",      0,
    "logLik",              "logLik",              2,
    "aic",                 "AIC",                 1,
    "bic",                 "BIC",                 1,
    "n_bootstraps",        "n_bootstraps",        0,
    "n_failed_bootstraps", "n_failed_bootstraps", 0,
    "converged",           "converged",           0
  )
  if (!is.null(extra)) {
    base <- dplyr::bind_rows(base, tibble::as_tibble(extra))
  }
  base
}
