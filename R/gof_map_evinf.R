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
#' @export
#'
#' @examples
#' gof_map_evinf()
#' gof_map_evinf(extra = data.frame(raw = "n_em_steps", clean = "EM steps",
#'                                  fmt = 0))
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
