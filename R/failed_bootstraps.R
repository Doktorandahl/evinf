#' Inspect the bootstrap replicates that failed
#'
#' Bootstrap fits that error are kept as \code{try-error} objects rather than
#' discarded; this returns their error messages, which are useful diagnostics
#' (e.g. a resample with a singular design).
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model.
#' @return A tibble with columns \code{id} and \code{message} (zero rows when
#'   nothing failed).
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' failed_bootstraps(model)
#' }
failed_bootstraps <- function(object) {
  b <- object$bootstraps
  if (is.null(b)) {
    return(tibble::tibble(id = character(), message = character()))
  }
  failed <- vapply(b, inherits, logical(1), "try-error")
  tibble::tibble(
    id = if (is.null(names(b))) which(failed) else names(b)[failed],
    message = vapply(b[failed], function(e) trimws(as.character(e)), character(1))
  )
}
