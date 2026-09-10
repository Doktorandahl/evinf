#' Inspect the bootstrap replicates that failed or came out degenerate
#'
#' Bootstrap fits that error are kept as \code{try-error} objects rather than
#' discarded; fits that converge but whose extreme-value tail is so heavy that
#' summaries from them are effectively unbounded are flagged \emph{degenerate}
#' (see \code{alpha_floor} in \code{\link{evinf_control}}). This returns both,
#' which are useful diagnostics (e.g. a resample with a singular design, or one
#' whose Pareto shape collapsed).
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model.
#' @return A tibble with columns \code{id}, \code{type} (\code{"error"} or
#'   \code{"degenerate"}) and \code{message} (the error text, or the reason the
#'   replicate is degenerate); zero rows when every replicate is usable.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' failed_bootstraps(model)
#' }
failed_bootstraps <- function(object) {
  empty <- tibble::tibble(id = character(), type = character(),
                          message = character())
  b <- object$bootstraps
  if (is.null(b)) {
    return(empty)
  }
  ids <- if (is.null(names(b))) as.character(seq_along(b)) else names(b)
  is_err <- vapply(b, inherits, logical(1), "try-error")
  is_deg <- vapply(b, function(x)
    !inherits(x, "try-error") && isTRUE(x$degenerate), logical(1))

  err <- tibble::tibble(
    id = ids[is_err], type = "error",
    message = vapply(b[is_err], function(e) trimws(as.character(e)), character(1))
  )
  deg <- tibble::tibble(
    id = ids[is_deg], type = "degenerate",
    message = vapply(b[is_deg], function(x)
      x$degenerate_reason %||% NA_character_, character(1))
  )
  dplyr::bind_rows(err, deg)
}

#' The bootstrap replicates usable for a summary
#'
#' Drops the \code{try-error} replicates and, unless
#' \code{exclude_degenerate = FALSE}, the ones flagged degenerate.
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model.
#' @param exclude_degenerate Also drop replicates flagged degenerate (default
#'   \code{TRUE}).
#' @return A (possibly empty) list of bootstrap fits.
#' @keywords internal
evinf_usable_bootstraps <- function(object, exclude_degenerate = TRUE) {
  b <- object$bootstraps
  if (is.null(b)) {
    return(list())
  }
  keep <- !vapply(b, inherits, logical(1), "try-error")
  if (isTRUE(exclude_degenerate)) {
    keep <- keep & !vapply(b, function(x)
      !inherits(x, "try-error") && isTRUE(x$degenerate), logical(1))
  }
  b[keep]
}
