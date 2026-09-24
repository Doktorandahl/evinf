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
#' @return A tibble with columns \code{id}, \code{type} (\code{"error"},
#'   \code{"degenerate"} or \code{"not_converged"} -- round10 G.3: a replicate
#'   that ran without erroring and isn't degenerate, but whose \code{converge}
#'   is \code{FALSE}, was previously invisible here), \code{message} (the
#'   error text, the degeneracy reason, or \code{NA} for
#'   \code{"not_converged"} -- see \code{n_em_steps}/\code{c_converged}/
#'   \code{c_warmup_capped} there instead) and, for \code{"not_converged"}
#'   rows, \code{n_em_steps}, \code{c_converged}, \code{c_warmup_capped}
#'   (\code{NA} for \code{"error"}/\code{"degenerate"} rows); zero rows when
#'   every replicate is usable and converged.
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
                          message = character(), n_em_steps = integer(),
                          c_converged = logical(), c_warmup_capped = logical())
  b <- object$bootstraps
  if (is.null(b)) {
    return(empty)
  }
  ids <- if (is.null(names(b))) as.character(seq_along(b)) else names(b)
  is_err <- vapply(b, inherits, logical(1), "try-error")
  is_deg <- vapply(b, function(x)
    !inherits(x, "try-error") && isTRUE(x$degenerate), logical(1))
  # round10 G.3: ran, isn't degenerate, but the EM itself didn't converge --
  # previously invisible to this function (neither an error nor degenerate).
  is_nc <- !is_err & !is_deg & vapply(b, function(x) !isTRUE(x$converge), logical(1))

  err <- tibble::tibble(
    id = ids[is_err], type = "error",
    message = vapply(b[is_err], function(e) trimws(as.character(e)), character(1)),
    n_em_steps = NA_integer_, c_converged = NA, c_warmup_capped = NA
  )
  deg <- tibble::tibble(
    id = ids[is_deg], type = "degenerate",
    message = vapply(b[is_deg], function(x)
      x$degenerate_reason %||% NA_character_, character(1)),
    n_em_steps = NA_integer_, c_converged = NA, c_warmup_capped = NA
  )
  nc <- tibble::tibble(
    id = ids[is_nc], type = "not_converged", message = NA_character_,
    n_em_steps = vapply(b[is_nc], function(x) x$n_em_steps %||% NA_integer_, integer(1)),
    c_converged = vapply(b[is_nc], function(x) isTRUE(x$c_converged), logical(1)),
    c_warmup_capped = vapply(b[is_nc], function(x) isTRUE(x$c_warmup_capped), logical(1))
  )
  dplyr::bind_rows(err, deg, nc)
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
