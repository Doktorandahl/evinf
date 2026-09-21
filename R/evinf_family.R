#' Model family for evzinb() / evinb() (round9 E.0)
#'
#' Specifies the distributional family of the count state (\code{count}) and
#' of the zero state (\code{zero}). One object, so later additions (a
#' truncated-count state, a geometric count) do not add more arguments to
#' \code{\link{evzinb}()} / \code{\link{evinb}()}.
#'
#' @param count \code{"nbinom"} (the default; a negative-binomial count
#'   state) or \code{"poisson"} (a Poisson count state -- \code{Alpha.NB} is
#'   dropped entirely: not in \code{par.all}, not in \code{coef()}/
#'   \code{vcov()}/\code{confint()}/\code{tidy()}, shown as absent (not
#'   \code{NA}) in \code{summary()}/\code{glance()}).
#' @param zero \code{"mixture"} (the default; the existing zero-inflation
#'   mixture, where a zero can come from either the zero state or the count
#'   state) or \code{"hurdle"} (the zero state owns every zero and the count
#'   state is zero-truncated; \code{evzinb()} only -- \code{evinb()} has no
#'   zero state to hurdle over and errors if asked for one).
#'
#' @return An object of class \code{"evinf_family"}: a list with elements
#'   \code{count} and \code{zero}.
#' @export
#'
#' @examples
#' evinf_family()
#' evinf_family(count = "poisson")
evinf_family <- function(count = c("nbinom", "poisson"),
                         zero = c("mixture", "hurdle")) {
  count <- match.arg(count)
  zero <- match.arg(zero)
  structure(list(count = count, zero = zero), class = "evinf_family")
}

# evzinb()/evinb()'s family = argument accepts an evinf_family() object or a
# plain string (shorthand for evinf_family(count = string)).
evinf_resolve_family <- function(family) {
  if (inherits(family, "evinf_family")) {
    return(family)
  }
  if (is.character(family) && length(family) == 1L) {
    return(evinf_family(count = family))
  }
  stop(
    "`family` must be an evinf_family() object or a single string naming ",
    "the count family (\"nbinom\" or \"poisson\").", call. = FALSE
  )
}

# Integer codes for the C++ side (round9 E.1/E.2): 0 = nbinom/mixture (the
# default, reproducing today's model exactly), 1 = poisson/hurdle.
evinf_family_count_code <- function(family) if (family$count == "poisson") 1L else 0L
evinf_family_zero_code <- function(family) if (family$zero == "hurdle") 1L else 0L

#' @export
print.evinf_family <- function(x, ...) {
  cat("<evinf_family>\n")
  cat("  count: ", x$count, "\n", sep = "")
  cat("  zero:  ", x$zero, "\n", sep = "")
  invisible(x)
}
