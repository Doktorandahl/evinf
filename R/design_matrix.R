# Design-matrix construction for the model components.
#
# The C++ estimation routines expect a numeric design matrix *without* an
# intercept column (they prepend the intercept themselves). Building these with
# stats::model.matrix() (rather than as.matrix(model.frame()[, -1])) means
# factors, interactions, poly()/ns() and in-formula transforms such as log(x)
# all work, and coefficient names come from colnames() of the matrix.
#
# The terms object and xlevels of each component are stored on the fitted object
# so that predictions on newdata use the same contrasts and factor levels.

# audit R0.8: an in-formula transformation (e.g. log(x) with x <= 0) can put
# -Inf / NaN into the design matrix and misbehave silently in the C++ code.
# Stop with an informative message naming the offending column(s) instead.
evinf_check_finite <- function(X, context = "design matrix") {
  if (!length(X)) {
    return(invisible(X))
  }
  bad <- colnames(X)[!apply(is.finite(X), 2L, all)]
  if (length(bad)) {
    stop(
      "Non-finite values in the ", context, " column(s): ",
      paste(bad, collapse = ", "),
      ". Check in-formula transformations (e.g. log() of non-positive values).",
      call. = FALSE
    )
  }
  invisible(X)
}

#' @noRd
evinf_design <- function(formula, data) {
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  mt <- attr(mf, "terms")
  X <- stats::model.matrix(mt, mf)
  int <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (int > 0L) {
    X <- X[, -int, drop = FALSE]
  }
  evinf_check_finite(X)
  list(
    X = X,
    terms = mt,
    xlevels = stats::.getXlevels(mt, mf)
  )
}

#' @noRd
evinf_design_newdata <- function(terms, xlevels, newdata) {
  terms <- stats::delete.response(terms)
  mf <- stats::model.frame(
    terms,
    newdata,
    xlev = xlevels,
    na.action = stats::na.pass
  )
  X <- stats::model.matrix(terms, mf, xlev = xlevels)
  int <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (int > 0L) {
    X <- X[, -int, drop = FALSE]
  }
  evinf_check_finite(X, context = "newdata design matrix")
  X
}
