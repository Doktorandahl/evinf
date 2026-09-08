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

#' @noRd
evinf_design <- function(formula, data) {
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  mt <- attr(mf, "terms")
  X <- stats::model.matrix(mt, mf)
  int <- match("(Intercept)", colnames(X), nomatch = 0L)
  if (int > 0L) {
    X <- X[, -int, drop = FALSE]
  }
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
  X
}
