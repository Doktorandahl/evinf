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

# audit0.10 §1.9: run_evzinb()/run_evinb() call this right after na.omit(), so
# the C++ side never sees a response it silently mishandles: a non-integer y
# is treated as ceiling(y) by the NB loop but exactly by the Pareto part, and
# a negative y makes em_c_candidates()'s while() condition error cryptically
# ("argument is of length zero").
evinf_check_response <- function(y) {
  bad_negative <- which(y < 0)
  if (length(bad_negative)) {
    stop(
      "The response contains ", length(bad_negative),
      " negative value(s); evzinb()/evinb() require a non-negative integer count.",
      call. = FALSE
    )
  }
  bad_noninteger <- which(!is.finite(y) | abs(y - round(y)) > sqrt(.Machine$double.eps))
  if (length(bad_noninteger)) {
    stop(
      "The response contains ", length(bad_noninteger),
      " non-finite or non-integer value(s); evzinb()/evinb() require a ",
      "non-negative integer count.",
      call. = FALSE
    )
  }
  invisible(y)
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
  off <- stats::model.offset(mf)
  list(
    X = X,
    terms = mt,
    xlevels = stats::.getXlevels(mt, mf),
    has_offset = !is.null(attr(mt, "offset")),
    offset = if (is.null(off)) rep(0, nrow(X)) else as.numeric(off)
  )
}

# Resolve a non-NB component formula (audit 4.4; round9 D.1):
#  * NULL           -> the NB formula, with any offset() term stripped (no
#                      inheritance -- an offset applies only where written
#                      explicitly);
#  * has offset()    -> kept when allow_offset is TRUE (formula_zi/formula_evi);
#                       an error when FALSE (formula_pareto: an offset on a
#                       shape parameter has no clear reading);
#  * one-sided       -> made two-sided with the NB response.
evinf_component_formula <- function(f, nb_formula, arg_name, allow_offset = TRUE) {
  has_offset <- function(x) !is.null(attr(stats::terms(x), "offset"))
  strip_offset <- function(x) {
    tt <- stats::terms(x)
    tl <- attr(tt, "term.labels")
    rhs <- if (length(tl)) paste(tl, collapse = " + ") else "1"
    if (attr(tt, "intercept") == 0L) {
      rhs <- paste(rhs, "- 1")
    }
    lhs <- if (length(x) == 3L) paste(deparse(x[[2]]), collapse = " ") else ""
    stats::as.formula(paste(lhs, "~", rhs), env = environment(x))
  }

  if (is.null(f)) {
    f <- nb_formula
    if (has_offset(f)) {
      f <- strip_offset(f)
    }
    return(f)
  }

  if (has_offset(f) && !allow_offset) {
    stop("offset() terms are not supported in ", arg_name, ": an offset on ",
         "a shape parameter has no clear reading. offset() is supported in ",
         "formula_nb, formula_zi and formula_evi.", call. = FALSE)
  }

  if (length(f) == 2L) {
    rhs <- paste(deparse(f[[2]]), collapse = " ")
    f <- stats::as.formula(
      paste(paste(deparse(nb_formula[[2]]), collapse = " "), "~", rhs),
      env = environment(f)
    )
  }
  f
}

# Offset vector for newdata from a stored component terms object (audit 4.4).
# NULL when the component formula has no offset() term; informative error when it
# does but the offset variable is absent from newdata.
evinf_offset_newdata <- function(terms, newdata) {
  if (is.null(attr(terms, "offset"))) {
    return(NULL)
  }
  terms <- stats::delete.response(terms)
  mf <- tryCatch(
    stats::model.frame(terms, newdata, na.action = stats::na.pass),
    error = function(e) {
      stop("The offset variable is missing from `newdata`.", call. = FALSE)
    }
  )
  off <- stats::model.offset(mf)
  if (is.null(off)) {
    return(NULL)
  }
  as.numeric(off)
}

# A component's design matrix never includes the offset, so drop the offset()
# term before building it from newdata (otherwise model.frame() would demand the
# offset variable be present even for design-only calls).
evinf_drop_offset_terms <- function(tt) {
  if (is.null(attr(tt, "offset"))) {
    return(tt)
  }
  tl <- attr(tt, "term.labels")
  rhs <- if (length(tl)) paste(tl, collapse = " + ") else "1"
  if (attr(tt, "intercept") == 0L) {
    rhs <- paste(rhs, "- 1")
  }
  stats::terms(stats::as.formula(paste("~", rhs)))
}

#' @noRd
evinf_design_newdata <- function(terms, xlevels, newdata) {
  terms <- evinf_drop_offset_terms(stats::delete.response(terms))
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
