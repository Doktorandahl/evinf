# round9 D.2 (audit §5.6): weights = for evzinb() / evinb().

# Resolve the `weights` argument, before na.omit(). Accepts, in order of
# preference (same anti-shadowing principle as evinf_block_name(), audit N5
# -- a bare symbol matching a data column is checked BEFORE evaluation, so a
# same-named object in the calling environment cannot shadow the column):
#   * NULL                                        -> NULL column name
#   * a bare symbol that names a column of `data` -> that column name
#   * a string naming a column of `data`          -> that column name
#   * a numeric vector                             -> injected into `data`
#     under a reserved column name, so it survives the same na.omit() step
#     and boot_id resampling as every other model variable
#
# Returns a list(data, weights_col): `data` is the input, or with the
# reserved column added when a raw vector was supplied; `weights_col` is the
# column name to use downstream (including for update()'s call reconstruction,
# same as `block`), or NULL when weights = NULL. `data` is returned even when
# unchanged so callers always reassign it uniformly.
evinf_resolve_weights <- function(quo, env = rlang::caller_env(), data) {
  if (rlang::quo_is_null(quo) || rlang::quo_is_missing(quo)) {
    return(list(data = data, weights_col = NULL))
  }
  expr <- rlang::quo_get_expr(quo)
  if (rlang::is_symbol(expr) && rlang::as_name(expr) %in% names(data)) {
    return(list(data = data, weights_col = rlang::as_name(expr)))
  }
  val <- tryCatch(rlang::eval_tidy(quo, env = env), error = function(e) NULL)
  if (is.character(val) && length(val) == 1L && val %in% names(data)) {
    return(list(data = data, weights_col = val))
  }
  if (is.numeric(val)) {
    evinf_check_weights(val, nrow(data))
    data[[".evinf_weights"]] <- val
    return(list(data = data, weights_col = ".evinf_weights"))
  }
  stop(
    "`weights` must be NULL, a bare column name, a string naming a column ",
    "of `data`, or a numeric vector.", call. = FALSE
  )
}

# Frequency/analytic weights must be positive and finite (non-integer values
# are allowed -- analytic weights, not just frequency counts).
evinf_check_weights <- function(w, n) {
  if (length(w) != n) {
    stop(
      "`weights` must have length equal to the number of rows of `data` (",
      n, "), not ", length(w), ".", call. = FALSE
    )
  }
  bad <- which(!is.finite(w) | w <= 0)
  if (length(bad)) {
    stop(
      "`weights` must be positive and finite; ", length(bad),
      " value(s) are not.", call. = FALSE
    )
  }
  invisible(w)
}
