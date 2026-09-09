# audit 4.13 / N5 - `block` accepts a bare column name as well as a string.

# Resolve the `block` argument of evzinb() / evinb() to a single column-name
# string (or NULL). Accepts, in order of preference:
#   * NULL                                        -> NULL
#   * a string literal ("id")                     -> "id"
#   * a bare symbol that names a column of `data` -> that column name
#     (checked BEFORE evaluation, so a same-named object in the calling
#      environment does not shadow the column -- audit N5)
#   * a string held in a variable  (b <- "id"; block = b)  -> that string
#   * anything else                               -> handed back for the usual
#                                                     "single string" error
evinf_block_name <- function(quo, env = rlang::caller_env(), data = NULL) {
  if (rlang::quo_is_null(quo) || rlang::quo_is_missing(quo)) {
    return(NULL)
  }
  expr <- rlang::quo_get_expr(quo)
  if (is.character(expr) && length(expr) == 1L) {
    return(expr)
  }
  if (rlang::is_symbol(expr) &&
      rlang::as_name(expr) %in% names(data)) {
    return(rlang::as_name(expr))
  }
  val <- tryCatch(rlang::eval_tidy(quo, env = env), error = function(e) NULL)
  if (is.character(val) && length(val) == 1L) {
    return(val)
  }
  if (is.null(val) && rlang::is_symbol(expr)) {
    return(rlang::as_name(expr))
  }
  val
}
