# audit 4.13 - `block` accepts a bare column name as well as a string.

# Resolve the `block` argument of evzinb() / evinb() to a single column-name
# string (or NULL). Accepts, in order of preference:
#   * NULL                      -> NULL
#   * a string / string variable ("id", or  b <- "id"; block = b)  -> that string
#   * a bare column name         (block = id)                       -> "id"
evinf_block_name <- function(quo, env = rlang::caller_env()) {
  if (rlang::quo_is_null(quo) || rlang::quo_is_missing(quo)) {
    return(NULL)
  }
  expr <- rlang::quo_get_expr(quo)
  if (is.character(expr) && length(expr) == 1L) {
    return(expr)
  }
  val <- tryCatch(rlang::eval_tidy(quo, env = env), error = function(e) NULL)
  if (is.character(val) && length(val) == 1L) {
    return(val)
  }
  if (is.null(val) && rlang::is_symbol(expr)) {
    return(rlang::as_name(expr))
  }
  # Anything else (e.g. block = 1:5): hand it back unchanged so evzinb() /
  # run_evzinb() report it through the usual "single string" validation.
  val
}
