# audit 2.7: one canonical vocabulary for model components.
#
# Canonical names: "count", "zero", "evi", "pareto" (+ "all" where relevant).
# The pre-0.9.4 names are accepted as deprecated aliases:
#   "nb"   -> "count"
#   "zi"   -> "zero"
#   "evinf"-> "evi"

.evinf_component_aliases <- c(
  count  = "count", nb    = "count",
  zero   = "zero",  zi    = "zero",
  evi    = "evi",   evinf = "evi",
  pareto = "pareto",
  all    = "all"
)

#' Normalise a component argument to the canonical vocabulary
#'
#' @param component A single string; one of the canonical names or a deprecated alias.
#' @param choices The canonical names allowed in this context.
#' @return The canonical component name.
#' @noRd
normalize_component <- function(component, choices = c("count", "zero", "evi", "pareto", "all")) {
  component <- component[1]
  if (!component %in% names(.evinf_component_aliases)) {
    stop("Unknown component ", sQuote(component), ". Use one of: ",
         paste(sQuote(choices), collapse = ", "), ".", call. = FALSE)
  }
  canonical <- unname(.evinf_component_aliases[[component]])
  if (!identical(component, canonical)) {
    warning("`component = \"", component, "\"` is deprecated; use \"", canonical,
            "\" instead.", call. = FALSE)
  }
  if (!canonical %in% choices) {
    stop("Component ", sQuote(canonical), " is not available here. Use one of: ",
         paste(sQuote(choices), collapse = ", "), ".", call. = FALSE)
  }
  canonical
}

# predict(type=) keeps its own vocabulary but accepts "zero" for "zi" and
# "evi" for "evinf".
normalize_predict_type <- function(type, choices) {
  type <- type[1]
  type <- switch(type, zero = "zi", evi = "evinf", type)
  match.arg(type, choices)
}
