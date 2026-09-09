.onLoad <- function(libname, pkgname) {
  # Register the model classes with marginaleffects (audit 4.8) so
  # marginaleffects::avg_slopes(model) etc. dispatch on the get_*/set_coef
  # methods this package provides.
  cur <- getOption("marginaleffects_model_classes", default = character())
  add <- setdiff(c("evzinb", "evinb"), cur)
  if (length(add)) {
    options(marginaleffects_model_classes = c(cur, add))
  }
  invisible()
}
