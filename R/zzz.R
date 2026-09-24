# nocov start
.onLoad <- function(libname, pkgname) {
  # Register the model classes with marginaleffects (audit 4.8) so
  # marginaleffects::avg_slopes(model) etc. dispatch on the get_*/set_coef
  # methods this package provides.
  cur <- getOption("marginaleffects_model_classes", default = character())
  add <- setdiff(c("evzinb", "evinb"), cur)
  if (length(add)) {
    options(marginaleffects_model_classes = c(cur, add))
  }

  # round10 I.3 (audit §5.9): delayed S3 registration of
  # recover_data.evzinb()/emm_basis.evzinb() (R/emmeans_evzinb.R) with
  # emmeans, the pattern emmeans' own "extending emmeans" vignette
  # recommends for a package that only Suggests emmeans.
  if (requireNamespace("emmeans", quietly = TRUE)) {
    emmeans::.emm_register(c("evzinb", "evinb"), pkgname)
  }

  invisible()
}
# nocov end
