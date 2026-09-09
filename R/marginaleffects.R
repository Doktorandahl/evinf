# audit 4.8 - marginaleffects compatibility.
#
# marginaleffects dispatches on get_coef / set_coef / get_vcov / get_predict for
# a custom class. With these, marginaleffects::avg_slopes(model) etc. work with
# delta-method standard errors from the bootstrap vcov().

evinf_set_flat_coef <- function(model, coefs) {
  comps <- evinf_components(model)
  for (cn in comps) {
    slot <- .evinf_coef_slots[[cn]]
    nm <- names(model$coef[[slot]])
    model$coef[[slot]] <- stats::setNames(unname(coefs[paste0(cn, "_", nm)]), nm)
  }
  model$coef$Alpha.NB <- unname(coefs[["alpha_nb"]])
  model$coef$C <- unname(coefs[["c_ev"]])
  model
}

#' @exportS3Method marginaleffects::get_coef evzinb
get_coef.evzinb <- function(model, ...) coef(model, "all")
#' @exportS3Method marginaleffects::get_coef evinb
get_coef.evinb <- function(model, ...) coef(model, "all")

#' @exportS3Method marginaleffects::set_coef evzinb
set_coef.evzinb <- function(model, coefs, ...) evinf_set_flat_coef(model, coefs)
#' @exportS3Method marginaleffects::set_coef evinb
set_coef.evinb <- function(model, coefs, ...) evinf_set_flat_coef(model, coefs)

#' @exportS3Method marginaleffects::get_vcov evzinb
get_vcov.evzinb <- function(model, ...) vcov(model)
#' @exportS3Method marginaleffects::get_vcov evinb
get_vcov.evinb <- function(model, ...) vcov(model)

evinf_get_predict <- function(model, newdata, type) {
  ok <- c("harmonic", "explog", "counts", "pareto_alpha")
  if (length(type) != 1L || is.na(type) || !type %in% ok) {
    type <- "harmonic"
  }
  p <- stats::predict(model, newdata = newdata, type = type)
  data.frame(rowid = seq_len(NROW(p)), estimate = as.numeric(p))
}

#' @exportS3Method marginaleffects::get_predict evzinb
get_predict.evzinb <- function(model, newdata, type = "harmonic", ...) {
  evinf_get_predict(model, newdata, type)
}
#' @exportS3Method marginaleffects::get_predict evinb
get_predict.evinb <- function(model, newdata, type = "harmonic", ...) {
  evinf_get_predict(model, newdata, type)
}

# insight (used by marginaleffects) tries to recover the model frame from the
# fitting environment first and warns when it cannot. The full model frame,
# with the response and every raw covariate, is stored on the object, so hand
# it straight back (issue 3.3).
#' @exportS3Method insight::get_data evzinb
get_data.evzinb <- function(x, ...) x$data$data
#' @exportS3Method insight::get_data evinb
get_data.evinb <- function(x, ...) x$data$data
