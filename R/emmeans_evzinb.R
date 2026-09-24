# round10 I.3 (audit sec5.9, decision D6): emmeans support for the count
# component's linear predictor (log link), via delayed S3 registration --
# recover_data.evzinb()/emm_basis.evzinb() are plain (unexported) functions
# here, registered with emmeans by emmeans::.emm_register() in .onLoad()
# (R/zzz.R) rather than via NAMESPACE S3method() entries, exactly the
# pattern emmeans' own "extending emmeans" vignette recommends for a package
# that only Suggests emmeans.
#
# Not implemented: a multinomial-logit basis over the three (zero/count/evi)
# state probabilities, the other half of the round10 brief for this item.
# emmeans' own multinomial support (emm_basis.multinom(), for nnet::multinom)
# relies on an unexported "clr"-transform post-grid hook
# (emmeans:::.multinom.postGrid) specific to nnet::multinom's k x (p+1)
# coefficient matrix -- evzinb's three separately-formula'd components
# (Beta.multinom.ZC / Beta.multinom.PL, each its own formula and design
# matrix) don't map onto that representation, and reverse-engineering an
# unexported, version-fragile internal to fabricate one risks a silently
# wrong delta-method Jacobian for a stats package's own confidence
# intervals -- worse than not shipping it (this is the round10 brief's own
# explicit fallback: "if the multinomial basis turns out to be
# disproportionately hard, ship the count component only and report it").
# predict(type = "states", confint = TRUE) already gives bootstrap-based
# state-probability intervals; use that instead.

# `component` is accepted (not just silently ignored) so a caller who tries
# component = "states"/"zero"/"evi" gets a clear pointer to the fallback
# above rather than emm_basis()'s generic "'arg' should be one of" error.
evinf_emmeans_component <- function(component) {
  if (!identical(component, "count")) {
    stop(
      "emmeans support for evzinb/evinb models covers the count component ",
      "only (component = \"count\"); a state-probability basis is not ",
      "implemented (round10 I.3) -- see ?emmeans-evzinb. Use ",
      "predict(type = \"states\", confint = TRUE) for state-probability ",
      "intervals instead.", call. = FALSE
    )
  }
  "count"
}

#' emmeans support for evzinb / evinb models
#'
#' On load, if \pkg{emmeans} is installed, the package registers
#' \code{recover_data.evzinb()}/\code{emm_basis.evzinb()} (and the
#' \code{evinb} equivalents) with \pkg{emmeans} via
#' \code{emmeans::.emm_register()}, so \code{emmeans::emmeans(model, ...)}
#' works for the count component's linear predictor (a log link;
#' \code{type = "response"}/\code{regrid()} back-transforms it the usual
#' \pkg{emmeans} way). Standard errors come from the count-component block of
#' the bootstrap covariance matrix (\code{\link{vcov.evzinb}}).
#'
#' @details Only the count component is supported (\code{component =
#'   "count"}, the default and, for now, the only accepted value) -- see the
#'   source comments in \code{R/emmeans_evzinb.R} for why a
#'   state-probability (zero/count/evi) multinomial basis is not
#'   implemented. \code{predict(type = "states", confint = TRUE)} covers
#'   that case instead, with bootstrap rather than delta-method intervals.
#'
#' @name emmeans-evzinb
#' @keywords internal
NULL

recover_data.evzinb <- function(object, component = "count", ...) {
  # round10 I.3: the component = "not count" error is raised from
  # emm_basis() only, not here -- ref_grid() calls recover_data() inside
  # try(), which discards this function's own stop() message in favour of a
  # generic "Perhaps a 'data' or 'params' argument is needed"; emm_basis()
  # is called unwrapped, so an error there reaches the user with its actual
  # message intact. `component` is still a formal argument here so it is
  # accepted (not "unused argument") when the reference grid forwards it.
  emmeans::recover_data(object$call, object$terms$nb, object$na.action, ...)
}

recover_data.evinb <- recover_data.evzinb

emm_basis.evzinb <- function(object, trms, xlev, grid, component = "count", ...) {
  component <- evinf_emmeans_component(component)

  m <- stats::model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X <- stats::model.matrix(trms, m)

  bhat <- coef(object, component = "count")
  X <- X[, names(bhat), drop = FALSE]
  bhat <- as.numeric(bhat)

  nm <- paste0("count_", names(coef(object, component = "count")))
  V <- vcov(object)[nm, nm, drop = FALSE]
  dimnames(V) <- NULL

  fam <- object$family %||% evinf_family()
  misc <- if (identical(fam$count, "poisson")) {
    list(tran = "log", inv.lbl = "rate")
  } else {
    list(tran = "log", inv.lbl = "count")
  }

  list(
    X = X, bhat = bhat, nbasis = estimability::all.estble, V = V,
    dffun = function(k, dfargs) Inf, dfargs = list(), misc = misc
  )
}

emm_basis.evinb <- emm_basis.evzinb
