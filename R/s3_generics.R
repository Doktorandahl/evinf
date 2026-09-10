# Standard S3 generics for evzinb / evinb objects (audit 4.1).
#
# coef(), vcov(), confint(), logLik(), AIC()/BIC() (via the logLik method),
# nobs(), formula(), terms(), model.frame(), fitted(), residuals(), simulate(),
# update(). These plug the models into modelsummary / marginaleffects / texreg /
# sandwich, which dispatch on exactly these.

#' @importFrom stats coef vcov confint nobs fitted residuals simulate update
NULL

# --- shared helpers ---------------------------------------------------------

# Component -> $coef slot name.
.evinf_coef_slots <- c(count = "Beta.NB", zero = "Beta.multinom.ZC",
                       evi = "Beta.multinom.PL", pareto = "Beta.PL")

# Flatten one $coef list into a named vector: <component>_<term>, then alpha_nb,
# c_ev. `coefs` is object$coef (full model) or a bootstrap's $coef.
evinf_flatten_coef <- function(coefs, components) {
  pieces <- lapply(components, function(cn) {
    v <- coefs[[.evinf_coef_slots[[cn]]]]
    if (is.null(v)) {
      return(NULL)
    }
    stats::setNames(as.numeric(v), paste0(cn, "_", names(v)))
  })
  flat <- unlist(pieces, use.names = TRUE)
  c(flat,
    alpha_nb = as.numeric(coefs$Alpha.NB),
    c_ev = as.numeric(coefs$C))
}

# The canonical component list for an object (count/zero/evi/pareto, no zero for
# evinb).
evinf_components <- function(object) {
  if (inherits(object, "evzinb")) {
    c("count", "zero", "evi", "pareto")
  } else {
    c("count", "evi", "pareto")
  }
}

# n x p matrix of bootstrap coefficient draws, columns in coef() order/names.
evinf_boot_coef_matrix <- function(object, exclude_degenerate = TRUE) {
  if (is.null(object$bootstraps)) {
    stop("This requires a model fitted with bootstrap = TRUE.", call. = FALSE)
  }
  boots <- evinf_usable_bootstraps(object, exclude_degenerate)
  if (!length(boots)) {
    stop("No usable bootstrap fits (all errored or were degenerate); cannot ",
         "compute a bootstrap covariance. See failed_bootstraps().",
         call. = FALSE)
  }
  comps <- evinf_components(object)
  m <- do.call(rbind, lapply(boots, function(b) evinf_flatten_coef(b$coef, comps)))
  colnames(m) <- names(evinf_flatten_coef(object$coef, comps))
  m
}

evinf_nobs <- function(object) nrow(object$data$x.nb)


# --- coef / vcov / confint -------------------------------------------------

#' Coefficients of an evzinb / evinb model
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model.
#' @param component One of \code{"all"} (the default), \code{"count"},
#'   \code{"zero"}, \code{"evi"}, \code{"pareto"}. Deprecated aliases \code{"nb"},
#'   \code{"zi"}, \code{"evinf"} are accepted.
#' @param ... Unused.
#' @return For \code{"all"}, a named numeric vector \code{<component>_<term>},
#'   \code{alpha_nb}, \code{c_ev}; for a single component the plain named vector.
#' @export
coef.evzinb <- function(object, component = "all", ...) {
  component <- normalize_component(component,
    c("all", "count", "zero", "evi", "pareto"))
  if (component == "all") {
    return(evinf_flatten_coef(object$coef, evinf_components(object)))
  }
  v <- object$coef[[.evinf_coef_slots[[component]]]]
  if (is.null(v)) {
    stop("Component ", sQuote(component), " is not present in this model.",
         call. = FALSE)
  }
  stats::setNames(as.numeric(v), names(v))
}

#' @rdname coef.evzinb
#' @export
coef.evinb <- function(object, component = "all", ...) {
  component <- normalize_component(component, c("all", "count", "evi", "pareto"))
  coef.evzinb(object, component = component, ...)
}

#' Bootstrap covariance matrix of an evzinb / evinb model
#'
#' @param object A fitted model with bootstraps.
#' @param exclude_degenerate Drop bootstrap replicates flagged degenerate
#'   (default \code{TRUE}); see \code{\link{evinf_control}}.
#' @param ... Unused.
#' @return The covariance matrix of the bootstrap coefficient draws, with the
#'   same names and order as \code{coef(object)} (includes \code{alpha_nb} and
#'   \code{c_ev}).
#' @details \code{c_ev} is estimated on the grid of unique observed values of
#'   the response, not by a smooth optimiser, so it sits on a discrete scale.
#'   Its bootstrap variance (and a percentile \code{confint()} on it) is
#'   meaningful, but delta-method standard errors that perturb \code{c_ev}
#'   continuously -- as \code{marginaleffects::avg_slopes()} etc. do -- should
#'   be read with that in mind. For uncertainty on covariate \emph{effects},
#'   prefer the bootstrap route, \code{\link{marginal_effects}()}.
#' @export
vcov.evzinb <- function(object, exclude_degenerate = TRUE, ...) {
  stats::cov(evinf_boot_coef_matrix(object, exclude_degenerate))
}

#' @rdname vcov.evzinb
#' @export
vcov.evinb <- function(object, exclude_degenerate = TRUE, ...) {
  stats::cov(evinf_boot_coef_matrix(object, exclude_degenerate))
}

#' Confidence intervals for an evzinb / evinb model
#'
#' @param object A fitted model with bootstraps.
#' @param parm Which parameters (names as in \code{coef(object)}); default all.
#' @param level Confidence level.
#' @param type \code{"percentile"} (bootstrap percentile intervals, the default)
#'   or \code{"approx"} (\code{estimate +/- qnorm() * sqrt(diag(vcov))}).
#' @param exclude_degenerate Drop bootstrap replicates flagged degenerate
#'   (default \code{TRUE}); see \code{\link{evinf_control}}.
#' @param ... Unused.
#' @return A two-column matrix.
#' @examples
#' \donttest{
#' data(genevzinb2)
#' m <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 25)
#' # the Pareto shape (alpha_nb) and the threshold (c_ev) are in coef()/confint()
#' confint(m, parm = c("count_x1", "alpha_nb", "c_ev"))
#' }
#' @export
confint.evzinb <- function(object, parm, level = 0.95,
                           type = c("percentile", "approx"),
                           exclude_degenerate = TRUE, ...) {
  type <- match.arg(type)
  est <- coef(object, "all")
  m <- evinf_boot_coef_matrix(object, exclude_degenerate)
  if (missing(parm) || is.null(parm)) {
    parm <- names(est)
  }
  a <- (1 - level) / 2
  labs <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")

  if (type == "percentile") {
    ci <- t(apply(m[, parm, drop = FALSE], 2L, stats::quantile,
                  probs = c(a, 1 - a), names = FALSE))
  } else {
    se <- sqrt(diag(stats::cov(m)))[parm]
    z <- stats::qnorm(1 - a)
    ci <- cbind(est[parm] - z * se, est[parm] + z * se)
  }
  colnames(ci) <- labs
  rownames(ci) <- parm
  ci
}

#' @rdname confint.evzinb
#' @export
confint.evinb <- function(object, parm, level = 0.95,
                          type = c("percentile", "approx"),
                          exclude_degenerate = TRUE, ...) {
  confint.evzinb(object, parm = parm, level = level, type = match.arg(type),
                 exclude_degenerate = exclude_degenerate, ...)
}


# --- logLik / nobs / formula / terms / model.frame ------------------------

#' logLik / nobs / model.frame / formula / terms for evzinb / evinb models
#'
#' @param object A fitted model.
#' @param formula A fitted model (for the \code{formula} method).
#' @param component Which component (\code{"count"} by default) for
#'   \code{formula()} / \code{terms()}.
#' @param x A fitted model (for the \code{terms} method).
#' @param ... Unused.
#' @return As for the corresponding \pkg{stats} generic. \code{logLik()} carries
#'   \code{df = length(object$par.all)} and \code{nobs}, so \code{AIC()} /
#'   \code{BIC()} work and match \code{object$AIC} / \code{object$BIC}.
#' @name evinf-s3-accessors
#' @export
logLik.evzinb <- function(object, ...) {
  structure(object$log.lik,
            df = length(object$par.all),
            nobs = evinf_nobs(object),
            class = "logLik")
}
#' @rdname evinf-s3-accessors
#' @export
logLik.evinb <- logLik.evzinb

#' @rdname evinf-s3-accessors
#' @export
nobs.evzinb <- function(object, ...) evinf_nobs(object)
#' @rdname evinf-s3-accessors
#' @export
nobs.evinb <- function(object, ...) evinf_nobs(object)

#' @rdname evinf-s3-accessors
#' @export
model.frame.evzinb <- function(formula, ...) formula$data$data
#' @rdname evinf-s3-accessors
#' @export
model.frame.evinb <- function(formula, ...) formula$data$data

.evinf_formula_slot <- c(count = "formula_nb", zero = "formula_zi",
                         evi = "formula_evi", pareto = "formula_pareto")

#' @rdname evinf-s3-accessors
#' @export
formula.evzinb <- function(x, component = "count", ...) {
  component <- normalize_component(component, c("count", "zero", "evi", "pareto"))
  x$formulas[[.evinf_formula_slot[[component]]]]
}
#' @rdname evinf-s3-accessors
#' @export
formula.evinb <- function(x, component = "count", ...) {
  component <- normalize_component(component, c("count", "evi", "pareto"))
  x$formulas[[.evinf_formula_slot[[component]]]]
}

#' @rdname evinf-s3-accessors
#' @export
terms.evzinb <- function(x, component = "count", ...) {
  component <- normalize_component(component, c("count", "zero", "evi", "pareto"))
  x$terms[[switch(component, count = "nb", zero = "zi", evi = "evi",
                  pareto = "pareto")]]
}
#' @rdname evinf-s3-accessors
#' @export
terms.evinb <- function(x, component = "count", ...) {
  component <- normalize_component(component, c("count", "evi", "pareto"))
  x$terms[[switch(component, count = "nb", evi = "evi", pareto = "pareto")]]
}


# --- fitted / residuals / simulate ---------------------------------------

#' Fitted values, residuals and simulations from an evzinb / evinb model
#'
#' @param object A fitted model.
#' @param type For \code{fitted()}, one of \code{"harmonic"}, \code{"explog"},
#'   \code{"counts"}, \code{"pareto_alpha"}. For \code{residuals()},
#'   \code{"response"} (\eqn{y - } harmonic prediction) or \code{"quantile"}
#'   (randomized quantile residuals from the mixture CDF).
#' @param seed Optional RNG seed for the randomized quantile residuals /
#'   \code{simulate()}.
#' @param nsim Number of simulated response vectors.
#' @param newdata Optional data to simulate for.
#' @param ... Unused.
#' @return \code{fitted()} / \code{residuals()} return a numeric vector;
#'   \code{simulate()} a data frame with columns \code{sim_1}, \code{sim_2}, ...
#' @name evinf-s3-predict
#' @export
fitted.evzinb <- function(object,
                          type = c("harmonic", "explog", "counts", "pareto_alpha"),
                          ...) {
  type <- match.arg(type)
  stats::predict(object, type = type)
}
#' @rdname evinf-s3-predict
#' @export
fitted.evinb <- fitted.evzinb

#' @rdname evinf-s3-predict
#' @export
residuals.evzinb <- function(object, type = c("response", "quantile"),
                             seed = NULL, ...) {
  type <- match.arg(type)
  y <- object$data$y
  if (type == "response") {
    return(y - stats::predict(object, type = "harmonic"))
  }
  probs <- object$props
  if (ncol(probs) == 2L) {          # evinb: prepend a zero-state column
    probs <- cbind(0, probs)
  }
  mu <- object$fitted$mu.nb
  alph <- object$fitted$alpha.pl
  Fy <- mixture_p(y, alph, object$coef$C, mu, object$coef$Alpha.NB, probs)
  Fy1 <- ifelse(y <= 0, 0,
                mixture_p(y - 1, alph, object$coef$C, mu, object$coef$Alpha.NB, probs))
  Fy <- pmin(pmax(Fy, 0), 1)
  Fy1 <- pmin(pmax(Fy1, 0), Fy)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  u <- stats::runif(length(y), Fy1, Fy)
  stats::qnorm(u)
}
#' @rdname evinf-s3-predict
#' @export
residuals.evinb <- residuals.evzinb

#' @rdname evinf-s3-predict
#' @export
simulate.evzinb <- function(object, nsim = 1, seed = NULL, newdata = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  draws <- revzinb_fit(object, newdata = newdata, n_draws = nsim)
  if (nsim == 1L) {
    draws <- list(draws)
  }
  out <- as.data.frame(draws)
  names(out) <- paste0("sim_", seq_len(nsim))
  out
}
#' @rdname evinf-s3-predict
#' @export
simulate.evinb <- function(object, nsim = 1, seed = NULL, newdata = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  draws <- revinb_fit(object, newdata = newdata, n_draws = nsim)
  if (nsim == 1L) {
    draws <- list(draws)
  }
  out <- as.data.frame(draws)
  names(out) <- paste0("sim_", seq_len(nsim))
  out
}


# --- update --------------------------------------------------------------

#' Update and re-fit an evzinb / evinb model
#'
#' @param object A fitted model (must carry \code{object$call}).
#' @param formula_nb.,formula_zi.,formula_evi.,formula_pareto. Optional
#'   \code{\link[stats]{update.formula}}-style changes per component, e.g.
#'   \code{formula_pareto. = . ~ . - x3}.
#' @param ... Other arguments of \code{\link{evzinb}} / \code{\link{evinb}} to
#'   change.
#' @param evaluate If \code{TRUE} (default) re-fit; otherwise return the updated call.
#' @return The updated fit, or the call.
#' @export
update.evzinb <- function(object, formula_nb., formula_zi., formula_evi.,
                          formula_pareto., ..., evaluate = TRUE) {
  evinf_update_impl(object,
                    if (!missing(formula_nb.)) formula_nb.,
                    if (!missing(formula_zi.)) formula_zi.,
                    if (!missing(formula_evi.)) formula_evi.,
                    if (!missing(formula_pareto.)) formula_pareto.,
                    list(...), evaluate, parent.frame())
}

#' @rdname update.evzinb
#' @export
update.evinb <- function(object, formula_nb., formula_evi., formula_pareto.,
                         ..., evaluate = TRUE) {
  evinf_update_impl(object,
                    if (!missing(formula_nb.)) formula_nb.,
                    NULL,
                    if (!missing(formula_evi.)) formula_evi.,
                    if (!missing(formula_pareto.)) formula_pareto.,
                    list(...), evaluate, parent.frame())
}

evinf_update_impl <- function(object, f_nb, f_zi, f_evi, f_pareto,
                              extras, evaluate, env) {
  cl <- object$call
  if (is.null(cl)) {
    stop("This model carries no stored call; refit with evinf >= 0.10.0.",
         call. = FALSE)
  }
  is_zinb <- inherits(object, "evzinb")
  cl[[1L]] <- if (is_zinb) quote(evinf::evzinb) else quote(evinf::evinb)

  # data: keep the user's original expression and evaluate the call in a child
  # of `env`; fall back to the embedded data frame (bound to `.evinf_data`) only
  # when that expression can no longer be resolved (audit N6). Either way the
  # refitted object's $call$data is an expression, never a data frame.
  eval_env <- new.env(parent = env)
  eval_env$.evinf_data <- object$data$data
  data_expr <- cl$data
  can_eval <- !is.null(data_expr) &&
    is.data.frame(tryCatch(eval(data_expr, env), error = function(e) NULL))
  if (!can_eval) {
    cl$data <- quote(.evinf_data)
  }

  cl$control <- object$control
  cl$block <- object$block

  chg <- function(component, change) {
    old <- object$formulas[[.evinf_formula_slot[[component]]]]
    if (is.null(change)) old else stats::update.formula(old, change)
  }
  cl$formula_nb <- chg("count", f_nb)
  if (is_zinb) {
    cl$formula_zi <- chg("zero", f_zi)
  }
  cl$formula_evi <- chg("evi", f_evi)
  cl$formula_pareto <- chg("pareto", f_pareto)

  for (nm in names(extras)) {
    cl[[nm]] <- extras[[nm]]
  }
  if (evaluate) eval(cl, eval_env) else cl
}
