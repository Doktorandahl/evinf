#' Running an extreme value and zero inflated negative binomial model
#'
#' @param formula_nb,formula_zi,formula_evi,formula_pareto Component formulas.
#' @param data Data to run the model on.
#' @param control An \code{evinf_control()} object.
#' @param block Optional string naming a case-identifier column for block bootstrapping; included in the na.omit() so the returned block vector aligns with the model data.
#' @param weights Optional string naming a weight column (round9 D.2); included in the na.omit() so the returned weight vector aligns with the model data, same as block.
#' @param time Optional string naming a time-index column (round9 F); required
#'   for \code{bootstrap_scheme \%in\% c("moving_block", "stationary")}.
#' @param bootstrap_scheme One of \code{"iid"}, \code{"cluster"},
#'   \code{"moving_block"}, \code{"stationary"} (round9 F), or \code{NULL} to
#'   default to \code{"cluster"} when \code{block} is given, \code{"iid"}
#'   otherwise.
#' @param block_length Block length for \code{"moving_block"} /
#'   \code{"stationary"} (round9 F); \code{NULL} uses \code{ceiling(T^(1/3))}
#'   per unit.
#' @param family An \code{\link{evinf_family}()} object (round9 E.0/E.1).
#' @param verbose Should progress be printed for the first run of evzinb.
#'
#' @return An object of class 'evzinb'
#' @noRd
run_evzinb <- function(
  formula_nb,
  formula_zi = NULL,
  formula_evi = NULL,
  formula_pareto = NULL,
  data,
  control = evinf_control(),
  block = NULL,
  weights = NULL,
  time = NULL,
  bootstrap_scheme = NULL,
  block_length = NULL,
  family = evinf_family(),
  start_seed = NULL,
  multicore = NULL,
  ncores = NULL,
  verbose = TRUE
) {
  control <- validate_evinf_control(control)
  family <- evinf_resolve_family(family)
  if (!is.null(block) && !(is.character(block) && length(block) == 1L)) {
    stop("`block` must be NULL or a single string naming a column of `data`.",
         call. = FALSE)
  }
  if (!is.null(time) && !(is.character(time) && length(time) == 1L)) {
    stop("`time` must be NULL or a single string naming a column of `data`.",
         call. = FALSE)
  }
  bootstrap_scheme <- evinf_resolve_bootstrap_scheme(bootstrap_scheme, block, time)
  # round9 D.1 (audit §5.6): offset() is supported in the count, zero-inflation
  # and EVI-inflation components; not in the Pareto shape component (a
  # non-NB formula the user supplied with offset() there is an error). A
  # component that merely inherited formula_nb's offset by defaulting to it
  # has that offset stripped -- no inheritance, an offset applies only where
  # it is written explicitly.
  formula_evi <- evinf_component_formula(formula_evi, formula_nb, "formula_evi")
  formula_pareto <- evinf_component_formula(formula_pareto, formula_nb, "formula_pareto",
                                            allow_offset = FALSE)
  formula_zi <- evinf_component_formula(formula_zi, formula_nb, "formula_zi")

  # Restrict to the union of variables used by any component (plus the block
  # variable, audit R0.1) and drop incomplete rows once, so the design matrices
  # and the block vector all refer to the same rows (audit 1.7, R0.1).
  model_vars <- unique(c(
    all.vars(formula_nb),
    all.vars(formula_zi),
    all.vars(formula_evi),
    all.vars(formula_pareto),
    block,
    weights,
    time
  ))
  model_data <- data %>%
    dplyr::select(dplyr::all_of(model_vars)) %>%
    na.omit()

  # round10 0.4 (review §3): validate `time` against `block` once, here,
  # before any EM fitting -- see evinf_validate_time()'s comment in
  # R/resample.R.
  evinf_validate_time(
    nrow(model_data), bootstrap_scheme,
    if (!is.null(block)) model_data[[block]] else NULL,
    if (!is.null(time)) model_data[[time]] else NULL
  )

  d_nb <- evinf_design(formula_nb, model_data)
  d_zi <- evinf_design(formula_zi, model_data)
  d_evi <- evinf_design(formula_evi, model_data)
  d_pareto <- evinf_design(formula_pareto, model_data)
  offset_nb <- d_nb$offset
  offset_zc <- d_zi$offset
  offset_pl_mult <- d_evi$offset
  # round9 D.2: sum(weights_vec) below is the effective sample size used for
  # AIC/BIC/nobs() -- rep(1, n) here reproduces today's row count exactly.
  weights_vec <- if (!is.null(weights)) model_data[[weights]] else rep(1, nrow(model_data))
  evinf_check_weights(weights_vec, nrow(model_data))

  OBS.Y <- as.matrix(model.response(model.frame(formula_nb, model_data)))
  evinf_check_response(as.numeric(OBS.Y))

  # Resolve NULL c.lim / init.C against the data (audit 4.5). This is the single
  # resolution point; lr_test() refits pass a concrete c.lim so no message fires.
  control <- evinf_resolve_c(control, as.numeric(OBS.Y))
  Control <- control

  OBS.X.obj <- list()
  OBS.X.obj$X.multinom.ZC <- d_zi$X
  OBS.X.obj$X.multinom.PL <- d_evi$X
  OBS.X.obj$X.NB <- d_nb$X
  OBS.X.obj$X.PL <- d_pareto$X
  OBS.X.obj$offset.nb <- offset_nb
  OBS.X.obj$offset.zc <- offset_zc
  OBS.X.obj$offset.pl_mult <- offset_pl_mult
  OBS.X.obj$weights <- weights_vec
  # round10 0.5: `weights` here is still the R-level argument (a column name
  # or NULL), not weights_vec -- has_weights is TRUE only when the user
  # actually supplied weights, never merely because weights_vec is all 1s.
  OBS.X.obj$has_weights <- !is.null(weights)

  init.Beta.multinom.ZC <- control$init.Beta.multinom.ZC
  init.Beta.multinom.PL <- control$init.Beta.multinom.PL
  init.Beta.NB <- control$init.Beta.NB
  init.Beta.PL <- control$init.Beta.PL

  # Parameter counts include the intercept the C++ routines prepend.
  n_zc <- ncol(OBS.X.obj$X.multinom.ZC) + 1L
  n_pl_mult <- ncol(OBS.X.obj$X.multinom.PL) + 1L
  n_nb <- ncol(OBS.X.obj$X.NB) + 1L
  n_pl <- ncol(OBS.X.obj$X.PL) + 1L

  Ini.Val <- list()
  if (is.null(init.Beta.multinom.ZC) | length(init.Beta.multinom.ZC) != n_zc) {
    Ini.Val$Beta.multinom.ZC <- rep(0, n_zc)
  } else {
    Ini.Val$Beta.multinom.ZC <- init.Beta.multinom.ZC
  }

  if (is.null(init.Beta.multinom.PL) | length(init.Beta.multinom.PL) != n_pl_mult) {
    Ini.Val$Beta.multinom.PL <- rep(0, n_pl_mult)
  } else {
    Ini.Val$Beta.multinom.PL <- init.Beta.multinom.PL
  }

  if (is.null(init.Beta.NB) | length(init.Beta.NB) != n_nb) {
    Ini.Val$Beta.NB <- rep(0, n_nb)
  } else {
    Ini.Val$Beta.NB <- init.Beta.NB
  }

  if (is.null(init.Beta.PL) | length(init.Beta.PL) != n_pl) {
    Ini.Val$Beta.PL <- rep(0, n_pl)
  } else {
    Ini.Val$Beta.PL <- init.Beta.PL
  }

  Ini.Val$Alpha.NB <- control$init.Alpha.NB
  Ini.Val$C <- control$init.C

  # round10 G.1 (audit §5.10): n_starts <= 1 (the default) is this same
  # single em_fit() call, bit-identical to before evinf_run_starts() existed.
  # The starts pmap respects the same multicore/ncores as the bootstrap
  # dispatch below (evinf_with_plan(multicore = NULL, ...) is a no-op, so
  # this defaults to whatever future::plan() is already active, same as
  # every other parallel loop in the package).
  starts_result <- evinf_with_plan(multicore, ncores, {
    evinf_run_starts(OBS.Y, OBS.X.obj, Ini.Val, Control,
                     model = "evzinb", family = family,
                     verbose = verbose, start_seed = start_seed)
  })
  object <- starts_result$best
  if (verbose && isTRUE(object$loglik_recomputed)) {
    message(
      "The trace-maximum log-likelihood differed from the value at the ",
      "returned parameters; log.lik / AIC / BIC were recomputed from the ",
      "returned parameters."
    )
  }
  object$par.mat$Beta.multinom.ZC <- as.numeric(object$par.mat$Beta.multinom.ZC)
  object$par.mat$Beta.multinom.PL <- as.numeric(object$par.mat$Beta.multinom.PL)
  object$par.mat$Beta.NB <- as.numeric(object$par.mat$Beta.NB)
  object$par.mat$Beta.PL <- as.numeric(object$par.mat$Beta.PL)

  names(object$par.mat$Beta.NB) <- c('(Intercept)', colnames(d_nb$X))
  names(object$par.mat$Beta.multinom.ZC) <- c('(Intercept)', colnames(d_zi$X))
  names(object$par.mat$Beta.multinom.PL) <- c('(Intercept)', colnames(d_evi$X))
  names(object$par.mat$Beta.PL) <- c('(Intercept)', colnames(d_pareto$X))

  # round9 E.1: Alpha.NB is not a parameter of a Poisson count state -- drop
  # it (absent, not NA) so coef()/vcov()/confint()/tidy()/glance() all show
  # it as absent (evinf_flatten_coef(), glance.evzinb()).
  if (family$count == "poisson") {
    object$par.mat$Alpha.NB <- NULL
  }

  object$formulas <- list(
    formula_nb = formula_nb,
    formula_zi = formula_zi,
    formula_evi = formula_evi,
    formula_pareto = formula_pareto
  )
  object$terms <- list(
    nb = d_nb$terms,
    zi = d_zi$terms,
    evi = d_evi$terms,
    pareto = d_pareto$terms
  )
  object$xlevels <- list(
    nb = d_nb$xlevels,
    zi = d_zi$xlevels,
    evi = d_evi$xlevels,
    pareto = d_pareto$xlevels
  )
  object$control <- control
  object$c_lim_default <- isTRUE(control$c_lim_default)
  object$has_offset <- isTRUE(d_nb$has_offset)
  object$offset_nb <- offset_nb
  object$offset_zc <- offset_zc
  object$offset_pl_mult <- offset_pl_mult
  object$weights <- weights_vec
  object$has_weights <- !is.null(weights)
  object$data <- list()

  object$data$data <- model_data
  object$block <- block
  object$block_vec <- if (!is.null(block)) model_data[[block]] else NULL
  object$time <- time
  object$time_vec <- if (!is.null(time)) model_data[[time]] else NULL
  object$bootstrap_scheme <- bootstrap_scheme
  # round9 F: resolve block_length = NULL to a concrete per-unit vector once,
  # here, rather than once per bootstrap replicate -- evinf_resample_ids()'s
  # default-block-length message would otherwise fire (and recompute
  # identical values) on every one of potentially hundreds of replicates.
  object$block_length <- if (bootstrap_scheme %in% c("moving_block", "stationary")) {
    evinf_resolve_block_length(nrow(model_data), object$block_vec, block_length)
  } else {
    block_length
  }
  object$data$y <- as.numeric(object$y)
  object$y <- NULL
  object$data$x.nb <- object$x.nb
  object$x.nb <- NULL
  object$data$x.pl <- object$x.pl
  object$x.pl <- NULL
  object$data$x.multinom.zc <- object$x.multinom.zc
  object$x.multinom.zc <- NULL
  object$data$x.multinom.pl <- object$x.multinom.pl
  object$x.multinom.pl <- NULL

  object$props <- object$par.mat$Props
  object$par.mat$Props <- NULL
  object$coef <- object$par.mat
  object$par.mat <- NULL
  colnames(object$props) <- colnames(object$resp) <- c(
    "zero",
    'count',
    'evi'
  )

  # audit 4.5: threshold diagnostics
  object$c_profile <- tibble::as_tibble(object$c_profile)
  object$loglik_trace <- object$log.lik.vec.all
  object$n_above_c <- sum(object$data$y >= object$coef$C)
  object$n_em_steps <- length(object$log.lik.vec.all)
  # round10 G.2: object$n_c_iter_warmup / object$n_loglik_warmup (from
  # em_fit()'s return, already merged into object above) mark where warm-up
  # ends in c_trace / loglik_trace, for plot(type = "trace").

  # round10 G.1: $starts / $start_seed are NULL for the (default) single
  # start, to keep those objects the same size as before this feature existed.
  object$starts <- starts_result$starts
  object$start_seed <- starts_result$start_seed

  object$fitted <- list()
  # round10 0.6 (review §4/§7, breaking change): the legacy y.hat.pl_* point
  # predictions duplicated predict() (harmonic/explog/quantile) but were not
  # hurdle-aware (they mix with the raw prior state probabilities, not the
  # hurdle-adjusted ones) -- cor() with the corresponding predict() output
  # was as low as 0.9988 on a hurdle fit. No internal reader of
  # object$fitted$y.hat.pl_* remained (grepped before removing); fitted()
  # already routes through predict() rather than these fields.
  object$y.hat.plexpElogy <- NULL
  object$y.hat.pl.E.inv.y <- NULL
  object$y.hat.plmedian <- NULL
  object$y.hat.plmean <- NULL
  object$fitted$mu.nb <- object$mu.nb.vec
  object$mu.nb.vec <- NULL
  object$fitted$alpha.pl <- object$alpha.pl.vec
  object$alpha.pl.vec <- NULL
  object$fitted$pl_exp.E.log.y <- object$exp.E.log.y
  object$exp.E.log.y <- NULL
  object$fitted$pl_median <- object$median.pl.vec
  object$median.pl.vec <- NULL
  object$fitted$pl_mean <- object$mean.pl.vec
  object$mean.pl.vec <- NULL
  object$fitted$prob_zero <- object$props[, 1]
  object$fitted$prob_count <- object$props[, 2]
  object$fitted$prob_evi <- object$props[, 3]
  object$fitted$posterior_zero <- object$resp[, 1]
  object$fitted$posterior_count <- object$resp[, 2]
  object$fitted$posterior_evi <- object$resp[, 3]
  # Deprecated pre-0.9.4 names, kept as duplicates for one release (audit 2.7).
  object$fitted$prob_pareto <- object$props[, 3]
  object$fitted$posterior_pareto <- object$resp[, 3]

  class(object) <- 'evzinb'
  return(object)
}


#' Running an extreme value and zero inflated negative binomial model with bootstrapping
#'
#' @param formula_nb Formula for the negative binomial (count) component of the model.
#'   May include an \code{offset()} term (\eqn{\mu_{NB} = \exp(x'\beta + offset)}).
#' @param formula_zi Formula for the zero-inflation component of the model. If NULL
#'   taken as the same formula as nb, with any \code{offset()} term stripped (an
#'   offset applies only where it is written explicitly, never by inheritance).
#'   May include its own \code{offset()} term: since the zero-inflation logit is
#'   the log-odds of the zero state \emph{against} the count state, an offset there
#'   shifts that log-odds, e.g. \code{offset(log(exposure))} makes a larger
#'   exposure relatively less likely to land in the structural-zero state.
#' @param formula_evi Formula for the extreme-value inflation component of the model.
#'   If NULL taken as the same formula as nb, offset stripped as above. May include
#'   its own \code{offset()} term (e.g. \code{offset(log(population))} for a
#'   probability of an extreme event that scales with exposure), read the same way:
#'   it shifts the EVI log-odds against the count state.
#' @param formula_pareto Formula for the pareto (extreme value) component of the
#'   model. If NULL taken as the same formula as nb, offset stripped as above.
#'   \code{offset()} is \strong{not} supported here (errors if present): an offset
#'   on a shape parameter has no clear reading.
#' @param data data to run the model on
#' @param bootstrap Should bootstrapping be performed. Needed to obtain standard errors and p-values
#' @param n_bootstraps Number of bootstraps to run. For use of bootstrapped p-values, at least 1,000 bootstraps are recommended. For approximate p-values, a lower number can be sufficient
#' @param multicore Bootstrap parallelisation shortcut. The default
#'   (\code{NULL}) respects whatever \code{future} plan is currently set (see
#'   the \strong{Parallel processing} section). \code{TRUE} sets a temporary
#'   \code{\link[future]{multisession}} plan for the duration of the call;
#'   \code{FALSE} forces sequential execution. The previous plan is always
#'   restored on exit.
#' @param ncores Number of workers when \code{multicore = TRUE}. Default
#'   (\code{NULL}) is one less than the number of available cores. Ignored when
#'   \code{multicore} is \code{NULL} or \code{FALSE}.
#' @param block Optional case-identifier column for block bootstrapping, given
#'   either as a bare column name (\code{block = id}) or a string
#'   (\code{block = "id"}). Note that the bundled \code{\link{hks}} data contain
#'   no conflict identifier, so the conflict-level cluster bootstrap in Randahl
#'   and Vegelius (2024) cannot be reproduced from them directly (see
#'   \code{?hks}).
#' @param weights Optional observation weights (round9 D.2), given as a bare
#'   column name (\code{weights = wt}), a string naming a column
#'   (\code{weights = "wt"}), or a numeric vector. Must be positive and
#'   finite; need not be integers (analytic weights are allowed, not just
#'   frequency counts). \strong{Frequency-weight semantics}: every
#'   observation's contribution to the log-likelihood and to the EM/M-step
#'   accumulations is multiplied by its weight, and \code{nobs()} -- and
#'   therefore \code{AIC}, \code{BIC} and the approximate t-based p-values --
#'   use \code{sum(weights)}, not the row count (see \code{sum_weights} in
#'   \code{\link{glance.evzinb}}). This interpretation of AIC/BIC assumes the
#'   weights really are frequency weights (repeat-count equivalents); for
#'   analytic weights the information-criterion values are still computed
#'   this way but their usual interpretation is weaker. \strong{\code{weights}
#'   does not give design-based standard errors for survey data} -- a
#'   sampling weight changes the point estimate, not the variance under the
#'   sampling design; for that, resample primary sampling units with
#'   \code{block =} instead. The bootstrap resamples rows exactly as without
#'   weights and carries each drawn row's weight along (the resampling
#'   probabilities themselves are not reweighted).
#' @param time Optional time index for panel/time-series bootstrap resampling
#'   (round9 F), given as a bare column name (\code{time = t}) or a string
#'   (\code{time = "t"}); required for \code{bootstrap_scheme \%in\%
#'   c("moving_block", "stationary")}. Must be strictly increasing within
#'   every \code{block} unit's rows as they already appear in the data --
#'   this is never sorted for you; sort \code{data} by \code{(block, time)}
#'   first if it isn't already.
#' @param bootstrap_scheme One of \code{"iid"} (plain row resampling, the
#'   default when \code{block} is not given), \code{"cluster"} (resample
#'   whole \code{block} units, the default when \code{block} is given -- what
#'   \code{block =} has always done), \code{"moving_block"} or
#'   \code{"stationary"} (block-resample each unit's own time series; see
#'   Kunsch 1989 / Politis and Romano 1994). The two block schemes need
#'   \code{time}; \code{block} is optional for them (the whole data is
#'   treated as one unit when omitted). With overlapping blocks the
#'   out-of-bag set is smaller and more temporally correlated than under iid
#'   resampling, so out-of-bag error (\code{\link{oob_evaluation}}) is
#'   optimistic relative to genuine forecasting performance; each bootstrap
#'   replicate's realised out-of-bag fraction is stored as
#'   \code{$oob_fraction} next to \code{$boot_id}.
#' @param block_length Block length for \code{bootstrap_scheme \%in\%
#'   c("moving_block", "stationary")}; \code{NULL} (the default) uses
#'   \code{ceiling(T^(1/3))} for each unit's own length \eqn{T} (a message
#'   names the value(s) used, once, at the original fit -- not on every
#'   bootstrap replicate).
#' @param family An \code{\link{evinf_family}()} object, or a string as
#'   shorthand for its \code{count} argument (round9 E.0/E.1/E.2), e.g.
#'   \code{family = "poisson"}. The default reproduces today's
#'   negative-binomial, mixture-zero model exactly. \code{count = "poisson"}
#'   drops \code{Alpha.NB} entirely (not merely fixes it): it is absent from
#'   \code{par.all}, \code{coef()}, \code{vcov()}, \code{confint()} and
#'   \code{tidy()}, and shown as absent (not \code{NA}) in \code{summary()}
#'   and \code{glance()}. \code{zero = "hurdle"} makes the zero state own
#'   every zero (rather than competing with the count state for them) and
#'   zero-truncates the count state; verified to match
#'   \code{pscl::hurdle()}'s coefficients and log-likelihood to numerical
#'   precision when the extreme-value state is unreachable.
#' @param boot_seed Optional bootstrap seed for reproducibility. When supplied
#'   it is used as-is; when \code{NULL} a seed is drawn and recorded, so
#'   \code{object$boot_seeds} is always populated for a bootstrapped model
#'   (see \code{\link{add_bootstraps}}).
#' @param start_seed Optional seed for the perturbed starts when
#'   \code{control$n_starts > 1} (round10 G.1); recorded as
#'   \code{object$start_seed} (\code{NULL} for the default \code{n_starts =
#'   1}, same as \code{boot_seed} for an unbootstrapped fit). Unused
#'   otherwise.
#' @param control An \code{\link{evinf_control}()} object holding the EM tuning
#'   settings (tolerances, candidate range for \eqn{C_{EV}}, BFGS steps, starting
#'   values, ...).
#' @param max.diff.par,max.no.em.steps,max.no.em.steps.warmup,c.lim,prune.c.range,max.upd.par.zc.multinomial,max.upd.par.pl.multinomial,max.upd.par.nb,max.upd.par.pl,no.m.bfgs.steps.multinomial,no.m.bfgs.steps.nb,no.m.bfgs.steps.pl,pdf.pl.type,eta.int,init.Beta.multinom.ZC,init.Beta.multinom.PL,init.Beta.NB,init.Beta.PL,init.Alpha.NB,init.C
#'   \strong{Deprecated.} These EM tuning arguments still work but should be
#'   passed through \code{control = evinf_control(...)}; supplying one directly
#'   overrides the corresponding \code{control} element and emits a warning. See
#'   \code{\link{evinf_control}} for their meaning. Will be removed in 0.11.0.
#' @param verbose Logical: should progress of the full run of the model be tracked?
#'
#' @section Parallel processing:
#' Bootstrap fits (and the per-bootstrap work in \code{\link{add_bootstraps}},
#' \code{\link{lr_test}}, \code{\link{compare_models}},
#' \code{\link{predict.evzinb}} and \code{\link{marginal_effects}}) are
#' dispatched with \pkg{furrr} on top of a \pkg{future} plan. Set the plan
#' once for your session and leave \code{multicore} at its default:
#' \preformatted{
#' future::plan(future::multisession, workers = 8)
#' progressr::handlers(global = TRUE)   # opt in to progress bars
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 1000)
#' }
#' Any \code{future} backend works (\code{multisession}, \code{cluster},
#' \code{callr}, a HPC batchtools plan, ...). As a convenience
#' \code{multicore = TRUE} sets a temporary \code{multisession} plan for the
#' single call and restores the previous plan on exit. Large models may need
#' \code{options(future.globals.maxSize = <bytes>)} to raise the default limit
#' on the data shipped to each worker.
#'
#' @section Reproducibility:
#' The bootstrap \strong{resample indices} (\code{object$bootstraps[[i]]$boot_id})
#' are fully determined by \code{boot_seed} and are independent of the
#' \pkg{future} backend, the number of workers, and chunking --- a fixed seed
#' always draws the same resamples.
#'
#' The \strong{fitted coefficients} are reproducible only to within
#' floating-point noise. BLAS operations are not bit-reproducible across
#' processes, so a sequential run and a \code{multisession} run of the same
#' seed can differ in the last few digits; and because the EM log-likelihood
#' can have close local optima, an occasional replicate converges to a
#' different one under a different backend. For replication material, record
#' the \code{future} plan and \code{utils::sessionInfo()} alongside
#' \code{boot_seed}, and produce the canonical set of estimates under a single
#' fixed plan.
#'
#' @return An object of class 'evzinb'
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' }
#'
#' # Peacekeeping and one-sided violence, with an in-formula log1p() transform
#' # for the Pareto (extreme-value) component (see `?hks`).
#' \dontrun{
#' data(hks)
#' hks_mod <- evzinb(
#'   osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
#'     lntpop + brv_AllLag_log + osvAllLagDum + incomp,
#'   formula_pareto = ~ log1p(troopLag),
#'   data = hks, n_bootstraps = 5, multicore = FALSE
#' )
#' summary(hks_mod)
#' glance(hks_mod)
#' predict(hks_mod, type = "harmonic")
#' }
evzinb <- function(
  formula_nb,
  formula_zi = NULL,
  formula_evi = NULL,
  formula_pareto = NULL,
  data,
  bootstrap = TRUE,
  n_bootstraps = 100,
  multicore = NULL,
  ncores = NULL,
  block = NULL,
  weights = NULL,
  time = NULL,
  bootstrap_scheme = NULL,
  block_length = NULL,
  family = evinf_family(),
  boot_seed = NULL,
  start_seed = NULL,
  control = evinf_control(),
  max.diff.par, max.no.em.steps, max.no.em.steps.warmup, c.lim, prune.c.range,
  max.upd.par.zc.multinomial, max.upd.par.pl.multinomial, max.upd.par.nb,
  max.upd.par.pl, no.m.bfgs.steps.multinomial, no.m.bfgs.steps.nb,
  no.m.bfgs.steps.pl, pdf.pl.type, eta.int, init.Beta.multinom.ZC,
  init.Beta.multinom.PL, init.Beta.NB, init.Beta.PL, init.Alpha.NB, init.C,
  verbose = FALSE
) {
  block <- evinf_block_name(rlang::enquo(block), parent.frame(), data)
  # round9 F: `time` accepts a bare column name / string, same resolution as
  # `block` (evinf_block_name() is generic despite its name).
  time <- evinf_block_name(rlang::enquo(time), parent.frame(), data)
  # round9 D.2: weights = accepts a bare column name, a string naming a
  # column, or a numeric vector; a raw vector is injected into `data` under a
  # reserved name so it survives na.omit() exactly like every other model
  # variable (evinf_resolve_weights()).
  weights_resolved <- evinf_resolve_weights(rlang::enquo(weights), parent.frame(), data)
  data <- weights_resolved$data
  weights_col <- weights_resolved$weights_col
  # round9 E.0: family = accepts an evinf_family() object or a plain string
  # (shorthand for the count family).
  family <- evinf_resolve_family(family)
  mc <- match.call()
  ctrl <- resolve_evinf_control(control, mc, environment(), fn = "evzinb")

  # A self-contained call for update(): the bootstrap-related arguments are
  # inlined by value and `data` is kept as the *expression* the user passed
  # (never the data frame itself, audit N6); update.evzinb() fills in
  # formulas / control / block from the fitted object.
  stored_call <- as.call(c(quote(evinf::evzinb), list(
    bootstrap = bootstrap, n_bootstraps = n_bootstraps, multicore = multicore,
    ncores = ncores, boot_seed = boot_seed, start_seed = start_seed,
    family = family, verbose = verbose
  )))
  stored_call$data <- mc$data

  # NULL component formulas are resolved in run_evzinb() (audit 4.4).
  full_run <- run_evzinb(
    formula_nb = formula_nb,
    formula_zi = formula_zi,
    formula_evi = formula_evi,
    formula_pareto = formula_pareto,
    data = data,
    control = ctrl,
    block = block,
    weights = weights_col,
    time = time,
    bootstrap_scheme = bootstrap_scheme,
    block_length = block_length,
    family = family,
    start_seed = start_seed,
    multicore = multicore,
    ncores = ncores,
    verbose = verbose
  )
  full_run$weights_col <- weights_col
  full_run$call <- stored_call

  block2 <- full_run$block_vec
  time2 <- full_run$time_vec

  if (bootstrap) {
    # Always record the seed actually used (draw one when the user passed NULL)
    # so object$boot_seeds is complete and add_bootstraps() can guard against a
    # reused seed (audit N4).
    if (is.null(boot_seed)) {
      boot_seed <- sample.int(.Machine$integer.max, 1L)
    }
    full_run$boot_seeds <- list(boot_seed)

    # Ship only a compact spec to the workers, not the whole fitted object.
    boot_spec <- evinf_boot_spec(full_run)
    boots <- evinf_with_plan(multicore, ncores, {
      # round10 J.1 (audit §5.11): project total runtime from the first
      # min(4, n_bootstraps) replicates actually completing, not the
      # full-sample fit's own time x n_bootstraps (the old proxy -- it
      # ignored parallelism, and a replicate's fit time can differ
      # substantially from the full-sample fit's).
      evinf_estimate_boot_runtime(
        bootrun_evzinb, boot_spec, block2, time2, n_bootstraps, boot_seed, verbose
      )
      evinf_pmap(
        seq_len(n_bootstraps),
        function(i, spec, blk, tv) try(bootrun_evzinb(spec, blk, tv)),
        spec = boot_spec, blk = block2, tv = time2,
        seed = boot_seed, label = "bootstrap", verbose = verbose,
        chunk_size = ctrl$chunk_size
      )
    })
    names(boots) <- paste0("bootstrap_", seq_along(boots))
    n_c_bnd <- evinf_warn_c_boundary(full_run, boots)
    out <- c(full_run, list(bootstraps = boots, n_c_on_boundary = n_c_bnd))
  } else {
    out <- full_run
  }
  class(out) <- 'evzinb'
  return(out)
}

#' A single bootstrap run on an evzinb object
#'
#' @param object The evzinb object to run the bootstrap on
#' @param block Optional string specifying varible for block bootstrapping
#' @param time_vec Optional time-index vector (round9 F), required for
#'   \code{object$bootstrap_scheme \%in\% c("moving_block", "stationary")}.
#' @param timing Should time be kept
#'
#' @return A bootstrapped evzinb object
#'
#' @noRd
bootrun_evzinb <- function(
  object,
  block = NULL,
  time_vec = NULL,
  timing = TRUE
) {
  tim <- Sys.time()
  # round9 F: the single resample-index generator behind every
  # bootstrap_scheme; "iid"/"cluster" reproduce the pre-F behaviour exactly.
  boot_id <- evinf_resample_ids(
    nrow(object$data$x.nb),
    scheme = object$bootstrap_scheme %||% (if (is.null(block)) "iid" else "cluster"),
    block_vec = block,
    time_vec = time_vec,
    block_length = object$block_length
  )
  OBS.Y <- object$data$y[boot_id]

  OBS.X.obj <- list()
  OBS.X.obj$X.multinom.ZC <- object$data$x.multinom.zc[boot_id, , drop = FALSE]
  OBS.X.obj$X.multinom.PL <- object$data$x.multinom.pl[boot_id, , drop = FALSE]
  OBS.X.obj$X.NB <- object$data$x.nb[boot_id, , drop = FALSE]
  OBS.X.obj$X.PL <- object$data$x.pl[boot_id, , drop = FALSE]
  OBS.X.obj$offset.nb <- if (is.null(object$offset_nb)) rep(0, length(boot_id)) else
    object$offset_nb[boot_id]
  OBS.X.obj$offset.zc <- if (is.null(object$offset_zc)) rep(0, length(boot_id)) else
    object$offset_zc[boot_id]
  OBS.X.obj$offset.pl_mult <- if (is.null(object$offset_pl_mult)) rep(0, length(boot_id)) else
    object$offset_pl_mult[boot_id]
  # round9 D.2: resample rows as above, carrying each row's weight along --
  # do not reweight the resampling probabilities themselves.
  OBS.X.obj$weights <- if (is.null(object$weights)) rep(1, length(boot_id)) else
    object$weights[boot_id]
  OBS.X.obj$has_weights <- isTRUE(object$has_weights)  # round10 0.5
  Control <- object$control

  Ini.Val <- list()
  Ini.Val$Beta.multinom.ZC <- as.numeric(object$coef$Beta.multinom.ZC)
  Ini.Val$Beta.multinom.PL <- as.numeric(object$coef$Beta.multinom.PL)
  Ini.Val$Beta.NB <- as.numeric(object$coef$Beta.NB)
  Ini.Val$Beta.PL <- as.numeric(object$coef$Beta.PL)
  # round9 E.1: object$coef$Alpha.NB is NULL for a Poisson count state (it
  # isn't a parameter); update_bfgs_fun() still needs *some* placeholder
  # double there (it's simply never updated/used), so fall back to the
  # control default.
  Ini.Val$Alpha.NB <- object$coef$Alpha.NB %||% object$control$init.Alpha.NB %||% 0.01
  Ini.Val$C <- object$coef$C
  capture.output(
    evzinb_boot <- em_fit(OBS.Y, OBS.X.obj, Ini.Val, Control, model = "evzinb",
                          full_sample = FALSE, family = object$family %||% evinf_family())
  )
  evzinb_boot$par.mat$Beta.multinom.ZC <- as.numeric(
    evzinb_boot$par.mat$Beta.multinom.ZC
  )
  evzinb_boot$par.mat$Beta.multinom.PL <- as.numeric(
    evzinb_boot$par.mat$Beta.multinom.PL
  )
  evzinb_boot$par.mat$Beta.NB <- as.numeric(evzinb_boot$par.mat$Beta.NB)
  evzinb_boot$par.mat$Beta.PL <- as.numeric(evzinb_boot$par.mat$Beta.PL)

  names(evzinb_boot$par.mat$Beta.NB) <- c(
    '(Intercept)',
    colnames(object$data$x.nb)
  )
  names(evzinb_boot$par.mat$Beta.multinom.ZC) <- c(
    '(Intercept)',
    colnames(object$data$x.multinom.zc)
  )
  names(evzinb_boot$par.mat$Beta.multinom.PL) <- c(
    '(Intercept)',
    colnames(object$data$x.multinom.pl)
  )
  names(evzinb_boot$par.mat$Beta.PL) <- c(
    '(Intercept)',
    colnames(object$data$x.pl)
  )

  evzinb_boot$props <- evzinb_boot$par.mat$Props

  evzinb_boot$par.mat$Props <- NULL
  # round9 E.1: mirror run_evzinb()'s drop, so a bootstrap replicate's coef
  # is absent Alpha.NB exactly like the full-sample fit's.
  if (isTRUE(evzinb_boot$family$count == "poisson")) {
    evzinb_boot$par.mat$Alpha.NB <- NULL
  }
  evzinb_boot$coef <- evzinb_boot$par.mat
  evzinb_boot$par.mat <- NULL
  evzinb_boot$resp <- NULL

  evzinb_boot$formulas <- object$formulas
  evzinb_boot$terms <- object$terms
  evzinb_boot$xlevels <- object$xlevels

  evzinb_boot$y.hat.plexpElogy <- NULL
  evzinb_boot$y.hat.pl.E.inv.y <- NULL
  evzinb_boot$y.hat.plmedian <- NULL
  evzinb_boot$y.hat.plmean <- NULL
  evzinb_boot$mu.nb.vec <- NULL
  evzinb_boot$alpha.pl.vec <- NULL
  evzinb_boot$exp.E.log.y <- NULL
  evzinb_boot$median.pl.vec <- NULL
  evzinb_boot$mean.pl.vec <- NULL
  evzinb_boot$x.nb <- NULL
  evzinb_boot$x.pl <- NULL
  evzinb_boot$x.multinom.pl <- NULL
  evzinb_boot$x.multinom.zc <- NULL
  evzinb_boot$c_profile <- NULL  # keep bootstraps small; c_trace is enough (4.5)
  # round10 G.3 (audit §5.10): converge/c_converged/c_warmup_capped already
  # survive here as-is (em_fit()'s own return, never nulled below); n_em_steps
  # is the one derived field run_evzinb() computes that bootstraps didn't get.
  evzinb_boot$n_em_steps <- length(evzinb_boot$log.lik.vec.all)

  evzinb_boot$data <- NULL
  evzinb_boot <- evinf_flag_degenerate(evzinb_boot, OBS.X.obj$X.PL, Control)
  evzinb_boot$boot_id <- boot_id
  # round9 F: with overlapping blocks (moving_block/stationary) the not-drawn
  # set is both smaller and more correlated than under iid resampling, so an
  # OOB error computed from it is optimistic -- report the realised fraction
  # alongside boot_id rather than leaving it implicit.
  evzinb_boot$oob_fraction <- length(setdiff(seq_len(nrow(object$data$x.nb)), unique(boot_id))) /
    nrow(object$data$x.nb)
  if (timing) {
    evzinb_boot$time <- difftime(Sys.time(), tim, units = 'secs')
  }
  return(evzinb_boot)
}
