#' Running an extreme value inflated negative binomial model
#'
#' @param formula_nb,formula_evi,formula_pareto Component formulas.
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
#'   \code{zero = "hurdle"} errors: \code{evinb()} has no zero state to
#'   hurdle over.
#' @param verbose Should progress be printed for the first run of evinb.
#'
#' @return An object of class 'evinb'
#' @noRd
run_evinb <- function(
  formula_nb,
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
  verbose = FALSE
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
  if (family$zero == "hurdle") {
    stop("evinb() has no zero state to hurdle over (family = ",
         "evinf_family(zero = \"hurdle\") is evzinb()-only).", call. = FALSE)
  }
  # round9 D.1 (audit §5.6): offset() is supported in the count and
  # EVI-inflation components; not in the Pareto shape component. evinb has no
  # zero-inflation component at all.
  formula_evi <- evinf_component_formula(formula_evi, formula_nb, "formula_evi")
  formula_pareto <- evinf_component_formula(formula_pareto, formula_nb, "formula_pareto",
                                            allow_offset = FALSE)

  # Restrict to the union of variables used by any component (plus the block
  # variable, audit R0.1) and drop incomplete rows once, so the design matrices
  # and the block vector all refer to the same rows (audit 1.7, R0.1).
  model_vars <- unique(c(
    all.vars(formula_nb),
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
  d_evi <- evinf_design(formula_evi, model_data)
  d_pareto <- evinf_design(formula_pareto, model_data)
  offset_nb <- d_nb$offset
  offset_pl_mult <- d_evi$offset
  weights_vec <- if (!is.null(weights)) model_data[[weights]] else rep(1, nrow(model_data))
  evinf_check_weights(weights_vec, nrow(model_data))

  OBS.Y <- as.matrix(model.response(model.frame(formula_nb, model_data)))
  evinf_check_response(as.numeric(OBS.Y))

  # Resolve NULL c.lim / init.C against the data (audit 4.5).
  control <- evinf_resolve_c(control, as.numeric(OBS.Y))
  # The zero-inflation component is switched off in evinb, so its step size
  # mirrors the extreme-value-inflation one.
  Control <- control
  Control$max.upd.par.zc.multinomial <- control$max.upd.par.pl.multinomial

  init.Beta.multinom.PL <- control$init.Beta.multinom.PL
  init.Beta.NB <- control$init.Beta.NB
  init.Beta.PL <- control$init.Beta.PL

  OBS.X.obj <- list()
  # The zero-inflation component is switched off in evinb (intercept fixed very
  # negative); its design just mirrors the NB design for shape.
  OBS.X.obj$X.multinom.ZC <- d_nb$X
  OBS.X.obj$X.multinom.PL <- d_evi$X
  OBS.X.obj$X.NB <- d_nb$X
  OBS.X.obj$X.PL <- d_pareto$X
  OBS.X.obj$offset.nb <- offset_nb
  OBS.X.obj$offset.pl_mult <- offset_pl_mult
  OBS.X.obj$weights <- weights_vec
  # round10 0.5: `weights` here is still the R-level argument (a column name
  # or NULL), not weights_vec -- has_weights is TRUE only when the user
  # actually supplied weights, never merely because weights_vec is all 1s.
  OBS.X.obj$has_weights <- !is.null(weights)

  # Parameter counts include the intercept the C++ routines prepend.
  n_nb <- ncol(d_nb$X) + 1L
  n_pl_mult <- ncol(d_evi$X) + 1L
  n_pl <- ncol(d_pareto$X) + 1L

  Ini.Val <- list()
  Ini.Val$Beta.multinom.ZC <- rep(0, n_nb)
  Ini.Val$Beta.multinom.ZC[1] <- -100

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
  starts_result <- evinf_with_plan(multicore, ncores, {
    evinf_run_starts(OBS.Y, OBS.X.obj, Ini.Val, Control,
                     model = "evinb", family = family,
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
  #names(object$par.mat$Beta.multinom.ZC) <- c('(Intercept)',all.vars(formula_zi)[-1])
  names(object$par.mat$Beta.multinom.PL) <- c('(Intercept)', colnames(d_evi$X))
  names(object$par.mat$Beta.PL) <- c('(Intercept)', colnames(d_pareto$X))

  object$formulas <- list(
    formula_nb = formula_nb,
    #formula_zi = formula_zi,
    formula_evi = formula_evi,
    formula_pareto = formula_pareto
  )
  object$terms <- list(
    nb = d_nb$terms,
    evi = d_evi$terms,
    pareto = d_pareto$terms
  )
  object$xlevels <- list(
    nb = d_nb$xlevels,
    evi = d_evi$xlevels,
    pareto = d_pareto$xlevels
  )
  object$c_lim_default <- isTRUE(control$c_lim_default)
  object$has_offset <- isTRUE(d_nb$has_offset)
  object$offset_nb <- offset_nb
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
  # round9 F: see the matching comment in run_evzinb().
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
  object$data$x.multinom.pl <- object$x.multinom.pl
  object$x.multinom.pl <- NULL
  object$resp <- object$resp[, 2:3]

  object$props <- object$par.mat$Props[, 2:3]
  colnames(object$props) <- colnames(object$resp) <- c('count', 'evi')
  object$par.mat$Props <- NULL
  # round9 E.1: mirror run_evzinb()'s drop, so coef()/vcov()/confint()/
  # tidy()/glance() all show Alpha.NB as absent (not NA) for a Poisson count
  # state.
  if (family$count == "poisson") {
    object$par.mat$Alpha.NB <- NULL
  }
  object$coef <- object$par.mat
  object$par.mat <- NULL
  object$coef$Beta.multinom.ZC <- NULL

  # The zero-inflation component is switched off (its intercept is fixed), so its
  # coefficients are not free parameters: drop them from par.all and recompute
  # the information criteria from the reduced count.
  n_zc <- length(Ini.Val$Beta.multinom.ZC)
  object$par.all <- object$par.all[-seq_len(n_zc)]
  n_par <- length(object$par.all)
  n_obs <- sum(object$weights)  # round9 D.2: sum(weights), row count when unweighted
  object$AIC <- 2 * n_par - 2 * object$log.lik
  object$BIC <- log(n_obs) * n_par - 2 * object$log.lik

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
  # round10 0.6 (review §4/§7, breaking change): see the matching comment in
  # run_evzinb() -- the legacy y.hat.pl_* point predictions duplicated
  # predict() but were not hurdle-aware, and no internal reader of
  # object$fitted$y.hat.pl_* remained.
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
  object$fitted$prob_count <- object$props[, 1]
  object$fitted$prob_evi <- object$props[, 2]
  object$fitted$posterior_count <- object$resp[, 1]
  object$fitted$posterior_evi <- object$resp[, 2]
  # Deprecated pre-0.9.4 names, kept as duplicates for one release (audit 2.7).
  object$fitted$prob_pareto <- object$props[, 2]
  object$fitted$posterior_pareto <- object$resp[, 2]

  class(object) <- 'evinb'
  return(object)
}


#' A single bootstrap run on an ezinb object
#'
#' @param object The evinb object to run the bootstrap on
#' @param block Optional string specifying varible for block bootstrapping
#' @param time_vec Optional time-index vector (round9 F), required for
#'   \code{object$bootstrap_scheme \%in\% c("moving_block", "stationary")}.
#' @param timing Should time be kept
#'
#' @return A bootstrapped evinb object
#'
#' @noRd
bootrun_evinb <- function(
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
  OBS.X.obj$X.multinom.ZC <- object$data$x.nb[boot_id, , drop = FALSE] # note: only for functionality, does not update
  OBS.X.obj$X.multinom.PL <- object$data$x.multinom.pl[boot_id, , drop = FALSE]
  OBS.X.obj$X.NB <- object$data$x.nb[boot_id, , drop = FALSE]
  OBS.X.obj$X.PL <- object$data$x.pl[boot_id, , drop = FALSE]
  OBS.X.obj$offset.nb <- if (is.null(object$offset_nb)) rep(0, length(boot_id)) else
    object$offset_nb[boot_id]
  OBS.X.obj$offset.pl_mult <- if (is.null(object$offset_pl_mult)) rep(0, length(boot_id)) else
    object$offset_pl_mult[boot_id]
  OBS.X.obj$weights <- if (is.null(object$weights)) rep(1, length(boot_id)) else
    object$weights[boot_id]
  OBS.X.obj$has_weights <- isTRUE(object$has_weights)  # round10 0.5
  Control <- object$control

  Ini.Val <- list()
  Ini.Val$Beta.multinom.ZC <- c(-100, rep(0, ncol(OBS.X.obj$X.multinom.ZC)))
  Ini.Val$Beta.multinom.PL <- as.numeric(object$coef$Beta.multinom.PL)
  Ini.Val$Beta.NB <- as.numeric(object$coef$Beta.NB)
  Ini.Val$Beta.PL <- as.numeric(object$coef$Beta.PL)
  # round9 E.1: see the matching comment in bootrun_evzinb().
  Ini.Val$Alpha.NB <- object$coef$Alpha.NB %||% object$control$init.Alpha.NB %||% 0.01
  Ini.Val$C <- object$coef$C
  capture.output(
    evinb_boot <- em_fit(OBS.Y, OBS.X.obj, Ini.Val, Control, model = "evinb",
                         full_sample = FALSE, family = object$family %||% evinf_family())
  )

  evinb_boot$par.mat$Beta.multinom.PL <- as.numeric(
    evinb_boot$par.mat$Beta.multinom.PL
  )
  evinb_boot$par.mat$Beta.NB <- as.numeric(evinb_boot$par.mat$Beta.NB)
  evinb_boot$par.mat$Beta.PL <- as.numeric(evinb_boot$par.mat$Beta.PL)

  names(evinb_boot$par.mat$Beta.NB) <- c(
    '(Intercept)',
    colnames(object$data$x.nb)
  )
  #names(evinb_boot$par.mat$Beta.multinom.ZC) <- c('(Intercept)',all.vars(object$formulas$formula_zi)[-1])
  names(evinb_boot$par.mat$Beta.multinom.PL) <- c(
    '(Intercept)',
    colnames(object$data$x.multinom.pl)
  )
  names(evinb_boot$par.mat$Beta.PL) <- c(
    '(Intercept)',
    colnames(object$data$x.pl)
  )

  evinb_boot$resp <- NULL
  evinb_boot$props <- evinb_boot$par.mat$Props[, 2:3]

  evinb_boot$par.mat$Props <- NULL
  if (isTRUE(evinb_boot$family$count == "poisson")) {
    evinb_boot$par.mat$Alpha.NB <- NULL
  }

  evinb_boot$coef <- evinb_boot$par.mat
  evinb_boot$par.mat <- NULL
  evinb_boot$coef$Beta.multinom.ZC <- NULL
  evinb_boot$formulas <- object$formulas
  evinb_boot$terms <- object$terms
  evinb_boot$xlevels <- object$xlevels

  evinb_boot$y.hat.plexpElogy <- NULL
  evinb_boot$y.hat.pl.E.inv.y <- NULL
  evinb_boot$y.hat.plmedian <- NULL
  evinb_boot$y.hat.plmean <- NULL
  evinb_boot$mu.nb.vec <- NULL
  evinb_boot$alpha.pl.vec <- NULL
  evinb_boot$exp.E.log.y <- NULL
  evinb_boot$median.pl.vec <- NULL
  evinb_boot$mean.pl.vec <- NULL
  evinb_boot$c_profile <- NULL  # keep bootstraps small; c_trace is enough (4.5)
  # round10 G.3 (audit §5.10): converge/c_converged/c_warmup_capped already
  # survive here as-is (em_fit()'s own return, never nulled below); n_em_steps
  # is the one derived field run_evinb() computes that bootstraps didn't get.
  evinb_boot$n_em_steps <- length(evinb_boot$log.lik.vec.all)

  evinb_boot$data <- NULL
  evinb_boot <- evinf_flag_degenerate(evinb_boot, OBS.X.obj$X.PL, Control)
  evinb_boot$boot_id <- boot_id
  # round9 F: see the matching comment in bootrun_evzinb().
  evinb_boot$oob_fraction <- length(setdiff(seq_len(nrow(object$data$x.nb)), unique(boot_id))) /
    nrow(object$data$x.nb)

  class(evinb_boot) <- 'evinb_boot'
  if (timing) {
    evinb_boot$time <- difftime(Sys.time(), tim, units = 'secs')
  }
  return(evinb_boot)
}

#' Running an extreme value inflated negative binomial model with bootstrapping
#'
#' @param formula_nb Formula for the negative binomial (count) component of the model.
#'   May include an \code{offset()} term (\eqn{\mu_{NB} = \exp(x'\beta + offset)}).
#' @param formula_evi Formula for the extreme-value inflation component of the model.
#'   If NULL taken as the same formula as nb, with any \code{offset()} term stripped
#'   (an offset applies only where it is written explicitly, never by inheritance).
#'   May include its own \code{offset()} term: the EVI logit is the log-odds of the
#'   extreme-value state against the count state, so an offset there shifts that
#'   log-odds, e.g. \code{offset(log(population))} for a probability of an extreme
#'   event that scales with exposure.
#' @param formula_pareto Formula for the pareto (extreme value) component of the
#'   model. If NULL taken as the same formula as nb, offset stripped as above.
#'   \code{offset()} is \strong{not} supported here (errors if present): an offset
#'   on a shape parameter has no clear reading.
#' @param data Data to run the model on
#' @param bootstrap Should bootstrapping be performed. Needed to obtain standard errors and p-values
#' @param n_bootstraps Number of bootstraps to run. For use of bootstrapped p-values, at least 1,000 bootstraps are recommended. For approximate p-values, a lower number can be sufficient
#' @inheritParams evzinb
#' @param block Optional case-identifier column for block bootstrapping, given
#'   either as a bare column name (\code{block = id}) or a string
#'   (\code{block = "id"}). Note that the bundled \code{\link{hks}} data contain
#'   no conflict identifier, so the conflict-level cluster bootstrap in Randahl
#'   and Vegelius (2024) cannot be reproduced from them directly (see
#'   \code{?hks}).
#' @param boot_seed Optional bootstrap seed for reproducibility. When supplied
#'   it is used as-is; when \code{NULL} a seed is drawn and recorded, so
#'   \code{object$boot_seeds} is always populated for a bootstrapped model
#'   (see \code{\link{add_bootstraps}}).
#' @param start_seed Optional seed for the perturbed starts when
#'   \code{control$n_starts > 1} (round10 G.1); recorded as
#'   \code{object$start_seed} (\code{NULL} for the default \code{n_starts =
#'   1}). Unused otherwise.
#' @param control An \code{\link{evinf_control}()} object holding the EM tuning
#'   settings.
#' @param max.diff.par,max.no.em.steps,max.no.em.steps.warmup,c.lim,prune.c.range,max.upd.par.pl.multinomial,max.upd.par.nb,max.upd.par.pl,no.m.bfgs.steps.multinomial,no.m.bfgs.steps.nb,no.m.bfgs.steps.pl,pdf.pl.type,eta.int,init.Beta.multinom.PL,init.Beta.NB,init.Beta.PL,init.Alpha.NB,init.C
#'   \strong{Deprecated.} Pass these through \code{control = evinf_control(...)};
#'   supplying one directly overrides the corresponding \code{control} element and
#'   emits a warning. See \code{\link{evinf_control}}. Will be removed in 0.11.0.
#' @param verbose Should progress be printed for the first run of evinb
#'
#' @inheritSection evzinb Parallel processing
#' @inheritSection evzinb Reproducibility
#'
#' @return An object of class 'evinb'
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' }
#'
#' \dontrun{
#' data(hks)
#' hks_mod <- evinb(
#'   osvAll ~ troopLag + policeLag + militaryobserversLag + epduration +
#'     lntpop + brv_AllLag_log + osvAllLagDum + incomp,
#'   formula_pareto = ~ log1p(troopLag),
#'   data = hks, n_bootstraps = 5, multicore = FALSE
#' )
#' }
evinb <- function(
  formula_nb,
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
  max.upd.par.pl.multinomial, max.upd.par.nb, max.upd.par.pl,
  no.m.bfgs.steps.multinomial, no.m.bfgs.steps.nb, no.m.bfgs.steps.pl,
  pdf.pl.type, eta.int, init.Beta.multinom.PL, init.Beta.NB, init.Beta.PL,
  init.Alpha.NB, init.C,
  verbose = FALSE
) {
  block <- evinf_block_name(rlang::enquo(block), parent.frame(), data)
  time <- evinf_block_name(rlang::enquo(time), parent.frame(), data)
  weights_resolved <- evinf_resolve_weights(rlang::enquo(weights), parent.frame(), data)
  data <- weights_resolved$data
  weights_col <- weights_resolved$weights_col
  family <- evinf_resolve_family(family)
  mc <- match.call()
  ctrl <- resolve_evinf_control(control, mc, environment(), fn = "evinb")

  stored_call <- as.call(c(quote(evinf::evinb), list(
    bootstrap = bootstrap, n_bootstraps = n_bootstraps, multicore = multicore,
    ncores = ncores, boot_seed = boot_seed, start_seed = start_seed,
    family = family, verbose = verbose
  )))
  stored_call$data <- mc$data  # the expression, not the data frame (audit N6)

  # NULL component formulas are resolved in run_evinb() (audit 4.4).
  t1 <- Sys.time()
  full_run <- run_evinb(
    formula_nb = formula_nb,
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
  runtime <- difftime(Sys.time(), t1)

  block2 <- full_run$block_vec
  time2 <- full_run$time_vec

  if (bootstrap) {
    # Always record the seed actually used (audit N4).
    if (is.null(boot_seed)) {
      boot_seed <- sample.int(.Machine$integer.max, 1L)
    }
    full_run$boot_seeds <- list(boot_seed)

    if (verbose) {
      cat(
        "\n ======",
        "Approximate runtime for bootstraps is",
        runtime * n_bootstraps,
        attributes(runtime)$units,
        ". Note: This is a very rough estimate of the runtime (sequential)."
      )
    }

    boot_spec <- evinf_boot_spec(full_run)
    boots <- evinf_with_plan(multicore, ncores, {
      evinf_pmap(
        seq_len(n_bootstraps),
        function(i, spec, blk, tv) try(bootrun_evinb(spec, blk, tv)),
        spec = boot_spec, blk = block2, tv = time2,
        seed = boot_seed, label = "bootstrap", verbose = verbose
      )
    })
    names(boots) <- paste0("bootstrap_", seq_along(boots))
    n_c_bnd <- evinf_warn_c_boundary(full_run, boots)
    out <- c(full_run, list(bootstraps = boots, n_c_on_boundary = n_c_bnd))
  } else {
    out <- full_run
  }

  class(out) <- 'evinb'
  return(out)
}
