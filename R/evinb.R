#' Running an extreme value inflated negative binomial model
#'
#' @param formula_nb,formula_evi,formula_pareto Component formulas.
#' @param data Data to run the model on.
#' @param control An \code{evinf_control()} object.
#' @param block Optional string naming a case-identifier column for block bootstrapping; included in the na.omit() so the returned block vector aligns with the model data.
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
  verbose = FALSE
) {
  control <- validate_evinf_control(control)
  if (!is.null(block) && !(is.character(block) && length(block) == 1L)) {
    stop("`block` must be NULL or a single string naming a column of `data`.",
         call. = FALSE)
  }
  # audit 4.4: offsets are only supported in the count component.
  formula_evi <- evinf_component_formula(formula_evi, formula_nb, "formula_evi")
  formula_pareto <- evinf_component_formula(formula_pareto, formula_nb, "formula_pareto")

  # Restrict to the union of variables used by any component (plus the block
  # variable, audit R0.1) and drop incomplete rows once, so the design matrices
  # and the block vector all refer to the same rows (audit 1.7, R0.1).
  model_vars <- unique(c(
    all.vars(formula_nb),
    all.vars(formula_evi),
    all.vars(formula_pareto),
    block
  ))
  model_data <- data %>%
    dplyr::select(dplyr::all_of(model_vars)) %>%
    na.omit()

  d_nb <- evinf_design(formula_nb, model_data)
  d_evi <- evinf_design(formula_evi, model_data)
  d_pareto <- evinf_design(formula_pareto, model_data)
  offset_nb <- d_nb$offset

  OBS.Y <- as.matrix(model.response(model.frame(formula_nb, model_data)))

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

  if (verbose) {
    object <- em_fit(OBS.Y, OBS.X.obj, Ini.Val, Control, model = "evinb")
  } else {
    capture.output(
      object <- em_fit(OBS.Y, OBS.X.obj, Ini.Val, Control, model = "evinb")
    )
  }
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
  object$data <- list()

  object$data$data <- model_data
  object$block <- block
  object$block_vec <- if (!is.null(block)) model_data[[block]] else NULL

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
  object$coef <- object$par.mat
  object$par.mat <- NULL
  object$coef$Beta.multinom.ZC <- NULL

  # The zero-inflation component is switched off (its intercept is fixed), so its
  # coefficients are not free parameters: drop them from par.all and recompute
  # the information criteria from the reduced count.
  n_zc <- length(Ini.Val$Beta.multinom.ZC)
  object$par.all <- object$par.all[-seq_len(n_zc)]
  n_par <- length(object$par.all)
  n_obs <- length(object$data$y)
  object$AIC <- 2 * n_par - 2 * object$log.lik
  object$BIC <- log(n_obs) * n_par - 2 * object$log.lik

  # audit 4.5: threshold diagnostics
  object$c_profile <- tibble::as_tibble(object$c_profile)
  object$loglik_trace <- object$log.lik.vec.all
  object$n_above_c <- sum(object$data$y >= object$coef$C)
  object$n_em_steps <- length(object$log.lik.vec.all)

  object$fitted <- list()
  object$fitted$y.hat.pl_exp.E.logy <- object$y.hat.plexpElogy
  object$y.hat.plexpElogy <- NULL
  object$fitted$y.hat.pl_E.inv.y <- object$y.hat.pl.E.inv.y
  object$y.hat.pl.E.inv.y <- NULL
  object$fitted$y.hat.pl_median <- object$y.hat.plmedian
  object$y.hat.plmedian <- NULL
  object$fitted$y.hat.pl_mean <- object$y.hat.plmean
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
#' @param timing Should time be kept
#'
#' @return A bootstrapped evinb object
#'
#' @noRd
bootrun_evinb <- function(
  object,
  block = NULL,
  timing = TRUE
) {
  tim <- Sys.time()
  if (is.null(block)) {
    boot_id <- sample(
      1:nrow(object$data$x.nb),
      nrow(object$data$x.nb),
      replace = T
    )
  } else {
    uniques <- unique(block)
    boot_block_id <- sample(uniques, length(uniques), replace = T)
    boot_id <- boot_block_id %>%
      purrr::map(~ which(block == .x)) %>%
      purrr::reduce(c)
  }
  OBS.Y <- object$data$y[boot_id]

  OBS.X.obj <- list()
  OBS.X.obj$X.multinom.ZC <- object$data$x.nb[boot_id, , drop = FALSE] # note: only for functionality, does not update
  OBS.X.obj$X.multinom.PL <- object$data$x.multinom.pl[boot_id, , drop = FALSE]
  OBS.X.obj$X.NB <- object$data$x.nb[boot_id, , drop = FALSE]
  OBS.X.obj$X.PL <- object$data$x.pl[boot_id, , drop = FALSE]
  OBS.X.obj$offset.nb <- if (is.null(object$offset_nb)) rep(0, length(boot_id)) else
    object$offset_nb[boot_id]
  Control <- object$control

  Ini.Val <- list()
  Ini.Val$Beta.multinom.ZC <- c(-100, rep(0, ncol(OBS.X.obj$X.multinom.ZC)))
  Ini.Val$Beta.multinom.PL <- as.numeric(object$coef$Beta.multinom.PL)
  Ini.Val$Beta.NB <- as.numeric(object$coef$Beta.NB)
  Ini.Val$Beta.PL <- as.numeric(object$coef$Beta.PL)
  Ini.Val$Alpha.NB <- object$coef$Alpha.NB
  Ini.Val$C <- object$coef$C
  capture.output(
    evinb_boot <- em_fit(OBS.Y, OBS.X.obj, Ini.Val, Control, model = "evinb")
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

  evinb_boot$data <- NULL
  evinb_boot$boot_id <- boot_id

  class(evinb_boot) <- 'evinb_boot'
  if (timing) {
    evinb_boot$time <- difftime(Sys.time(), tim, units = 'secs')
  }
  return(evinb_boot)
}

#' Running an extreme value inflated negative binomial model with bootstrapping
#'
#' @param formula_nb Formula for the negative binomial (count) component of the model
#' @param formula_evi Formula for the extreme-value inflation component of the model. If NULL taken as the same formula as nb
#' @param formula_pareto Formula for the pareto (extreme value) component of the model. If NULL taken as the same formula as nb
#' @param data Data to run the model on
#' @param bootstrap Should bootstrapping be performed. Needed to obtain standard errors and p-values
#' @param n_bootstraps Number of bootstraps to run. For use of bootstrapped p-values, at least 1,000 bootstraps are recommended. For approximate p-values, a lower number can be sufficient
#' @param multicore Should multiple cores be used?
#' @param ncores Number of cores if multicore is used. Default (NULL) is one less than the available number of cores
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
#' @param control An \code{\link{evinf_control}()} object holding the EM tuning
#'   settings.
#' @param max.diff.par,max.no.em.steps,max.no.em.steps.warmup,c.lim,prune.c.range,max.upd.par.pl.multinomial,max.upd.par.nb,max.upd.par.pl,no.m.bfgs.steps.multinomial,no.m.bfgs.steps.nb,no.m.bfgs.steps.pl,pdf.pl.type,eta.int,init.Beta.multinom.PL,init.Beta.NB,init.Beta.PL,init.Alpha.NB,init.C
#'   \strong{Deprecated.} Pass these through \code{control = evinf_control(...)};
#'   supplying one directly overrides the corresponding \code{control} element and
#'   emits a warning. See \code{\link{evinf_control}}.
#' @param verbose Should progress be printed for the first run of evinb
#'
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
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
  multicore = FALSE,
  ncores = NULL,
  block = NULL,
  boot_seed = NULL,
  control = evinf_control(),
  max.diff.par, max.no.em.steps, max.no.em.steps.warmup, c.lim, prune.c.range,
  max.upd.par.pl.multinomial, max.upd.par.nb, max.upd.par.pl,
  no.m.bfgs.steps.multinomial, no.m.bfgs.steps.nb, no.m.bfgs.steps.pl,
  pdf.pl.type, eta.int, init.Beta.multinom.PL, init.Beta.NB, init.Beta.PL,
  init.Alpha.NB, init.C,
  verbose = FALSE
) {
  i <- 'temp_iter'
  block <- evinf_block_name(rlang::enquo(block), parent.frame(), data)
  mc <- match.call()
  ctrl <- resolve_evinf_control(control, mc, environment(), fn = "evinb")

  stored_call <- as.call(c(quote(evinf::evinb), list(
    bootstrap = bootstrap, n_bootstraps = n_bootstraps, multicore = multicore,
    ncores = ncores, boot_seed = boot_seed, verbose = verbose
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
    verbose = verbose
  )
  full_run$call <- stored_call
  runtime <- difftime(Sys.time(), t1)

  block2 <- full_run$block_vec

  if (bootstrap) {
    # Always record the seed actually used (audit N4).
    if (is.null(boot_seed)) {
      boot_seed <- sample.int(.Machine$integer.max, 1L)
    }
    full_run$boot_seeds <- list(boot_seed)

    be <- evinf_setup_backend(multicore, ncores)
    on.exit(be$stop(), add = TRUE)
    if (multicore) {
      ex_time <- runtime * n_bootstraps / be$ncores
    } else {
      ex_time <- runtime * n_bootstraps
    }
    if (verbose) {
      cat(
        "\n ======",
        "Approximate runtime for bootstraps is",
        ex_time,
        attributes(runtime)$units,
        ". Note: This is a very rough estimate of the runtime."
      )
    }

    boots <- evinf_progress_run(n_bootstraps, verbose, function(p) {
      foreach::foreach(
        i = 1:n_bootstraps,
        .options.RNG = boot_seed,
        .packages = 'evinf',
        .export = c('full_run', 'block2', 'n_bootstraps')
      ) %dorng% {
        res <- try(bootrun_evinb(full_run, block2))
        p()
        res
      }
    })

    names(boots) <- paste('bootstrap_', 1:length(boots), sep = "")
    out <- c(full_run, list(bootstraps = boots))
  } else {
    out <- full_run
  }

  class(out) <- 'evinb'
  return(out)
}
