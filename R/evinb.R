#' Running an extreme value and zero inflated negative binomial model
#'
#' @param formula_nb Formula for the negative binomial (count) component of the model
#' @param formula_evi Formula for the extreme-value inflation component of the model. If NULL taken as the same formula as nb
#' @param formula_pareto Formula for the pareto (extreme value) component of the model If NULL taken as the same formula as nb
#' @param data Data to run the model on
#' @param max.diff.par Tolerance for EM algorithm. Will be considered to have converged if the maximum absolute difference in the parameter estimates are lower than this value
#' @param max.no.em.steps Maximum number of EM steps to run. Will be considered to not have converged if this number is reached and convergence is not reached
#' @param max.no.em.steps.warmup Number of EM steps in the warmup rounds
#' @param c.lim Integer range defining the possible values of C
#' @param prune.c.range The pruning fraction for the C range. Useful if there are many unique values in c.lim and the model is slow to fit. If FALSE, no pruning is done. If numeric [0-1], the C range is pruned to approximately length(c.lim)*prune.c.range
#' @param max.upd.par.pl.multinomial Maximum parameter change step size in the extreme value inflation component
#' @param max.upd.par.nb Maximum parameter change step size in the count component
#' @param max.upd.par.pl Maximum parameter change step size in the pareto component
#' @param no.m.bfgs.steps.multinomial Number of BFGS steps for the multinomial model
#' @param no.m.bfgs.steps.nb Number of BFGS steps for the negative binomial model
#' @param no.m.bfgs.steps.pl Number of BFGS steps for the pareto model
#' @param pdf.pl.type Probability density function type for the pareto component. Either 'approx' or 'exact'. 'approx' is adviced in most cases
#' @param eta.int Initial values for eta
#' @param init.Beta.multinom.PL Initial values for beta parameters in the extreme value inflation component. Vector of same length as number of parameters in the extreme value inflation component or NULL (which gives starting values of 0)
#' @param init.Beta.NB Initial values for beta parameters in the count component. Vector of same length as number of parameters in the count component or NULL (which gives starting values of 0)
#' @param init.Beta.PL Initial values for beta parameters in the pareto component. Vector of same length as number of parameters in the pareto component or NULL (which gives starting values of 0)
#' @param init.Alpha.NB Initial value of Alpha NB, integer or NULL (giving a starting value of 0)
#' @param init.C Initial value of C. Integer which should be within the C_lim range.
#'
#' @return An object of class 'evinf' containing XX
#' @noRd
run_evinb <- function(
  formula_nb,
  formula_evi = NULL,
  formula_pareto = NULL,
  data,
  max.diff.par = 1e-3,
  max.no.em.steps = 200,
  max.no.em.steps.warmup = 5,
  c.lim = c(70, 300),
  prune.c.range = FALSE,
  max.upd.par.pl.multinomial = 0.5,
  max.upd.par.nb = 0.5,
  max.upd.par.pl = 0.5,
  no.m.bfgs.steps.multinomial = 3,
  no.m.bfgs.steps.nb = 3,
  no.m.bfgs.steps.pl = 3,
  pdf.pl.type = c("approx", 'exact'),
  eta.int = c(-1, 1),
  init.Beta.multinom.PL = NULL,
  init.Beta.NB = NULL,
  init.Beta.PL = NULL,
  init.Alpha.NB = 0.01,
  init.C = 200,
  verbose = FALSE
) {
  if (is.null(formula_evi)) {
    formula_evi <- formula_nb
  }
  if (is.null(formula_pareto)) {
    formula_pareto <- formula_nb
  }
  pdf.pl.type <- match.arg(pdf.pl.type, c("approx", 'exact'))

  mf_nb <- model.frame(formula_nb, data)
  mf_zi <- model.frame(formula_nb, data)
  mf_evi <- model.frame(formula_evi, data)
  mf_pareto <- model.frame(formula_pareto, data)

  OBS.Y <- as.matrix(model.response(mf_nb))

  OBS.X.obj <- list()
  OBS.X.obj$X.multinom.ZC <- as.matrix(mf_zi[, -1])
  OBS.X.obj$X.multinom.PL <- as.matrix(mf_evi[, -1])
  OBS.X.obj$X.NB <- as.matrix(mf_nb[, -1])
  OBS.X.obj$X.PL <- as.matrix(mf_pareto[, -1])
  Control <- list(
    max.diff.par = max.diff.par,
    max.no.em.steps = max.no.em.steps,
    max.no.em.steps.warmup = max.no.em.steps.warmup,
    c.lim = c.lim,
    max.upd.par.zc.multinomial = max.upd.par.pl.multinomial,
    max.upd.par.pl.multinomial = max.upd.par.pl.multinomial,
    max.upd.par.nb = max.upd.par.nb,
    max.upd.par.pl = max.upd.par.pl,
    no.m.bfgs.steps.multinomial = no.m.bfgs.steps.multinomial,
    no.m.bfgs.steps.nb = no.m.bfgs.steps.nb,
    no.m.bfgs.steps.pl = no.m.bfgs.steps.pl,
    pdf.pl.type = pdf.pl.type,
    eta.int = eta.int,
    prune.c.range = prune.c.range
  )

  Ini.Val <- list()
  Ini.Val$Beta.multinom.ZC <- rep(0, ncol(mf_zi))
  Ini.Val$Beta.multinom.ZC[1] <- -100

  if (
    is.null(init.Beta.multinom.PL) |
      length(init.Beta.multinom.PL) != ncol(mf_evi)
  ) {
    Ini.Val$Beta.multinom.PL <- rep(0, ncol(mf_evi))
  } else {
    Ini.Val$Beta.multinom.PL <- init.Beta.multinom.PL
  }

  if (is.null(init.Beta.NB) | length(init.Beta.NB) != ncol(mf_evi)) {
    Ini.Val$Beta.NB <- rep(0, ncol(mf_nb))
  } else {
    Ini.Val$Beta.NB <- init.Beta.NB
  }

  if (is.null(init.Beta.PL) | length(init.Beta.PL) != ncol(mf_pareto)) {
    Ini.Val$Beta.PL <- rep(0, ncol(mf_pareto))
  } else {
    Ini.Val$Beta.PL <- init.Beta.PL
  }

  Ini.Val$Alpha.NB <- init.Alpha.NB
  Ini.Val$C <- init.C

  if (verbose) {
    object <- plinfl.nb.regression.fun(OBS.Y, OBS.X.obj, Ini.Val, Control)
  } else {
    capture.output(
      object <- plinfl.nb.regression.fun(OBS.Y, OBS.X.obj, Ini.Val, Control)
    )
  }
  object$par.mat$Beta.multinom.ZC <- as.numeric(object$par.mat$Beta.multinom.ZC)
  object$par.mat$Beta.multinom.PL <- as.numeric(object$par.mat$Beta.multinom.PL)
  object$par.mat$Beta.NB <- as.numeric(object$par.mat$Beta.NB)
  object$par.mat$Beta.PL <- as.numeric(object$par.mat$Beta.PL)

  names(object$par.mat$Beta.NB) <- c('(Intercept)', all.vars(formula_nb)[-1])
  #names(object$par.mat$Beta.multinom.ZC) <- c('(Intercept)',all.vars(formula_zi)[-1])
  names(object$par.mat$Beta.multinom.PL) <- c(
    '(Intercept)',
    all.vars(formula_evi)[-1]
  )
  names(object$par.mat$Beta.PL) <- c(
    '(Intercept)',
    all.vars(formula_pareto)[-1]
  )

  object$formulas <- list(
    formula_nb = formula_nb,
    #formula_zi = formula_zi,
    formula_evi = formula_evi,
    formula_pareto = formula_pareto
  )
  object$data <- list()

  object$data$data <- data %>%
    dplyr::select(dplyr::all_of(unique(c(
      all.vars(formula_nb),
      #all.vars(formula_zi),
      all.vars(formula_evi),
      all.vars(formula_pareto)
    )))) %>%
    na.omit()

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
  colnames(object$props) <- colnames(object$resp) <- c('count', 'pareto')
  object$par.mat$Props <- NULL
  object$coef <- object$par.mat
  object$par.mat <- NULL
  object$coef$Beta.multinom.ZC <- NULL
  object$par_all <- object$par_all[
    length(Ini.Val$Beta.multinom.ZC) + 1:length(object$par_all)
  ]

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
  object$fitted$prob_pareto <- object$props[, 2]
  object$fitted$posterior_count <- object$resp[, 1]
  object$fitted$posterior_pareto <- object$resp[, 2]

  class(object) <- 'evinb'
  return(object)
}


#' A single bootstrap run on an ezinb object
#'
#' @param object The evinb object to run the bootstrap on
#' @param block Optional string specifying varible for block bootstrapping
#' @param timing Should time be kept
#' @param track_progress Should progress be tracked (experimental)
#' @param id Id when tracking progress
#' @param maxboot Number of bootstraps when tracking progress
#'
#' @return A bootstrapped evinb object
#'
#' @noRd
bootrun_evinb <- function(
  object,
  block = NULL,
  timing = TRUE,
  track_progress = FALSE,
  id = NULL,
  maxboot = NULL
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
  OBS.X.obj$X.multinom.ZC <- object$data$x.nb[boot_id, ] # note: only for functionality, does not update
  OBS.X.obj$X.multinom.PL <- object$data$x.multinom.pl[boot_id, ]
  OBS.X.obj$X.NB <- object$data$x.nb[boot_id, ]
  OBS.X.obj$X.PL <- object$data$x.pl[boot_id, ]
  Control <- object$control

  Ini.Val <- list()
  Ini.Val$Beta.multinom.ZC <- c(-100, rep(0, ncol(OBS.X.obj$X.multinom.ZC)))
  Ini.Val$Beta.multinom.PL <- as.numeric(object$coef$Beta.multinom.PL)
  Ini.Val$Beta.NB <- as.numeric(object$coef$Beta.NB)
  Ini.Val$Beta.PL <- as.numeric(object$coef$Beta.PL)
  Ini.Val$Alpha.NB <- object$coef$Alpha.NB
  Ini.Val$C <- object$coef$C
  capture.output(
    evinb_boot <- plinfl.nb.regression.fun(OBS.Y, OBS.X.obj, Ini.Val, Control)
  )

  evinb_boot$par.mat$Beta.multinom.PL <- as.numeric(
    evinb_boot$par.mat$Beta.multinom.PL
  )
  evinb_boot$par.mat$Beta.NB <- as.numeric(evinb_boot$par.mat$Beta.NB)
  evinb_boot$par.mat$Beta.PL <- as.numeric(evinb_boot$par.mat$Beta.PL)

  names(evinb_boot$par.mat$Beta.NB) <- c(
    '(Intercept)',
    all.vars(object$formulas$formula_nb)[-1]
  )
  #names(evinb_boot$par.mat$Beta.multinom.ZC) <- c('(Intercept)',all.vars(object$formulas$formula_zi)[-1])
  names(evinb_boot$par.mat$Beta.multinom.PL) <- c(
    '(Intercept)',
    all.vars(object$formulas$formula_evi)[-1]
  )
  names(evinb_boot$par.mat$Beta.PL) <- c(
    '(Intercept)',
    all.vars(object$formulas$formula_pareto)[-1]
  )

  evinb_boot$resp <- NULL
  evinb_boot$props <- evinb_boot$par.mat$Props[, 2:3]

  evinb_boot$par.mat$Props <- NULL

  evinb_boot$coef <- evinb_boot$par.mat
  evinb_boot$par.mat <- NULL
  evinb_boot$coef$Beta.multinom.ZC <- NULL
  evinb_boot$formulas <- object$formulas

  evinb_boot$y.hat.plexpElogy <- NULL
  evinb_boot$y.hat.pl.E.inv.y <- NULL
  evinb_boot$y.hat.plmedian <- NULL
  evinb_boot$y.hat.plmean <- NULL
  evinb_boot$mu.nb.vec <- NULL
  evinb_boot$alpha.pl.vec <- NULL
  evinb_boot$exp.E.log.y <- NULL
  evinb_boot$median.pl.vec <- NULL
  evinb_boot$mean.pl.vec <- NULL

  evinb_boot$data <- NULL
  evinb_boot$boot_id <- boot_id

  class(evinb_boot) <- 'evinb_boot'
  if (timing) {
    evinb_boot$time <- difftime(Sys.time(), tim, units = 'secs')
  }
  if (track_progress) {
    cat(
      "\n ======= Bootstrap ",
      id,
      " of ",
      maxboot,
      "done in ",
      round(evinb_boot$time, 1),
      "seconds =====  Converged:",
      evinb_boot$converge,
      '. Number of em steps:',
      length(evinb_boot$log.lik.vec.all),
      "===== \n"
    )
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
#' @param block Optional string indicating a case-identifier variable when using block bootstrapping
#' @param boot_seed Optional bootstrap seed to ensure reproducible results.
#' @param max.diff.par Tolerance for EM algorithm. Will be considered to have converged if the maximum absolute difference in the parameter estimates are lower than this value
#' @param max.no.em.steps Maximum number of EM steps to run. Will be considered to not have converged if this number is reached and convergence is not reached
#' @param max.no.em.steps.warmup Number of EM steps in the warmup rounds
#' @param c.lim Integer range defining the possible values of C
#' @param prune.c.range The pruning fraction for the C range. Useful if there are many unique values in c.lim and the model is slow to fit. If FALSE, no pruning is done. If numeric [0-1], the C range is pruned to approximately length(c.lim)*(1-prune.c.range)
#' @param max.upd.par.pl.multinomial Maximum parameter change step size in the extreme value inflation component
#' @param max.upd.par.nb Maximum parameter change step size in the count component
#' @param max.upd.par.pl Maximum parameter change step size in the pareto component
#' @param no.m.bfgs.steps.multinomial Number of BFGS steps for the multinomial model
#' @param no.m.bfgs.steps.nb Number of BFGS steps for the negative binomial model
#' @param no.m.bfgs.steps.pl Number of BFGS steps for the pareto model
#' @param pdf.pl.type Probability density function type for the pareto component. Either 'approx' or 'exact'. 'approx' is adviced in most cases
#' @param eta.int Initial values for eta
#' @param init.Beta.multinom.PL Initial values for beta parameters in the extreme value inflation component. Vector of same length as number of parameters in the extreme value inflation component or NULL (which gives starting values of 0)
#' @param init.Beta.NB Initial values for beta parameters in the count component. Vector of same length as number of parameters in the count component or NULL (which gives starting values of 0)
#' @param init.Beta.PL Initial values for beta parameters in the pareto component. Vector of same length as number of parameters in the pareto component or NULL (which gives starting values of 0)
#' @param init.Alpha.NB Initial value of Alpha NB, integer or NULL (giving a starting value of 0)
#' @param init.C Initial value of C. Integer which should be within the C_lim range.
#' @param verbose Should progress be printed for the first run of evinb
#'
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#'
#' @return An object of class 'evinb'
#' @export
#'
#' @examples
#' data(genevzinb2)
#' model <- evinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 10)
#'
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
  max.diff.par = 1e-2,
  max.no.em.steps = 500,
  max.no.em.steps.warmup = 5,
  c.lim = c(50, 1000),
  prune.c.range = FALSE,
  max.upd.par.pl.multinomial = 0.5,
  max.upd.par.nb = 0.5,
  max.upd.par.pl = 0.5,
  no.m.bfgs.steps.multinomial = 3,
  no.m.bfgs.steps.nb = 3,
  no.m.bfgs.steps.pl = 3,
  pdf.pl.type = "approx",
  eta.int = c(-1, 1),
  init.Beta.multinom.PL = NULL,
  init.Beta.NB = NULL,
  init.Beta.PL = NULL,
  init.Alpha.NB = 0.01,
  init.C = 200,
  verbose = FALSE
) {
  i <- 'temp_iter'

  if (is.null(formula_evi)) {
    formula_evi <- formula_nb
  }
  if (is.null(formula_pareto)) {
    formula_pareto <- formula_nb
  }
  t1 <- Sys.time()
  full_run <- run_evinb(
    formula_nb = formula_nb,
    formula_evi = formula_evi,
    formula_pareto = formula_pareto,
    data = data,
    max.diff.par = max.diff.par,
    max.no.em.steps = max.no.em.steps,
    max.no.em.steps.warmup = max.no.em.steps.warmup,
    c.lim = c.lim,
    prune.c.range = prune.c.range,
    max.upd.par.pl.multinomial = max.upd.par.pl.multinomial,
    max.upd.par.nb = max.upd.par.nb,
    max.upd.par.pl = max.upd.par.pl,
    no.m.bfgs.steps.multinomial = no.m.bfgs.steps.multinomial,
    no.m.bfgs.steps.nb = no.m.bfgs.steps.nb,
    no.m.bfgs.steps.pl = no.m.bfgs.steps.pl,
    pdf.pl.type = pdf.pl.type,
    eta.int = eta.int,
    init.Beta.multinom.PL = init.Beta.multinom.PL,
    init.Beta.NB = init.Beta.NB,
    init.Beta.PL = init.Beta.PL,
    init.Alpha.NB = init.Alpha.NB,
    init.C = init.C,
    verbose = verbose
  )
  runtime <- difftime(Sys.time(), t1)

  full_run$block <- block
  if (!is.null(block)) {
    if (is.character(block)) {
      block2 <- data %>%
        dplyr::select(dplyr::all_of(unique(c(
          all.vars(formula_nb),
          #all.vars(formula_zi),
          all.vars(formula_evi),
          all.vars(formula_pareto),
          block
        )))) %>%
        na.omit() %>%
        dplyr::select(dplyr::all_of(block)) %>%
        dplyr::pull()

      full_run$data$data <- dplyr::bind_cols(
        data %>%
          dplyr::select(dplyr::all_of(unique(c(
            all.vars(formula_nb),
            #all.vars(formula_zi),
            all.vars(formula_evi),
            all.vars(formula_pareto),
            block
          )))) %>%
          na.omit() %>%
          dplyr::select(dplyr::all_of(block)),
        full_run$data$data
      )
    }
  } else {
    block2 <- NULL
  }

  if (bootstrap) {
    if (multicore) {
      if (is.null(ncores)) {
        ncores <- parallel::detectCores() - 1
      }
      if (.Platform$OS.type == 'windows') {
        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
      } else {
        doParallel::registerDoParallel(cores = ncores)
      }
      ex_time <- runtime * n_bootstraps / ncores
    } else {
      ex_time <- runtime * n_bootstraps
    }
    cat(
      "\n ======",
      "Approximate runtime for bootstraps is",
      ex_time,
      attributes(runtime)$units,
      ". Note: This is a very rough estimate of the runtime."
    )

    if (verbose) {
      boots <- foreach::foreach(
        i = 1:n_bootstraps,
        .options.RNG = boot_seed,
        .packages = 'evinf'
      ) %dorng%
        try(bootrun_evinb(
          full_run,
          block2,
          track_progress = T,
          id = i,
          maxboot = n_bootstraps
        ))
    } else {
      boots <- foreach::foreach(
        i = 1:n_bootstraps,
        .options.RNG = boot_seed,
        .packages = 'evinf'
      ) %dorng%
        try(bootrun_evinb(
          full_run,
          block2,
          track_progress = F,
          id = i,
          maxboot = n_bootstraps
        ))
    }

    out <- c(full_run, list(bootstraps = boots))

    if (multicore) {
      if (.Platform$OS.type == 'windows') {
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
      } else {
        doParallel::stopImplicitCluster()
      }
    }
    names(boots) <- paste('bootstrap_', 1:length(boots), sep = "")
  } else {
    out <- full_run
  }

  class(out) <- 'evinb'
  return(out)
}
