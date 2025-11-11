#' Running an extreme value and zero inflated negative binomial model
#'
#' @param formula_nb Formula for the negative binomial (count) component of the model.
#' @param formula_zi Formula for the zero-inflation component of the model. If NULL taken as the same formula as nb
#' @param formula_evi Formula for the extreme-value inflation component of the model. If NULL taken as the same formula as nb
#' @param formula_pareto Formula for the pareto (extreme value) component of the model. If NULL taken as the same formula as nb
#' @param data Data to run the model on
#' @param max.diff.par Tolerance for EM algorithm. Will be considered to have converged if the maximum absolute difference in the parameter estimates are lower than this value
#' @param max.no.em.steps Maximum number of EM steps to run. Will be considered to not have converged if this number is reached and convergence is not reached
#' @param max.no.em.steps.warmup Number of EM steps in the warmup rounds
#' @param c.lim Integer range defining the possible values of C
#' @param prune.c.range The pruning fraction for the C range. Useful if there are many unique values in c.lim and the model is slow to fit. If FALSE, no pruning is done. If numeric [0-1], the C range is pruned to approximately length(c.lim)*(1-prune.c.range)
#' @param max.upd.par.zc.multinomial Maximum parameter change step size in the zero inflation component
#' @param max.upd.par.pl.multinomial Maximum parameter change step size in the extreme value inflation component
#' @param max.upd.par.nb Maximum parameter change step size in the count component
#' @param max.upd.par.pl Maximum parameter change step size in the pareto component
#' @param no.m.bfgs.steps.multinomial Number of BFGS steps for the multinomial model
#' @param no.m.bfgs.steps.nb Number of BFGS steps for the negative binomial model
#' @param no.m.bfgs.steps.pl Number of BFGS steps for the pareto model
#' @param pdf.pl.type Probability density function type for the pareto component. Either 'approx' or 'exact'. 'approx' is adviced in most cases
#' @param eta.int Initial values for eta
#' @param init.Beta.multinom.ZC Initial values for beta parameters in the zero value inflation component. Vector of same length as number of parameters in the zero value inflation component or NULL (which gives starting values of 0)
#' @param init.Beta.multinom.PL Initial values for beta parameters in the extreme value inflation component. Vector of same length as number of parameters in the extreme value inflation component or NULL (which gives starting values of 0)
#' @param init.Beta.NB Initial values for beta parameters in the count component. Vector of same length as number of parameters in the count component or NULL (which gives starting values of 0)
#' @param init.Beta.PL Initial values for beta parameters in the pareto component. Vector of same length as number of parameters in the pareto component or NULL (which gives starting values of 0)
#' @param init.Alpha.NB Initial value of Alpha NB, integer or NULL (giving a starting value of 0)
#' @param init.C Initial value of C. Integer which should be within the C_lim range.
#' @param verbose Should progress be printed for the first run of evzinb
#'
#' @return An object of class 'evzinb'
#' @noRd
run_evzinb <- function(
  formula_nb,
  formula_zi = NULL,
  formula_evi = NULL,
  formula_pareto = NULL,
  data,
  max.diff.par = 1e-3,
  max.no.em.steps = 200,
  max.no.em.steps.warmup = 5,
  c.lim = c(50, 1000),
  prune.c.range = FALSE,
  max.upd.par.zc.multinomial = 0.5,
  max.upd.par.pl.multinomial = 0.5,
  max.upd.par.nb = 0.5,
  max.upd.par.pl = 0.5,
  no.m.bfgs.steps.multinomial = 3,
  no.m.bfgs.steps.nb = 3,
  no.m.bfgs.steps.pl = 3,
  pdf.pl.type = "approx",
  eta.int = c(-1, 1),
  init.Beta.multinom.ZC = NULL,
  init.Beta.multinom.PL = NULL,
  init.Beta.NB = NULL,
  init.Beta.PL = NULL,
  init.Alpha.NB = 0.01,
  init.C = 200,
  verbose = TRUE
) {
  if (is.null(formula_evi)) {
    formula_evi <- formula_nb
  }
  if (is.null(formula_pareto)) {
    formula_pareto <- formula_nb
  }
  if (is.null(formula_zi)) {
    formula_zi <- formula_nb
  }

  mf_nb <- model.frame(formula_nb, data)
  mf_zi <- model.frame(formula_zi, data)
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
    max.upd.par.zc.multinomial = max.upd.par.zc.multinomial,
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
  if (
    is.null(init.Beta.multinom.ZC) |
      length(init.Beta.multinom.ZC) != ncol(mf_zi)
  ) {
    Ini.Val$Beta.multinom.ZC <- rep(0, ncol(mf_zi))
  } else {
    Ini.Val$Beta.multinom.ZC <- init.Beta.multinom.ZC
  }

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
    object <- zerinfl.nb.pl.regression.fun(OBS.Y, OBS.X.obj, Ini.Val, Control)
  } else {
    capture.output(
      object <- zerinfl.nb.pl.regression.fun(OBS.Y, OBS.X.obj, Ini.Val, Control)
    )
  }
  object$par.mat$Beta.multinom.ZC <- as.numeric(object$par.mat$Beta.multinom.ZC)
  object$par.mat$Beta.multinom.PL <- as.numeric(object$par.mat$Beta.multinom.PL)
  object$par.mat$Beta.NB <- as.numeric(object$par.mat$Beta.NB)
  object$par.mat$Beta.PL <- as.numeric(object$par.mat$Beta.PL)

  names(object$par.mat$Beta.NB) <- c('(Intercept)', all.vars(formula_nb)[-1])
  names(object$par.mat$Beta.multinom.ZC) <- c(
    '(Intercept)',
    all.vars(formula_zi)[-1]
  )
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
    formula_zi = formula_zi,
    formula_evi = formula_evi,
    formula_pareto = formula_pareto
  )
  object$data <- list()
  full_data <- data %>%
    dplyr::select(dplyr::all_of(unique(c(
      all.vars(formula_nb),
      all.vars(formula_zi),
      all.vars(formula_evi),
      all.vars(formula_pareto)
    )))) %>%
    na.omit()

  object$data$data <- full_data
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
    'pareto'
  )

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
  object$fitted$prob_zero <- object$props[, 1]
  object$fitted$prob_count <- object$props[, 2]
  object$fitted$prob_pareto <- object$props[, 3]
  object$fitted$posterior_zero <- object$resp[, 1]
  object$fitted$posterior_count <- object$resp[, 2]
  object$fitted$posterior_pareto <- object$resp[, 3]

  class(object) <- 'evzinb'
  return(object)
}


#' Running an extreme value and zero inflated negative binomial model with bootstrapping
#'
#' @param formula_nb Formula for the negative binomial (count) component of the model
#' @param formula_zi Formula for the zero-inflation component of the model. If NULL taken as the same formula as nb
#' @param formula_evi Formula for the extreme-value inflation component of the model. If NULL taken as the same formula as nb
#' @param formula_pareto Formula for the pareto (extreme value) component of the model. If NULL taken as the same formula as nb
#' @param data data to run the model on
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
#' @param prune.c.range The pruning fraction for the C range. Useful if there are many unique values in c.lim and the model is slow to fit. If FALSE, no pruning is done. If numeric [0-1], the C range is pruned to approximately length(c.lim)*prune.c.range
#' @param max.upd.par.zc.multinomial Maximum parameter change step size in the zero inflation component
#' @param max.upd.par.pl.multinomial Maximum parameter change step size in the extreme value inflation component
#' @param max.upd.par.nb Maximum parameter change step size in the count component
#' @param max.upd.par.pl Maximum parameter change step size in the pareto component
#' @param no.m.bfgs.steps.multinomial Number of BFGS steps for the multinomial model
#' @param no.m.bfgs.steps.nb Number of BFGS steps for the negative binomial model
#' @param no.m.bfgs.steps.pl Number of BFGS steps for the pareto model
#' @param pdf.pl.type Probability density function type for the pareto component. Either 'approx' or 'exact'. 'approx' is adviced in most cases
#' @param eta.int Initial values for eta
#' @param init.Beta.multinom.ZC Initial values for beta parameters in the zero value inflation component. Vector of same length as number of parameters in the zero value inflation component or NULL (which gives starting values of 0)
#' @param init.Beta.multinom.PL Initial values for beta parameters in the extreme value inflation component. Vector of same length as number of parameters in the extreme value inflation component or NULL (which gives starting values of 0)
#' @param init.Beta.NB Initial values for beta parameters in the count component. Vector of same length as number of parameters in the count component or NULL (which gives starting values of 0)
#' @param init.Beta.PL Initial values for beta parameters in the pareto component. Vector of same length as number of parameters in the pareto component or NULL (which gives starting values of 0)
#' @param init.Alpha.NB Initial value of Alpha NB, integer or NULL (giving a starting value of 0)
#' @param init.C Initial value of C. Integer which should be within the C_lim range.
#' @param verbose Logical: should progress of the full run of the model be tracked?
#'
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
#'
#' @return An object of class 'evzinb'
#' @export
#'
#' @examples
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 10)
#'
evzinb <- function(
  formula_nb,
  formula_zi = NULL,
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
  max.upd.par.zc.multinomial = 0.5,
  max.upd.par.pl.multinomial = 0.5,
  max.upd.par.nb = 0.5,
  max.upd.par.pl = 0.5,
  no.m.bfgs.steps.multinomial = 3,
  no.m.bfgs.steps.nb = 3,
  no.m.bfgs.steps.pl = 3,
  pdf.pl.type = "approx",
  eta.int = c(-1, 1),
  init.Beta.multinom.ZC = NULL,
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
  if (is.null(formula_zi)) {
    formula_zi <- formula_nb
  }
  t1 <- Sys.time()
  full_run <- run_evzinb(
    formula_nb = formula_nb,
    formula_zi = formula_zi,
    formula_evi = formula_evi,
    formula_pareto = formula_pareto,
    data = data,
    max.diff.par = max.diff.par,
    max.no.em.steps = max.no.em.steps,
    max.no.em.steps.warmup = max.no.em.steps.warmup,
    c.lim = c.lim,
    prune.c.range = prune.c.range,
    max.upd.par.zc.multinomial = max.upd.par.zc.multinomial,
    max.upd.par.pl.multinomial = max.upd.par.pl.multinomial,
    max.upd.par.nb = max.upd.par.nb,
    max.upd.par.pl = max.upd.par.pl,
    no.m.bfgs.steps.multinomial = no.m.bfgs.steps.multinomial,
    no.m.bfgs.steps.nb = no.m.bfgs.steps.nb,
    no.m.bfgs.steps.pl = no.m.bfgs.steps.pl,
    pdf.pl.type = pdf.pl.type,
    eta.int = eta.int,
    init.Beta.multinom.ZC = init.Beta.multinom.ZC,
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
          all.vars(formula_zi),
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
            all.vars(formula_zi),
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
        parallel::clusterExport(
          cl,
          c(
            'bootrun_evzinb',
            "zerinfl.nb.pl.regression.fun",
            "zerinfl.nb.pl.reg.cond.c.fun",
            "log_lik_fun",
            "_evinf_log_lik_fun"
          ),
          envir = environment(bootrun_evzinb)
        )
        parallel::clusterExport(
          cl,
          c('full_run', 'n_bootstraps'),
          envir = environment()
        )
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
        try(bootrun_evzinb(
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
        try(bootrun_evzinb(
          full_run,
          block2,
          track_progress = F,
          id = i,
          maxboot = n_bootstraps
        ))
    }
    if (multicore) {
      if (.Platform$OS.type == 'windows') {
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
      } else {
        doParallel::stopImplicitCluster()
      }
    }
    names(boots) <- paste('bootstrap_', 1:length(boots), sep = "")
    out <- c(full_run, list(bootstraps = boots))
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
#' @param timing Should time be kept
#' @param track_progress Should progress be tracked (experimental)
#' @param id Id when tracking progress
#' @param maxboot Number of bootstraps when tracking progress
#'
#' @return A bootstrapped evzinb object
#'
#' @noRd
bootrun_evzinb <- function(
  object,
  block = NULL,
  timing = TRUE,
  track_progress = TRUE,
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
  OBS.X.obj$X.multinom.ZC <- object$data$x.multinom.zc[boot_id, ]
  OBS.X.obj$X.multinom.PL <- object$data$x.multinom.pl[boot_id, ]
  OBS.X.obj$X.NB <- object$data$x.nb[boot_id, ]
  OBS.X.obj$X.PL <- object$data$x.pl[boot_id, ]
  Control <- object$control

  Ini.Val <- list()
  Ini.Val$Beta.multinom.ZC <- as.numeric(object$coef$Beta.multinom.ZC)
  Ini.Val$Beta.multinom.PL <- as.numeric(object$coef$Beta.multinom.PL)
  Ini.Val$Beta.NB <- as.numeric(object$coef$Beta.NB)
  Ini.Val$Beta.PL <- as.numeric(object$coef$Beta.PL)
  Ini.Val$Alpha.NB <- object$coef$Alpha.NB
  Ini.Val$C <- object$coef$C
  capture.output(
    evzinb_boot <- zerinfl.nb.pl.regression.fun(
      OBS.Y,
      OBS.X.obj,
      Ini.Val,
      Control
    )
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
    all.vars(object$formulas$formula_nb)[-1]
  )
  names(evzinb_boot$par.mat$Beta.multinom.ZC) <- c(
    '(Intercept)',
    all.vars(object$formulas$formula_zi)[-1]
  )
  names(evzinb_boot$par.mat$Beta.multinom.PL) <- c(
    '(Intercept)',
    all.vars(object$formulas$formula_evi)[-1]
  )
  names(evzinb_boot$par.mat$Beta.PL) <- c(
    '(Intercept)',
    all.vars(object$formulas$formula_pareto)[-1]
  )

  evzinb_boot$props <- evzinb_boot$par.mat$Props

  evzinb_boot$par.mat$Props <- NULL
  evzinb_boot$coef <- evzinb_boot$par.mat
  evzinb_boot$par.mat <- NULL
  evzinb_boot$resp <- NULL

  evzinb_boot$formulas <- object$formulas

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

  evzinb_boot$data <- NULL
  evzinb_boot$boot_id <- boot_id
  if (timing) {
    evzinb_boot$time <- difftime(Sys.time(), tim, units = 'secs')
  }
  if (track_progress) {
    cat(
      "\n ======= Bootstrap ",
      id,
      " of ",
      maxboot,
      "done in ",
      round(evzinb_boot$time, 1),
      "seconds. Converged:",
      evzinb_boot$converge,
      '. Number of em steps:',
      length(evzinb_boot$log.lik.vec.all),
      "===== \n"
    )
  }
  return(evzinb_boot)
}
