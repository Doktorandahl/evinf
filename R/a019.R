## Internal functions for the EVZINB and EVINB functions
zerinfl.nb.pl.reg.cond.c.fun <- function(y, x.obj, ini.val, control) {
  #Initialize parameters and data
  x.multinom.zc <- x.obj$X.multinom.ZC
  x.multinom.pl <- x.obj$X.multinom.PL

  x.nb <- x.obj$X.NB
  x.pl <- x.obj$X.PL

  n <- max(nrow(x.nb), nrow(x.pl), nrow(x.multinom.zc), nrow(x.multinom.pl))

  if (is.null(dim(x.multinom.zc))) {
    x.multinom.zc.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.multinom.zc.extended <- cbind(1, x.multinom.zc)
  }
  if (is.null(dim(x.multinom.pl))) {
    x.multinom.pl.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.multinom.pl.extended <- cbind(1, x.multinom.pl)
  }
  if (is.null(dim(x.nb))) {
    x.nb.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.nb.extended <- cbind(1, x.nb)
  }
  if (is.null(dim(x.pl))) {
    x.pl.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.pl.extended <- cbind(1, x.pl)
  }

  n.beta.zc.multinomial <- dim(x.multinom.zc.extended)[2]
  n.beta.pl.multinomial <- dim(x.multinom.pl.extended)[2]
  n.theta.multinomial <- n.beta.zc.multinomial + n.beta.pl.multinomial
  n.beta.nb <- dim(x.nb.extended)[2]
  n.beta.pl <- dim(x.pl.extended)[2]

  beta.zc.multinomial.old <- ini.val$Beta.multinom.ZC
  beta.pl.multinomial.old <- ini.val$Beta.multinom.PL
  beta.nb.old <- ini.val$Beta.NB
  alpha.nb.old <- ini.val$Alpha.NB
  theta.nb.old <- c(beta.nb.old, alpha.nb.old)
  beta.pl.old <- ini.val$Beta.PL
  c.pl <- ini.val$C

  #Initial log likelihood value

  #func.val.initial <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)
  func.val.initial <- log_lik_fun(
    beta.zc.multinomial.old,
    beta.pl.multinomial.old,
    beta.nb.old,
    alpha.nb.old,
    beta.pl.old,
    c.pl,
    x.multinom.zc.extended,
    x.multinom.pl.extended,
    x.nb.extended,
    x.pl.extended,
    y
  )

  # The maximum change after one EM step needs to be small in order for the algorithm to stop
  max.abs.par.diff <- 100

  #Number of nas produced
  na.beta.nb <- 0
  na.alpha.nb <- 0
  na.beta.pl <- 0
  na.beta.mult.zc <- 0
  na.beta.mult.pl <- 0

  ################
  #Start EM algorithm
  ###################
  i.em <- 1
  func.val.vec <- c()
  func.val.old <- -1e50
  max.no.em.steps <- control$max.no.em.steps
  max.diff.par <- control$max.diff.par
  max.upd.par <- control$max.upd.par.nb
  no.m.bfgs.steps <- control$no.m.bfgs.steps.nb

  #max.no.em.steps <- 21

  while (i.em < max.no.em.steps & max.abs.par.diff > max.diff.par) {
    beta.zc.multinomial.start.em <- beta.zc.multinomial.old
    beta.pl.multinomial.start.em <- beta.pl.multinomial.old
    beta.nb.start.em <- beta.nb.old
    alpha.nb.start.em <- abs(alpha.nb.old)
    beta.pl.start.em <- beta.pl.old

    par.start.em <- c(
      beta.zc.multinomial.start.em,
      beta.pl.multinomial.start.em,
      beta.nb.start.em,
      alpha.nb.start.em,
      beta.pl.start.em
    )

    #Update parameters with BFGS
    upd.obj <- update_bfgs_fun(
      beta.zc.multinomial.old,
      beta.pl.multinomial.old,
      beta.nb.old,
      alpha.nb.old,
      beta.pl.old,
      c.pl,
      x.multinom.zc.extended,
      x.multinom.pl.extended,
      x.nb.extended,
      x.pl.extended,
      y,
      max.upd.par,
      control$no.m.bfgs.steps.nb
    )

    if (sum(is.na(upd.obj$beta_nb_old)) == 0) {
      beta.nb.new <- upd.obj$beta_nb_old
    } else {
      beta.nb.new <- beta.nb.old
      na.beta.nb <- na.beta.nb + length(is.na(upd.obj$beta_nb_old))
    }

    if (sum(is.na(upd.obj$alpha_nb_old)) == 0) {
      alpha.nb.new <- upd.obj$alpha_nb_old
    } else {
      alpha.nb.new <- alpha.nb.old
      na.alpha.nb <- na.alpha.nb + length(is.na(upd.obj$alpha_nb_old))
    }

    if (sum(is.na(upd.obj$beta_pl_old)) == 0) {
      beta.pl.new <- upd.obj$beta_pl_old
    } else {
      beta.pl.new <- beta.pl.old
      na.beta.pl <- na.beta.pl + length(is.na(upd.obj$beta_pl_old))
    }

    if (sum(is.na(upd.obj$gamma_z_old)) == 0) {
      beta.zc.multinomial.new <- upd.obj$gamma_z_old
    } else {
      beta.zc.multinomial.new <- beta.zc.multinomial.old
      na.beta.mult.zc <- na.beta.mult.zc + length(is.na(upd.obj$gamma_z_old))
    }

    if (sum(is.na(upd.obj$gamma_pl_old)) == 0) {
      beta.pl.multinomial.new <- upd.obj$gamma_pl_old
    } else {
      beta.pl.multinomial.new <- beta.pl.multinomial.old
      na.beta.mult.pl <- na.beta.mult.pl + length(is.na(upd.obj$gamma_pl_old))
    }

    # alpha.nb.new <- upd.obj$alpha_nb_old
    # beta.pl.new <- upd.obj$beta_pl_old
    # beta.zc.multinomial.new <- upd.obj$gamma_z_old
    # beta.pl.multinomial.new <- upd.obj$gamma_pl_old

    par.end.em <- c(
      beta.zc.multinomial.new,
      beta.pl.multinomial.new,
      beta.nb.new,
      alpha.nb.new,
      beta.pl.new
    )
    par.diff.em <- par.end.em - par.start.em

    ##############Make sure the likelihood increases

    ################################################
    ##################Optimize theta.nb.old
    nb.log.lik.optim <- function(eta) {
      #Go back eta times the step that was already made
      theta.nb.old <- c(beta.nb.old, alpha.nb.old) +
        eta * upd.obj$change_nb_bfgs # change.nb.bfgs
      beta.nb.tmp <- theta.nb.old[1:n.beta.nb]
      alpha.nb.tmp <- theta.nb.old[n.beta.nb + 1]
      return(
        -1.0 *
          log_lik_fun(
            beta.zc.multinomial.old,
            beta.pl.multinomial.old,
            beta.nb.tmp,
            alpha.nb.tmp,
            beta.pl.old,
            c.pl,
            x.multinom.zc.extended,
            x.multinom.pl.extended,
            x.nb.extended,
            x.pl.extended,
            y
          )
      )
    }

    if (is.na(upd.obj$func_val_after_nb) == FALSE) {
      if (upd.obj$func_val_after_nb < upd.obj$func_val_before_bfgs) {
        if (sum(is.na(upd.obj$change_nb_bfgs)) == 0) {
          ###########################
          change.nb.obj <- optimise(
            f = nb.log.lik.optim,
            interval = control$eta.int
          )
          eta.nb <- change.nb.obj$minimum

          theta.nb.old <- c(beta.nb.old, alpha.nb.old)
          theta.nb.after.optim <- theta.nb.old + eta.nb * upd.obj$change_nb_bfgs
          beta.nb.after.optim <- theta.nb.after.optim[1:n.beta.nb]
          alpha.nb.after.optim <- theta.nb.after.optim[n.beta.nb + 1]

          func.val.after.nb.optim <- log_lik_fun(
            beta.zc.multinomial.old,
            beta.pl.multinomial.old,
            beta.nb.after.optim,
            alpha.nb.after.optim,
            beta.pl.old,
            c.pl,
            x.multinom.zc.extended,
            x.multinom.pl.extended,
            x.nb.extended,
            x.pl.extended,
            y
          )

          theta.nb.old <- theta.nb.after.optim
          beta.nb.new <- theta.nb.old[1:n.beta.nb]
          alpha.nb.new <- theta.nb.old[n.beta.nb + 1]
        }
      }
    }
    ##################################
    ##############################################

    ################################################
    ##################Optimize theta.pl.old
    pl.log.lik.optim <- function(eta) {
      #Go back eta times the step that was already made
      beta.pl.tmp <- beta.pl.old + eta * upd.obj$change_pl_bfgs
      return(
        -1.0 *
          log_lik_fun(
            beta.zc.multinomial.old,
            beta.pl.multinomial.old,
            beta.nb.old,
            alpha.nb.old,
            beta.pl.tmp,
            c.pl,
            x.multinom.zc.extended,
            x.multinom.pl.extended,
            x.nb.extended,
            x.pl.extended,
            y
          )
      )
    }

    if (is.na(upd.obj$func_val_after_pl) == FALSE) {
      if (upd.obj$func_val_after_pl < upd.obj$func_val_before_bfgs) {
        ###########################
        #if(det(d2Qdbeta2.pl)>0){
        change.pl.obj <- optimise(
          f = pl.log.lik.optim,
          interval = control$eta.int
        )
        eta.pl <- change.pl.obj$minimum

        beta.pl.after.optim <- beta.pl.old + eta.pl * upd.obj$change_pl_bfgs

        func.val.after.pl.optim <- log_lik_fun(
          beta.zc.multinomial.old,
          beta.pl.multinomial.old,
          beta.nb.old,
          alpha.nb.old,
          beta.pl.after.optim,
          c.pl,
          x.multinom.zc.extended,
          x.multinom.pl.extended,
          x.nb.extended,
          x.pl.extended,
          y
        )

        beta.pl.new <- beta.pl.after.optim
      }
    }

    ################################################
    ##################Optimize theta.pl.old
    zc.multinomial.log.lik.optim <- function(eta) {
      #Go back eta times the step that was already made
      beta.zc.multinomial.tmp <- beta.zc.multinomial.old +
        eta * upd.obj$change_mult_z_bfgs
      return(
        -1.0 *
          log_lik_fun(
            beta.zc.multinomial.tmp,
            beta.pl.multinomial.old,
            beta.nb.old,
            alpha.nb.old,
            beta.pl.old,
            c.pl,
            x.multinom.zc.extended,
            x.multinom.pl.extended,
            x.nb.extended,
            x.pl.extended,
            y
          )
      )
    }

    if (is.na(upd.obj$func_val_after_mult_z) == FALSE) {
      if (upd.obj$func_val_after_mult_z < upd.obj$func_val_before_bfgs) {
        ###########################
        change.zc.multinomial.obj <- optimise(
          f = zc.multinomial.log.lik.optim,
          interval = control$eta.int
        )
        eta.zc.multinomial <- change.zc.multinomial.obj$minimum

        beta.zc.multinomial.after.optim <- beta.zc.multinomial.old +
          eta.zc.multinomial * upd.obj$change_mult_z_bfgs

        func.val.after.zc.multinomial.optim <- log_lik_fun(
          beta.zc.multinomial.after.optim,
          beta.pl.multinomial.old,
          beta.nb.old,
          alpha.nb.old,
          beta.pl.old,
          c.pl,
          x.multinom.zc.extended,
          x.multinom.pl.extended,
          x.nb.extended,
          x.pl.extended,
          y
        )

        beta.zc.multinomial.new <- beta.zc.multinomial.after.optim
      }
    }

    ##################################
    ##############################################

    ################################################
    ##################Optimize theta.pl.old
    pl.multinomial.log.lik.optim <- function(eta) {
      #Go back eta times the step that was already made
      beta.pl.multinomial.tmp <- beta.pl.multinomial.old +
        eta * upd.obj$change_mult_pl_bfgs
      return(
        -1.0 *
          log_lik_fun(
            beta.zc.multinomial.old,
            beta.pl.multinomial.tmp,
            beta.nb.old,
            alpha.nb.old,
            beta.pl.old,
            c.pl,
            x.multinom.zc.extended,
            x.multinom.pl.extended,
            x.nb.extended,
            x.pl.extended,
            y
          )
      )
    }

    ###########################
    if (is.na(upd.obj$func_val_after_mult_pl) == FALSE) {
      if (upd.obj$func_val_after_mult_pl < upd.obj$func_val_before_bfgs) {
        change.pl.multinomial.obj <- optimise(
          f = pl.multinomial.log.lik.optim,
          interval = control$eta.int
        )
        eta.pl.multinomial <- change.pl.multinomial.obj$minimum

        beta.pl.multinomial.after.optim <- beta.pl.multinomial.old +
          eta.pl.multinomial * upd.obj$change_mult_pl_bfgs

        func.val.after.pl.multinomial.optim <- log_lik_fun(
          beta.zc.multinomial.old,
          beta.pl.multinomial.after.optim,
          beta.nb.old,
          alpha.nb.old,
          beta.pl.old,
          c.pl,
          x.multinom.zc.extended,
          x.multinom.pl.extended,
          x.nb.extended,
          x.pl.extended,
          y
        )

        beta.pl.multinomial.new <- beta.pl.multinomial.after.optim
      }
    }
    ##################################
    ##############################################

    ########################
    ############################## Finished updating parameters

    par.end.em <- c(
      beta.zc.multinomial.new,
      beta.pl.multinomial.new,
      beta.nb.new,
      alpha.nb.new,
      beta.pl.new
    )

    par.diff <- par.end.em - par.start.em
    max.abs.par.diff <- max(abs(par.diff))

    func.val.after.bfgs <- log_lik_fun(
      beta.zc.multinomial.new,
      beta.pl.multinomial.new,
      beta.nb.new,
      alpha.nb.new,
      beta.pl.new,
      c.pl,
      x.multinom.zc.extended,
      x.multinom.pl.extended,
      x.nb.extended,
      x.pl.extended,
      y
    )

    par.all.tmp <- c(
      beta.zc.multinomial.old,
      beta.pl.multinomial.old,
      beta.nb.old,
      alpha.nb.old,
      beta.pl.old
    )

    beta.nb.old <- beta.nb.new
    alpha.nb.old <- alpha.nb.new
    beta.pl.old <- beta.pl.new
    beta.zc.multinomial.old <- beta.zc.multinomial.new
    beta.pl.multinomial.old <- beta.pl.multinomial.new

    # if(func.val.after.bfgs<upd.obj$func_val_before_bfgs){
    #  print("The function decreased")
    # par.diff.mod <- 0.01*par.diff
    # par.all.tmp <- par.start.em + par.diff.mod
    # beta.nb.multinomial.old <- par.all.tmp[1:n.beta.nb.multinomial]
    # par.all.tmp <- par.all.tmp[-c(1:n.beta.nb.multinomial)]
    # beta.pl.multinomial.old <- par.all.tmp[1:n.beta.pl.multinomial]
    # par.all.tmp <- par.all.tmp[-c(1:n.beta.pl.multinomial)]
    # beta.nb.old <- par.all.tmp[1:n.beta.nb]
    # par.all.tmp <- par.all.tmp[-c(1:n.beta.nb)]
    # alpha.nb.old <- par.all.tmp[1]
    # par.all.tmp <- par.all.tmp[-1]
    # beta.pl.old <- par.all.tmp
    # rm(par.all.tmp)
    # }

    #func.val.old <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)

    #func.val.old <- func.val

    cat(
      'Iteration ',
      i.em,
      ': The max abs diff in parameters is ',
      max.abs.par.diff,
      '. The function value is ',
      round(func.val.after.bfgs, 4),
      '\n',
      sep = ''
    )

    func.val.vec[i.em] <- func.val.after.bfgs
    i.em <- i.em + 1
    # print(beta.nb.multinomial.old)
    # print(beta.pl.multinomial.old)
    # print(beta.nb.old)
    # print(beta.pl.old)
    # print(alpha.nb.old)
  }

  #Gather results after finishing the estimation
  res <- list()
  par.mat <- list()
  par.mat$Props <- upd.obj$prop
  par.mat$Beta.multinom.ZC <- beta.zc.multinomial.old
  par.mat$Beta.multinom.PL <- beta.pl.multinomial.old
  par.mat$Beta.NB <- beta.nb.old
  par.mat$Alpha.NB <- alpha.nb.old
  par.mat$Beta.PL <- as.numeric(beta.pl.old)
  par.mat$C <- c.pl

  res$par.mat <- par.mat

  res$par.all <- par.end.em

  res$resp <- upd.obj$resp

  res$log.lik.vec <- func.val.vec
  res$log.lik <- func.val.after.bfgs

  if (i.em < control$max.no.em.steps) {
    res$converge <- TRUE
  } else {
    res$converge <- FALSE
  }
  return(res)
}


# x.obj <- OBS.X.obj.r
# y <- OBS.Y.r
# control <- Control.r
# ini.val <- Ini.Val.r
# # ini.val <- TruePar

zerinfl.nb.pl.regression.fun <- function(y, x.obj, ini.val, control) {
  #Initialize data and parameters

  x.multinom.zc <- x.obj$X.multinom.ZC
  x.multinom.pl <- x.obj$X.multinom.PL
  x.nb <- x.obj$X.NB
  x.pl <- x.obj$X.PL

  n <- max(nrow(x.nb), nrow(x.pl), nrow(x.multinom.zc), nrow(x.multinom.pl))

  if (is.null(dim(x.multinom.zc))) {
    x.multinom.zc.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.multinom.zc.extended <- cbind(1, x.multinom.zc)
  }
  if (is.null(dim(x.multinom.pl))) {
    x.multinom.pl.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.multinom.pl.extended <- cbind(1, x.multinom.pl)
  }
  if (is.null(dim(x.nb))) {
    x.nb.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.nb.extended <- cbind(1, x.nb)
  }
  if (is.null(dim(x.pl))) {
    x.pl.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.pl.extended <- cbind(1, x.pl)
  }

  n.beta.zc.multinomial <- dim(x.multinom.zc.extended)[2]
  n.beta.pl.multinomial <- dim(x.multinom.pl.extended)[2]
  n.theta.multinomial <- n.beta.zc.multinomial + n.beta.pl.multinomial
  n.beta.nb <- dim(x.nb.extended)[2]
  n.beta.pl <- dim(x.pl.extended)[2]

  prel.val <- ini.val

  n <- max(nrow(x.nb), nrow(x.pl), nrow(x.multinom.zc), nrow(x.multinom.pl))

  #Determine which values of c to be investigated
  c.range <- unique(sort(y))[
    unique(sort(y)) >= control$c.lim[1] & unique(sort(y)) <= control$c.lim[2]
  ]

  ###When c is not well known we only run short versions of zerinfl.nb.pl.reg.cond.c.fun() in order to save time
  control.warmup <- control
  control.warmup$max.no.em.steps <- control$max.no.em.steps.warmup

  ###Initializing the warm-up phase when the update stops after control$max.no.em.steps.warmup EM steps
  #If the update of c is less than 1, stop the warmup
  c.abs.diff <- 100
  log.lik.vec.all <- NULL
  cat('Begin warm-up', '\n', sep = '')
  while (c.abs.diff > 0) {
    log.lik.vec <- c()
    est.obj <- zerinfl.nb.pl.reg.cond.c.fun(y, x.obj, prel.val, control.warmup)
    log.lik.vec.all <- c(log.lik.vec.all, est.obj$log.lik.vec)
    prel.val <- est.obj$par.mat
    props.old <- prel.val$Props
    beta.nb.old <- prel.val$Beta.NB
    alpha.nb.old <- prel.val$Alpha.NB
    beta.pl.old <- prel.val$Beta.PL
    beta.zc.multinomial.old <- prel.val$Beta.multinom.ZC
    beta.pl.multinomial.old <- prel.val$Beta.multinom.PL
    c.pl <- prel.val$C

    for (k in 1:length(c.range)) {
      log.lik.vec[k] <- log_lik_fun(
        beta.zc.multinomial.old,
        beta.pl.multinomial.old,
        beta.nb.old,
        alpha.nb.old,
        beta.pl.old,
        c.range[k],
        x.multinom.zc.extended,
        x.multinom.pl.extended,
        x.nb.extended,
        x.pl.extended,
        y
      )
    }

    c.pl.new <- c.range[which(log.lik.vec == max(log.lik.vec))]

    if (is.infinite(max(log.lik.vec))) {
      stop(
        'The log-likelihood is infinite for all tested values of c in the c-range. Try expanding the c-range with the c.lim argument of the function'
      )
    }

    c.abs.diff <- abs(c.pl.new - prel.val$C)

    prel.val$C <- c.pl.new
    func.val <- max(log.lik.vec)

    log.lik.vec.all <- c(log.lik.vec.all, func.val)

    cat(
      'The new c is ',
      c.pl.new,
      '. The function value is ',
      round(func.val, 4),
      '\n',
      sep = ''
    )
    prel.val$C <- c.pl.new
  }

  ###The end phase, when c is presumably rather accurate
  c.abs.diff <- 100
  cat('End warm-up. Run until convergence', '\n', sep = '')
  while (c.abs.diff > 0) {
    log.lik.vec <- c()

    est.obj <- zerinfl.nb.pl.reg.cond.c.fun(y, x.obj, prel.val, control)
    log.lik.vec.all <- c(log.lik.vec.all, est.obj$log.lik.vec)
    prel.val <- est.obj$par.mat
    props.old <- prel.val$Props
    beta.nb.old <- prel.val$Beta.NB
    alpha.nb.old <- prel.val$Alpha.NB
    beta.pl.old <- prel.val$Beta.PL
    beta.zc.multinomial.old <- est.obj$par.mat$Beta.multinom.ZC
    beta.pl.multinomial.old <- est.obj$par.mat$Beta.multinom.PL
    c.pl <- prel.val$C

    for (k in 1:length(c.range)) {
      log.lik.vec[k] <- log_lik_fun(
        beta.zc.multinomial.old,
        beta.pl.multinomial.old,
        beta.nb.old,
        alpha.nb.old,
        beta.pl.old,
        c.range[k],
        x.multinom.zc.extended,
        x.multinom.pl.extended,
        x.nb.extended,
        x.pl.extended,
        y
      )
    }

    c.pl.new <- c.range[which(log.lik.vec == max(log.lik.vec))]

    if (is.infinite(max(log.lik.vec))) {
      stop(
        'The log-likelihood is infinite for all tested values of c in the c-range. Try expanding the c-range with the c.lim argument of the function'
      )
    }

    c.abs.diff <- abs(c.pl.new - prel.val$C)

    prel.val$C <- c.pl.new
    func.val <- max(log.lik.vec.all)
    log.lik.vec.all <- c(log.lik.vec.all, func.val)

    cat(
      'The new c is ',
      c.pl.new,
      '. The function value is ',
      round(func.val, 4),
      '\n',
      sep = ''
    )
    prel.val$C <- c.pl.new
  }

  final.val <- prel.val

  #Mean conditional on negative binomial component
  #Pareto shape parameter
  mu.nb.vec <- c()
  alpha.pl.vec <- c()
  mean.pl.vec <- c()
  exp.E.log.y <- c() #exp(E(log(y))), a third alternative for prediction
  median.pl.vec <- c()
  E.inv.y <- c()
  y.hat.plmedian <- c()
  y.hat.plmean <- c()
  y.hat.plexpElogy <- c()
  y.hat.pl.E.inv.y <- c()
  for (i in 1:n) {
    mu.nb.vec[i] <- exp(x.nb.extended[i, ] %*% prel.val$Beta.NB)
    alpha.pl.vec[i] <- exp(x.pl.extended[i, ] %*% prel.val$Beta.PL)
    exp.E.log.y[i] <- c.pl.new * exp(1 / alpha.pl.vec[i])
    E.inv.y[i] <- alpha.pl.vec[i] / (c.pl.new * (alpha.pl.vec[i] + 1))
    if (alpha.pl.vec[i] > 1) {
      mean.pl.vec[i] <- alpha.pl.vec[i] * c.pl.new / (alpha.pl.vec[i] - 1)
    } else {
      mean.pl.vec[i] <- NA
    }
    median.pl.vec[i] <- median.pl.vec[i] <- c.pl.new * (2)^(1 / alpha.pl.vec[i])
    y.hat.plmedian[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] * median.pl.vec[i]
    y.hat.plmean[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] * mean.pl.vec[i]
    y.hat.plexpElogy[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] * exp.E.log.y[i]
    y.hat.pl.E.inv.y[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] / E.inv.y[i]
  }

  par.all <- c(
    final.val$Beta.multinom.ZC,
    final.val$Beta.multinom.PL,
    final.val$Beta.NB,
    final.val$Alpha.NB,
    final.val$Beta.PL,
    final.val$C
  )
  n.par <- length(par.all)
  BIC <- log(n) * n.par - 2 * func.val
  AIC <- 2 * n.par - 2 * func.val

  #Gather results from estimation
  out <- list()

  out$control <- control
  out$par.mat <- final.val
  out$log.lik.vec.all <- log.lik.vec.all
  out$log.lik <- log.lik.vec.all[length(log.lik.vec.all)]
  out$resp <- est.obj$resp
  out$converge <- est.obj$converge
  out$ini.val <- ini.val
  out$x.nb <- x.nb
  out$x.pl <- x.pl
  out$x.multinom.zc <- x.multinom.zc
  out$x.multinom.pl <- x.multinom.pl
  out$median.pl.vec <- median.pl.vec
  out$mean.pl.vec <- mean.pl.vec
  out$y <- y
  out$y.hat.plmedian <- y.hat.plmedian
  out$y.hat.plmean <- y.hat.plmean
  out$y.hat.plexpElogy <- y.hat.plexpElogy
  out$mu.nb.vec <- mu.nb.vec
  out$alpha.pl.vec <- alpha.pl.vec
  out$mean.pl.vec <- mean.pl.vec
  out$median.pl.vec <- median.pl.vec
  out$exp.E.log.y <- exp.E.log.y
  out$y.hat.pl.E.inv.y <- y.hat.pl.E.inv.y
  out$par.all <- par.all
  out$BIC <- BIC
  out$AIC <- AIC

  return(out)
}

#est.par <- Est.Obj$par.mat
#x.obj <- OBS.X.obj

#Takes parameters and x as input and predicts y (using mean or median of the pl component)
prediction.znbpl.fun <- function(x.obj, est.par) {
  x.multinom.zc <- x.obj$X.multinom.ZC
  x.multinom.pl <- x.obj$X.multinom.PL
  x.nb <- x.obj$X.NB
  x.pl <- x.obj$X.PL

  n <- max(nrow(x.multinom.zc), nrow(x.multinom.pl), nrow(x.nb), nrow(x.pl))

  if (is.null(x.multinom.zc)) {
    x.multinom.zc.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.multinom.zc.extended <- cbind(1, x.multinom.zc)
  }
  if (is.null(x.multinom.pl)) {
    x.multinom.pl.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.multinom.pl.extended <- cbind(1, x.multinom.pl)
  }
  if (is.null(x.nb)) {
    x.nb.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.nb.extended <- cbind(1, x.nb)
  }
  if (is.null(x.pl)) {
    x.pl.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.pl.extended <- cbind(1, x.pl)
  }

  n.beta.zc.multinomial <- dim(x.multinom.zc.extended)[2]
  n.beta.pl.multinomial <- dim(x.multinom.pl.extended)[2]
  n.theta.multinomial <- n.beta.zc.multinomial + n.beta.pl.multinomial
  n.beta.nb <- dim(x.nb.extended)[2]
  n.beta.pl <- dim(x.pl.extended)[2]

  beta.zc.multinomial.old <- est.par$Beta.multinom.ZC
  beta.pl.multinomial.old <- est.par$Beta.multinom.PL
  beta.nb.old <- est.par$Beta.NB
  alpha.nb.old <- est.par$Alpha.NB
  theta.nb.old <- c(beta.nb.old, alpha.nb.old)
  beta.pl.old <- est.par$Beta.PL
  c.pl <- est.par$C

  props.old <- matrix(NA, nrow = n, ncol = 3)
  for (i in 1:n) {
    denominator <- 1 +
      exp(t(beta.zc.multinomial.old) %*% x.multinom.zc.extended[i, ]) +
      exp(t(beta.pl.multinomial.old) %*% x.multinom.pl.extended[i, ])
    props.old[i, 1] <- exp(
      t(beta.zc.multinomial.old) %*% x.multinom.zc.extended[i, ]
    ) /
      denominator
    props.old[i, 2] <- 1 / denominator
    props.old[i, 3] <- exp(
      t(beta.pl.multinomial.old) %*% x.multinom.pl.extended[i, ]
    ) /
      denominator
  }

  #Mean conditional on negative binomial component
  #Pareto shape parameter
  mu.nb.vec <- c()
  alpha.pl.vec <- c()
  mean.pl.vec <- c()
  median.pl.vec <- c()
  exp.E.log.y <- c()
  E.inv.y <- c()
  y.hat.plmedian <- c()
  y.hat.plmean <- c()
  y.hat.plexpElogy <- c()
  y.hat.pl.E.inv.y <- c()
  for (i in 1:n) {
    mu.nb.vec[i] <- exp(x.nb.extended[i, ] %*% est.par$Beta.NB)
    alpha.pl.vec[i] <- exp(x.pl.extended[i, ] %*% est.par$Beta.PL)
    exp.E.log.y[i] <- c.pl * exp(1 / alpha.pl.vec[i])
    E.inv.y[i] <- alpha.pl.vec[i] / (c.pl * (alpha.pl.vec[i] + 1))
    if (alpha.pl.vec[i] > 1) {
      mean.pl.vec[i] <- alpha.pl.vec[i] * c.pl / (alpha.pl.vec[i] - 1)
    } else {
      mean.pl.vec[i] <- NA
    }
    median.pl.vec[i] <- c.pl * (2)^(1 / alpha.pl.vec[i])
    y.hat.plmedian[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] * median.pl.vec[i]
    y.hat.plmean[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] * mean.pl.vec[i]
    y.hat.plexpElogy[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] * exp.E.log.y[i]
    y.hat.pl.E.inv.y[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] / E.inv.y[i]
  }
  out <- list()
  out$y.hat.plmean <- y.hat.plmean
  out$y.hat.plmedian <- y.hat.plmedian
  out$y.hat.plexpElogy <- y.hat.plexpElogy
  out$y.hat.pl.E.inv.y <- y.hat.pl.E.inv.y
  return(out)
}


#znb <- Znb.Obj
#x <- OBS.X.obj$X.multinom.ZC
prediction.znb.fun <- function(x, znb) {
  coefficients <- znb$coefficients
  theta <- znb$theta
  n <- dim(x)[1]

  if (is.null(x)) {
    x.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.extended <- cbind(1, x)
  }

  n.beta <- dim(x.extended)[2]

  beta.zc.multinomial.old <- coefficients$zero
  beta.nb.old <- coefficients$count
  alpha.nb.old <- 1 / theta #Not sure. Could be 1/theta
  theta.nb.old <- c(beta.nb.old, alpha.nb.old)

  props.old <- matrix(NA, nrow = n, ncol = 2)
  for (i in 1:n) {
    denominator <- 1 + exp(t(beta.zc.multinomial.old) %*% x.extended[i, ])
    props.old[i, 1] <- exp(t(beta.zc.multinomial.old) %*% x.extended[i, ]) /
      denominator
    props.old[i, 2] <- 1 / denominator
  }

  #Mean conditional on negative binomial component
  #Pareto shape parameter
  mu.nb.vec <- c()
  y.hat.mean <- c()

  for (i in 1:n) {
    mu.nb.vec[i] <- exp(x.extended[i, ] %*% beta.nb.old)
    y.hat.mean[i] <- props.old[i, 2] * mu.nb.vec[i]
  }
  return(y.hat.mean)
}


#nb <- nb.obj
#x <- OBS.X.ins
prediction.nb.fun <- function(x, nb) {
  coefficients <- nb$coefficients
  theta <- nb$theta
  n <- dim(x)[1]

  if (is.null(x)) {
    x.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.extended <- cbind(1, x)
  }

  n.beta <- dim(x.extended)[2]

  beta.nb.old <- coefficients
  alpha.nb.old <- 1 / theta #Not sure. Could be 1/theta
  theta.nb.old <- c(beta.nb.old, alpha.nb.old)

  #Mean conditional on negative binomial component
  #Pareto shape parameter
  mu.nb.vec <- c()
  y.hat.mean <- c()

  for (i in 1:n) {
    mu.nb.vec[i] <- exp(x.extended[i, ] %*% beta.nb.old)
    y.hat.mean[i] <- mu.nb.vec[i]
  }
  return(y.hat.mean)
}


marginal.effect.nb.fun <- function(sign.level, j, x, x.lim, dx, beta.nb) {
  n <- dim(x)[1]

  beta.nb <- as.matrix(beta.nb)

  if (is.null(x)) {
    x.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.extended <- cbind(1, x)
  }

  j <- j + 1 #To compensate for the intercept

  beta.nb.j <- beta.nb[j]

  #i <- 1  #individual by individual
  x.j.all <- seq(x.lim[1], x.lim[2], by = dx)
  nx.j <- length(x.j.all)

  margEff.nb <- matrix(NA, nrow = n, ncol = nx.j)
  for (l in 1:nx.j) {
    for (i in 1:n) {
      x.il.extended <- as.matrix(x.extended[i, ])
      x.il.extended[j] <- x.j.all[l]
      exp.beta.nb.x <- exp(t(beta.nb) %*% x.il.extended)
      margEff.nb[i, l] <- exp.beta.nb.x
    }
  }
  mean.margEff.nb <- colMeans(margEff.nb)
  sort.margEff.nb <- matrix(NA, nrow = n, ncol = nx.j)
  upper.conf.nb <- c()
  lower.conf.nb <- c()
  lower.index <- round(0.5 * sign.level * n)
  upper.index <- round((1 - 0.5 * sign.level) * n)
  for (l in 1:nx.j) {
    sort.margEff.nb[, l] <- sort(margEff.nb[, l])
    lower.conf.nb[l] <- sort.margEff.nb[lower.index, l]
    upper.conf.nb[l] <- sort.margEff.nb[upper.index, l]
  }
  out <- list()
  out$x <- x.j.all
  out$mean.margEff.nb <- mean.margEff.nb
  out$lower.conf.nb <- lower.conf.nb
  out$upper.conf.nb <- upper.conf.nb
  return(out)
}

#j is an integer representing the chosen covariate
#dx is the interval in the chosen covariate
#x.lim is the min and max of the chosen covariate
marginal.effect.znb.fun <- function(
  sign.level,
  j,
  x,
  x.lim,
  dx,
  gamma.zc,
  beta.nb
) {
  n <- dim(x)[1]

  beta.nb <- as.matrix(beta.nb)
  gamma.zc <- as.matrix(gamma.zc)

  if (is.null(x)) {
    x.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.extended <- cbind(1, x)
  }

  j <- j + 1 #To compensate for the intercept

  beta.nb.j <- beta.nb[j]
  gamma.zc.j <- gamma.zc[j]

  #i <- 1  #individual by individual
  x.j.all <- seq(x.lim[1], x.lim[2], by = dx)
  nx.j <- length(x.j.all)

  margEff.znb <- matrix(NA, nrow = n, ncol = nx.j)
  prop.z <- matrix(NA, nrow = n, ncol = nx.j)
  prop.nb <- matrix(NA, nrow = n, ncol = nx.j)
  for (l in 1:nx.j) {
    for (i in 1:n) {
      x.il.extended <- as.matrix(x.extended[i, ])
      x.il.extended[j] <- x.j.all[l]
      exp.gamma.zc.x <- exp(t(gamma.zc) %*% x.il.extended)
      exp.beta.nb.x <- exp(t(beta.nb) %*% x.il.extended)
      margEff.znb[i, l] <- exp.beta.nb.x / (1 + exp.gamma.zc.x)
      #margEff.znb[i,l] <- beta.nb.j - gamma.zc.j*exp.gamma.zc.x/(1+exp.gamma.zc.x)
      prop.nb[i, l] <- 1 / (exp.gamma.zc.x + 1)
      prop.z[i, l] <- 1 - prop.nb[i, l]
    }
  }
  mean.margEff.znb <- colMeans(margEff.znb)
  sort.margEff.znb <- matrix(NA, nrow = n, ncol = nx.j)
  mean.prop.z <- colMeans(prop.z)
  mean.prop.nb <- colMeans(prop.nb)
  upper.conf.znb <- c()
  lower.conf.znb <- c()
  upper.conf.prop.z <- c()
  lower.conf.prop.z <- c()
  upper.conf.prop.nb <- c()
  lower.conf.prop.nb <- c()
  sort.prop.z <- matrix(NA, nrow = n, ncol = nx.j)
  sort.prop.nb <- matrix(NA, nrow = n, ncol = nx.j)
  lower.index <- round(0.5 * sign.level * n)
  upper.index <- round((1 - 0.5 * sign.level) * n)
  for (l in 1:nx.j) {
    sort.margEff.znb[, l] <- sort(margEff.znb[, l])
    sort.prop.z[, l] <- sort(prop.z[, l])
    sort.prop.nb[, l] <- sort(prop.nb[, l])
    lower.conf.znb[l] <- sort.margEff.znb[lower.index, l]
    lower.conf.prop.z[l] <- sort.prop.z[lower.index, l]
    lower.conf.prop.nb[l] <- sort.prop.nb[lower.index, l]
    upper.conf.znb[l] <- sort.margEff.znb[upper.index, l]
    upper.conf.prop.z[l] <- sort.prop.z[upper.index, l]
    upper.conf.prop.nb[l] <- sort.prop.nb[upper.index, l]
  }
  out <- list()
  out$x <- x.j.all
  out$mean.margEff.znb <- mean.margEff.znb
  out$lower.conf.znb <- lower.conf.znb
  out$upper.conf.znb <- upper.conf.znb
  out$mean.prop.z <- mean.prop.z
  out$mean.prop.nb <- mean.prop.nb
  out$lower.conf.prop.z <- lower.conf.prop.z
  out$upper.conf.prop.z <- upper.conf.prop.z
  out$lower.conf.prop.nb <- lower.conf.prop.nb
  out$upper.conf.prop.nb <- upper.conf.prop.nb
  return(out)
}


# sign.level <- 0.05
# j <- 9 #is ilustrative
# #j <- 11
# dx <- 0.2
# x <- OBS.X
# x.lim <- c(min(x[,j]),max(x[,j]))
# #x.lim <- c(0,75)
# #x.lim <- c(min(x[,j]),50) #for j=7
# gamma.zc <- Est.Obj$par.mat$Beta.multinom.ZC
# gamma.pl <- Est.Obj$par.mat$Beta.multinom.PL[1:5]
# gamma.pl <- c(gamma.pl,0,Est.Obj$par.mat$Beta.multinom.PL[6:12])
# beta.nb <- Est.Obj$par.mat$Beta.NB
# beta.pl <- Est.Obj$par.mat$Beta.PL
# beta.pl <- c(beta.pl,rep(0,12))
# c.pl <- Est.Obj$par.mat$C

marginal.effect.znbpl.fun <- function(
  sign.level,
  j,
  x,
  x.lim,
  dx,
  gamma.zc,
  gamma.pl,
  beta.nb,
  beta.pl,
  c.pl
) {
  n <- dim(x)[1]

  gamma.zc <- as.matrix(gamma.zc)
  gamma.pl <- as.matrix(gamma.pl)
  beta.nb <- as.matrix(beta.nb)
  beta.pl <- as.matrix(beta.pl)

  if (is.null(x)) {
    x.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.extended <- cbind(1, x)
  }

  j <- j + 1 #To compensate for the intercept

  gamma.zc.j <- gamma.zc[j]
  gamma.pl.j <- gamma.pl[j]
  beta.nb.j <- beta.nb[j]
  beta.pl.j <- beta.pl[j]

  x.j.all <- seq(x.lim[1], x.lim[2], by = dx)
  nx.j <- length(x.j.all)

  margEff.znbpl <- matrix(NA, nrow = n, ncol = nx.j)
  prop.z <- matrix(NA, nrow = n, ncol = nx.j)
  prop.nb <- matrix(NA, nrow = n, ncol = nx.j)
  prop.p <- matrix(NA, nrow = n, ncol = nx.j)
  for (l in 1:nx.j) {
    for (i in 1:n) {
      x.il.extended <- as.matrix(x.extended[i, ])
      x.il.extended[j] <- x.j.all[l]
      gamma.plTx <- t(gamma.pl) %*% x.il.extended
      exp.gamma.zc.x <- exp(t(gamma.zc) %*% x.il.extended)
      exp.gamma.pl.x <- exp(t(gamma.pl) %*% x.il.extended)
      exp.beta.nb.x <- exp(t(beta.nb) %*% x.il.extended)
      exp.beta.pl.x <- exp(t(beta.pl) %*% x.il.extended)

      numerator <- exp.beta.nb.x +
        c.pl * exp.gamma.pl.x / exp.beta.pl.x +
        c.pl * exp.gamma.pl.x
      denominator <- 1 + exp.gamma.zc.x + exp.gamma.pl.x

      margEff.znbpl[i, l] <- numerator / denominator
      prop.nb[i, l] <- 1 / (exp.beta.pl.x + exp.gamma.zc.x + 1)
      prop.p[i, l] <- exp.beta.pl.x / (exp.beta.pl.x + exp.gamma.zc.x + 1)
      prop.z[i, l] <- exp.gamma.zc.x / (exp.beta.pl.x + exp.gamma.zc.x + 1)
    }
  }
  mean.margEff.znbpl <- colMeans(margEff.znbpl)
  upper.conf.znbpl <- c()
  lower.conf.znbpl <- c()
  sort.margEff.znbpl <- matrix(NA, nrow = n, ncol = nx.j)
  mean.prop.z <- colMeans(prop.z)
  mean.prop.nb <- colMeans(prop.nb)
  mean.prop.p <- colMeans(prop.p)
  upper.conf.prop.z <- c()
  lower.conf.prop.z <- c()
  upper.conf.prop.nb <- c()
  lower.conf.prop.nb <- c()
  upper.conf.prop.p <- c()
  lower.conf.prop.p <- c()
  sort.prop.z <- matrix(NA, nrow = n, ncol = nx.j)
  sort.prop.nb <- matrix(NA, nrow = n, ncol = nx.j)
  sort.prop.p <- matrix(NA, nrow = n, ncol = nx.j)
  lower.index <- round(0.5 * sign.level * n)
  upper.index <- round((1 - 0.5 * sign.level) * n)
  for (l in 1:nx.j) {
    sort.margEff.znbpl[, l] <- sort(margEff.znbpl[, l])
    lower.conf.znbpl[l] <- sort.margEff.znbpl[lower.index, l]
    upper.conf.znbpl[l] <- sort.margEff.znbpl[upper.index, l]
    sort.prop.z[, l] <- sort(prop.z[, l])
    sort.prop.nb[, l] <- sort(prop.nb[, l])
    sort.prop.p[, l] <- sort(prop.p[, l])
    lower.conf.prop.z[l] <- sort.prop.z[lower.index, l]
    lower.conf.prop.nb[l] <- sort.prop.nb[lower.index, l]
    lower.conf.prop.p[l] <- sort.prop.p[lower.index, l]
    upper.conf.prop.z[l] <- sort.prop.z[upper.index, l]
    upper.conf.prop.nb[l] <- sort.prop.nb[upper.index, l]
    upper.conf.prop.p[l] <- sort.prop.p[upper.index, l]
  }
  out <- list()
  out$x <- x.j.all
  out$mean.margEff.znbpl <- mean.margEff.znbpl
  out$lower.conf.znbpl <- lower.conf.znbpl
  out$upper.conf.znbpl <- upper.conf.znbpl
  out$mean.prop.z <- mean.prop.z
  out$mean.prop.nb <- mean.prop.nb
  out$mean.prop.p <- mean.prop.p
  out$lower.conf.prop.z <- lower.conf.prop.z
  out$upper.conf.prop.z <- upper.conf.prop.z
  out$lower.conf.prop.nb <- lower.conf.prop.nb
  out$upper.conf.prop.nb <- upper.conf.prop.nb
  out$lower.conf.prop.p <- lower.conf.prop.p
  out$upper.conf.prop.p <- upper.conf.prop.p
  return(out)
  #plot(x.j.all,mean.margEff.znbpl)
}


plinfl.nb.regression.fun <- function(y, x.obj, ini.val, control) {
  #Initialize data and parameters

  x.multinom.zc <- x.obj$X.multinom.ZC
  x.multinom.pl <- x.obj$X.multinom.PL
  x.nb <- x.obj$X.NB
  x.pl <- x.obj$X.PL

  n <- max(nrow(x.nb), nrow(x.pl), nrow(x.multinom.zc), nrow(x.multinom.pl))

  if (is.null(dim(x.multinom.zc))) {
    x.multinom.zc.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.multinom.zc.extended <- cbind(1, x.multinom.zc)
  }
  if (is.null(dim(x.multinom.pl))) {
    x.multinom.pl.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.multinom.pl.extended <- cbind(1, x.multinom.pl)
  }
  if (is.null(dim(x.nb))) {
    x.nb.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.nb.extended <- cbind(1, x.nb)
  }
  if (is.null(dim(x.pl))) {
    x.pl.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.pl.extended <- cbind(1, x.pl)
  }

  n.beta.zc.multinomial <- dim(x.multinom.zc.extended)[2]
  n.beta.pl.multinomial <- dim(x.multinom.pl.extended)[2]
  n.theta.multinomial <- n.beta.zc.multinomial + n.beta.pl.multinomial
  n.beta.nb <- dim(x.nb.extended)[2]
  n.beta.pl <- dim(x.pl.extended)[2]

  prel.val <- ini.val

  n <- max(nrow(x.nb), nrow(x.pl), nrow(x.multinom.zc), nrow(x.multinom.pl))

  #Determine which values of c to be investigated
  c.range <- unique(sort(y))[
    unique(sort(y)) >= control$c.lim[1] & unique(sort(y)) <= control$c.lim[2]
  ]

  ###When c is not well known we only run short versions of zerinfl.nb.pl.reg.cond.c.fun() in order to save time
  control.warmup <- control
  control.warmup$max.no.em.steps <- control$max.no.em.steps.warmup

  ###Initializing the warm-up phase when the update stops after control$max.no.em.steps.warmup EM steps
  #If the update of c is less than 1, stop the warmup
  c.abs.diff <- 100
  log.lik.vec.all <- NULL
  cat('Begin warm-up', '\n', sep = '')
  while (c.abs.diff > 0) {
    log.lik.vec <- c()
    est.obj <- zerinfl.nb.pl.reg.cond.c.fun(y, x.obj, prel.val, control.warmup)
    log.lik.vec.all <- c(log.lik.vec.all, est.obj$log.lik.vec)
    prel.val <- est.obj$par.mat
    props.old <- prel.val$Props
    beta.nb.old <- prel.val$Beta.NB
    alpha.nb.old <- prel.val$Alpha.NB
    beta.pl.old <- prel.val$Beta.PL
    beta.zc.multinomial.old <- ini.val$Beta.multinom.ZC
    beta.pl.multinomial.old <- prel.val$Beta.multinom.PL
    c.pl <- prel.val$C

    for (k in 1:length(c.range)) {
      log.lik.vec[k] <- log_lik_fun(
        beta.zc.multinomial.old,
        beta.pl.multinomial.old,
        beta.nb.old,
        alpha.nb.old,
        beta.pl.old,
        c.range[k],
        x.multinom.zc.extended,
        x.multinom.pl.extended,
        x.nb.extended,
        x.pl.extended,
        y
      )
    }

    c.pl.new <- c.range[which(log.lik.vec == max(log.lik.vec))]

    if (is.infinite(max(log.lik.vec))) {
      stop(
        'The log-likelihood is infinite for all tested values of c in the c-range. Try expanding the c-range with the c.lim argument of the function'
      )
    }

    c.abs.diff <- abs(c.pl.new - prel.val$C)

    prel.val$C <- c.pl.new
    func.val <- max(log.lik.vec)

    log.lik.vec.all <- c(log.lik.vec.all, func.val)

    cat(
      'The new c is ',
      c.pl.new,
      '. The function value is ',
      round(func.val, 4),
      '\n',
      sep = ''
    )
    prel.val$C <- c.pl.new
  }

  ###The end phase, when c is presumably rather accurate
  c.abs.diff <- 100
  cat('End warm-up. Run until convergence', '\n', sep = '')
  while (c.abs.diff > 0) {
    log.lik.vec <- c()

    est.obj <- plinfl.nb.reg.cond.c.fun(y, x.obj, prel.val, control)
    est.obj$par.mat$Beta.multinom.ZC <- ini.val$Beta.multinom.ZC
    log.lik.vec.all <- c(log.lik.vec.all, est.obj$log.lik.vec)
    prel.val <- est.obj$par.mat
    props.old <- prel.val$Props
    beta.nb.old <- prel.val$Beta.NB
    alpha.nb.old <- prel.val$Alpha.NB
    beta.pl.old <- prel.val$Beta.PL
    beta.zc.multinomial.old <- ini.val$Beta.multinom.ZC
    beta.pl.multinomial.old <- est.obj$par.mat$Beta.multinom.PL
    c.pl <- prel.val$C

    for (k in 1:length(c.range)) {
      log.lik.vec[k] <- log_lik_fun(
        beta.zc.multinomial.old,
        beta.pl.multinomial.old,
        beta.nb.old,
        alpha.nb.old,
        beta.pl.old,
        c.range[k],
        x.multinom.zc.extended,
        x.multinom.pl.extended,
        x.nb.extended,
        x.pl.extended,
        y
      )
    }

    c.pl.new <- c.range[which(log.lik.vec == max(log.lik.vec))]

    if (is.infinite(max(log.lik.vec))) {
      stop(
        'The log-likelihood is infinite for all tested values of c in the c-range. Try expanding the c-range with the c.lim argument of the function'
      )
    }

    c.abs.diff <- abs(c.pl.new - prel.val$C)

    prel.val$C <- c.pl.new
    func.val <- max(log.lik.vec.all)
    log.lik.vec.all <- c(log.lik.vec.all, func.val)

    cat(
      'The new c is ',
      c.pl.new,
      '. The function value is ',
      round(func.val, 4),
      '\n',
      sep = ''
    )
    prel.val$C <- c.pl.new
  }

  final.val <- prel.val

  #Mean conditional on negative binomial component
  #Pareto shape parameter
  mu.nb.vec <- c()
  alpha.pl.vec <- c()
  mean.pl.vec <- c()
  exp.E.log.y <- c() #exp(E(log(y))), a third alternative for prediction
  median.pl.vec <- c()
  E.inv.y <- c()
  y.hat.plmedian <- c()
  y.hat.plmean <- c()
  y.hat.plexpElogy <- c()
  y.hat.pl.E.inv.y <- c()
  for (i in 1:n) {
    mu.nb.vec[i] <- exp(x.nb.extended[i, ] %*% prel.val$Beta.NB)
    alpha.pl.vec[i] <- exp(x.pl.extended[i, ] %*% prel.val$Beta.PL)
    exp.E.log.y[i] <- c.pl.new * exp(1 / alpha.pl.vec[i])
    E.inv.y[i] <- alpha.pl.vec[i] / (c.pl.new * (alpha.pl.vec[i] + 1))
    if (alpha.pl.vec[i] > 1) {
      mean.pl.vec[i] <- alpha.pl.vec[i] * c.pl.new / (alpha.pl.vec[i] - 1)
    } else {
      mean.pl.vec[i] <- NA
    }
    median.pl.vec[i] <- median.pl.vec[i] <- c.pl.new * (2)^(1 / alpha.pl.vec[i])
    y.hat.plmedian[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] * median.pl.vec[i]
    y.hat.plmean[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] * mean.pl.vec[i]
    y.hat.plexpElogy[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] * exp.E.log.y[i]
    y.hat.pl.E.inv.y[i] <- props.old[i, 2] *
      mu.nb.vec[i] +
      props.old[i, 3] / E.inv.y[i]
  }

  par.all <- c(
    final.val$Beta.multinom.ZC,
    final.val$Beta.multinom.PL,
    final.val$Beta.NB,
    final.val$Alpha.NB,
    final.val$Beta.PL,
    final.val$C
  )
  n.par <- length(par.all)
  BIC <- log(n) * n.par - 2 * func.val
  AIC <- 2 * n.par - 2 * func.val

  #Gather results from estimation
  out <- list()

  out$control <- control
  out$par.mat <- final.val
  out$log.lik.vec.all <- log.lik.vec.all
  out$log.lik <- log.lik.vec.all[length(log.lik.vec.all)]
  out$resp <- est.obj$resp
  out$converge <- est.obj$converge
  out$ini.val <- ini.val
  out$x.nb <- x.nb
  out$x.pl <- x.pl
  out$x.multinom.zc <- x.multinom.zc
  out$x.multinom.pl <- x.multinom.pl
  out$median.pl.vec <- median.pl.vec
  out$mean.pl.vec <- mean.pl.vec
  out$y <- y
  out$y.hat.plmedian <- y.hat.plmedian
  out$y.hat.plmean <- y.hat.plmean
  out$y.hat.plexpElogy <- y.hat.plexpElogy
  out$mu.nb.vec <- mu.nb.vec
  out$alpha.pl.vec <- alpha.pl.vec
  out$mean.pl.vec <- mean.pl.vec
  out$median.pl.vec <- median.pl.vec
  out$exp.E.log.y <- exp.E.log.y
  out$y.hat.pl.E.inv.y <- y.hat.pl.E.inv.y
  out$par.all <- par.all
  out$BIC <- BIC
  out$AIC <- AIC

  return(out)
}


## Internal functions for the EVZINB and EVINB functions
plinfl.nb.reg.cond.c.fun <- function(y, x.obj, ini.val, control) {
  #Initialize parameters and data
  x.multinom.zc <- x.obj$X.multinom.ZC
  x.multinom.pl <- x.obj$X.multinom.PL

  x.nb <- x.obj$X.NB
  x.pl <- x.obj$X.PL

  n <- max(nrow(x.nb), nrow(x.pl), nrow(x.multinom.zc), nrow(x.multinom.pl))

  if (is.null(dim(x.multinom.zc))) {
    x.multinom.zc.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.multinom.zc.extended <- cbind(1, x.multinom.zc)
  }
  if (is.null(dim(x.multinom.pl))) {
    x.multinom.pl.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.multinom.pl.extended <- cbind(1, x.multinom.pl)
  }
  if (is.null(dim(x.nb))) {
    x.nb.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.nb.extended <- cbind(1, x.nb)
  }
  if (is.null(dim(x.pl))) {
    x.pl.extended <- matrix(1, nrow = n, ncol = 1)
  } else {
    x.pl.extended <- cbind(1, x.pl)
  }

  n.beta.zc.multinomial <- dim(x.multinom.zc.extended)[2]
  n.beta.pl.multinomial <- dim(x.multinom.pl.extended)[2]
  n.theta.multinomial <- n.beta.zc.multinomial + n.beta.pl.multinomial
  n.beta.nb <- dim(x.nb.extended)[2]
  n.beta.pl <- dim(x.pl.extended)[2]

  beta.zc.multinomial.old <- ini.val$Beta.multinom.ZC
  beta.pl.multinomial.old <- ini.val$Beta.multinom.PL
  beta.nb.old <- ini.val$Beta.NB
  alpha.nb.old <- ini.val$Alpha.NB
  theta.nb.old <- c(beta.nb.old, alpha.nb.old)
  beta.pl.old <- ini.val$Beta.PL
  c.pl <- ini.val$C

  #Initial log likelihood value

  #func.val.initial <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)
  func.val.initial <- log_lik_fun(
    beta.zc.multinomial.old,
    beta.pl.multinomial.old,
    beta.nb.old,
    alpha.nb.old,
    beta.pl.old,
    c.pl,
    x.multinom.zc.extended,
    x.multinom.pl.extended,
    x.nb.extended,
    x.pl.extended,
    y
  )

  # The maximum change after one EM step needs to be small in order for the algorithm to stop
  max.abs.par.diff <- 100

  #Number of nas produced
  na.beta.nb <- 0
  na.alpha.nb <- 0
  na.beta.pl <- 0
  na.beta.mult.zc <- 0
  na.beta.mult.pl <- 0

  ################
  #Start EM algorithm
  ###################
  i.em <- 1
  func.val.vec <- c()
  func.val.old <- -1e50
  max.no.em.steps <- control$max.no.em.steps
  max.diff.par <- control$max.diff.par
  max.upd.par <- control$max.upd.par.nb
  no.m.bfgs.steps <- control$no.m.bfgs.steps.nb

  #max.no.em.steps <- 21

  while (i.em < max.no.em.steps & max.abs.par.diff > max.diff.par) {
    beta.zc.multinomial.start.em <- ini.val$Beta.multinom.ZC
    beta.pl.multinomial.start.em <- beta.pl.multinomial.old
    beta.nb.start.em <- beta.nb.old
    alpha.nb.start.em <- abs(alpha.nb.old)
    beta.pl.start.em <- beta.pl.old

    par.start.em <- c(
      beta.zc.multinomial.start.em,
      beta.pl.multinomial.start.em,
      beta.nb.start.em,
      alpha.nb.start.em,
      beta.pl.start.em
    )

    #Update parameters with BFGS
    upd.obj <- update_bfgs_fun(
      ini.val$Beta.multinom.ZC,
      beta.pl.multinomial.old,
      beta.nb.old,
      alpha.nb.old,
      beta.pl.old,
      c.pl,
      x.multinom.zc.extended,
      x.multinom.pl.extended,
      x.nb.extended,
      x.pl.extended,
      y,
      max.upd.par,
      control$no.m.bfgs.steps.nb
    )

    if (sum(is.na(upd.obj$beta_nb_old)) == 0) {
      beta.nb.new <- upd.obj$beta_nb_old
    } else {
      beta.nb.new <- beta.nb.old
      na.beta.nb <- na.beta.nb + length(is.na(upd.obj$beta_nb_old))
    }

    if (sum(is.na(upd.obj$alpha_nb_old)) == 0) {
      alpha.nb.new <- upd.obj$alpha_nb_old
    } else {
      alpha.nb.new <- alpha.nb.old
      na.alpha.nb <- na.alpha.nb + length(is.na(upd.obj$alpha_nb_old))
    }

    if (sum(is.na(upd.obj$beta_pl_old)) == 0) {
      beta.pl.new <- upd.obj$beta_pl_old
    } else {
      beta.pl.new <- beta.pl.old
      na.beta.pl <- na.beta.pl + length(is.na(upd.obj$beta_pl_old))
    }

    if (sum(is.na(upd.obj$gamma_z_old)) == 0) {
      beta.zc.multinomial.new <- ini.val$Beta.multinom.ZC
    } else {
      beta.zc.multinomial.new <- ini.val$Beta.multinom.ZC
      na.beta.mult.zc <- na.beta.mult.zc + length(is.na(upd.obj$gamma_z_old))
    }

    if (sum(is.na(upd.obj$gamma_pl_old)) == 0) {
      beta.pl.multinomial.new <- upd.obj$gamma_pl_old
    } else {
      beta.pl.multinomial.new <- beta.pl.multinomial.old
      na.beta.mult.pl <- na.beta.mult.pl + length(is.na(upd.obj$gamma_pl_old))
    }

    # alpha.nb.new <- upd.obj$alpha_nb_old
    # beta.pl.new <- upd.obj$beta_pl_old
    # beta.zc.multinomial.new <- upd.obj$gamma_z_old
    # beta.pl.multinomial.new <- upd.obj$gamma_pl_old

    par.end.em <- c(
      beta.zc.multinomial.new,
      beta.pl.multinomial.new,
      beta.nb.new,
      alpha.nb.new,
      beta.pl.new
    )
    par.diff.em <- par.end.em - par.start.em

    ##############Make sure the likelihood increases

    ################################################
    ##################Optimize theta.nb.old
    nb.log.lik.optim <- function(eta) {
      #Go back eta times the step that was already made
      theta.nb.old <- c(beta.nb.old, alpha.nb.old) +
        eta * upd.obj$change_nb_bfgs # change.nb.bfgs
      beta.nb.tmp <- theta.nb.old[1:n.beta.nb]
      alpha.nb.tmp <- theta.nb.old[n.beta.nb + 1]
      return(
        -1.0 *
          log_lik_fun(
            ini.val$Beta.multinom.ZC,
            beta.pl.multinomial.old,
            beta.nb.tmp,
            alpha.nb.tmp,
            beta.pl.old,
            c.pl,
            x.multinom.zc.extended,
            x.multinom.pl.extended,
            x.nb.extended,
            x.pl.extended,
            y
          )
      )
    }

    if (is.na(upd.obj$func_val_after_nb) == FALSE) {
      if (upd.obj$func_val_after_nb < upd.obj$func_val_before_bfgs) {
        if (sum(is.na(upd.obj$change_nb_bfgs)) == 0) {
          ###########################
          change.nb.obj <- optimise(
            f = nb.log.lik.optim,
            interval = control$eta.int
          )
          eta.nb <- change.nb.obj$minimum

          theta.nb.old <- c(beta.nb.old, alpha.nb.old)
          theta.nb.after.optim <- theta.nb.old + eta.nb * upd.obj$change_nb_bfgs
          beta.nb.after.optim <- theta.nb.after.optim[1:n.beta.nb]
          alpha.nb.after.optim <- theta.nb.after.optim[n.beta.nb + 1]

          func.val.after.nb.optim <- log_lik_fun(
            ini.val$Beta.multinom.ZC,
            beta.pl.multinomial.old,
            beta.nb.after.optim,
            alpha.nb.after.optim,
            beta.pl.old,
            c.pl,
            x.multinom.zc.extended,
            x.multinom.pl.extended,
            x.nb.extended,
            x.pl.extended,
            y
          )

          theta.nb.old <- theta.nb.after.optim
          beta.nb.new <- theta.nb.old[1:n.beta.nb]
          alpha.nb.new <- theta.nb.old[n.beta.nb + 1]
        }
      }
    }
    ##################################
    ##############################################

    ################################################
    ##################Optimize theta.pl.old
    pl.log.lik.optim <- function(eta) {
      #Go back eta times the step that was already made
      beta.pl.tmp <- beta.pl.old + eta * upd.obj$change_pl_bfgs
      return(
        -1.0 *
          log_lik_fun(
            ini.val$Beta.multinom.ZC,
            beta.pl.multinomial.old,
            beta.nb.old,
            alpha.nb.old,
            beta.pl.tmp,
            c.pl,
            x.multinom.zc.extended,
            x.multinom.pl.extended,
            x.nb.extended,
            x.pl.extended,
            y
          )
      )
    }

    if (is.na(upd.obj$func_val_after_pl) == FALSE) {
      if (upd.obj$func_val_after_pl < upd.obj$func_val_before_bfgs) {
        ###########################
        #if(det(d2Qdbeta2.pl)>0){
        change.pl.obj <- optimise(
          f = pl.log.lik.optim,
          interval = control$eta.int
        )
        eta.pl <- change.pl.obj$minimum

        beta.pl.after.optim <- beta.pl.old + eta.pl * upd.obj$change_pl_bfgs

        func.val.after.pl.optim <- log_lik_fun(
          ini.val$Beta.multinom.ZC,
          beta.pl.multinomial.old,
          beta.nb.old,
          alpha.nb.old,
          beta.pl.after.optim,
          c.pl,
          x.multinom.zc.extended,
          x.multinom.pl.extended,
          x.nb.extended,
          x.pl.extended,
          y
        )

        beta.pl.new <- beta.pl.after.optim
      }
    }

    ################################################
    ##################Optimize theta.pl.old
    # zc.multinomial.log.lik.optim <- function(eta){
    #   #Go back eta times the step that was already made
    #   beta.zc.multinomial.tmp <- beta.zc.multinomial.old + eta*upd.obj$change_mult_z_bfgs
    #   return(-1.0*log_lik_fun(beta.zc.multinomial.tmp,beta.pl.multinomial.old,beta.nb.old,alpha.nb.old,beta.pl.old,c.pl,x.multinom.zc.extended,x.multinom.pl.extended,x.nb.extended,x.pl.extended,y))
    # }

    # if(is.na(upd.obj$func_val_after_mult_z)==FALSE){
    #   if(upd.obj$func_val_after_mult_z<upd.obj$func_val_before_bfgs){
    #     ###########################
    #     change.zc.multinomial.obj <- optimise(f=zc.multinomial.log.lik.optim,interval=control$eta.int)
    #     eta.zc.multinomial <- change.zc.multinomial.obj$minimum
    #
    #     beta.zc.multinomial.after.optim <- beta.zc.multinomial.old + eta.zc.multinomial*upd.obj$change_mult_z_bfgs
    #
    #     func.val.after.zc.multinomial.optim <- log_lik_fun(beta.zc.multinomial.after.optim,beta.pl.multinomial.old,beta.nb.old,alpha.nb.old,beta.pl.old,c.pl,x.multinom.zc.extended,x.multinom.pl.extended,x.nb.extended,x.pl.extended,y)
    #
    #     beta.zc.multinomial.new <- beta.zc.multinomial.after.optim
    #   }
    #
    # }

    ##################################
    ##############################################

    ################################################
    ##################Optimize theta.pl.old
    pl.multinomial.log.lik.optim <- function(eta) {
      #Go back eta times the step that was already made
      beta.pl.multinomial.tmp <- beta.pl.multinomial.old +
        eta * upd.obj$change_mult_pl_bfgs
      return(
        -1.0 *
          log_lik_fun(
            beta.zc.multinomial.old,
            beta.pl.multinomial.tmp,
            beta.nb.old,
            alpha.nb.old,
            beta.pl.old,
            c.pl,
            x.multinom.zc.extended,
            x.multinom.pl.extended,
            x.nb.extended,
            x.pl.extended,
            y
          )
      )
    }

    ###########################
    if (is.na(upd.obj$func_val_after_mult_pl) == FALSE) {
      if (upd.obj$func_val_after_mult_pl < upd.obj$func_val_before_bfgs) {
        change.pl.multinomial.obj <- optimise(
          f = pl.multinomial.log.lik.optim,
          interval = control$eta.int
        )
        eta.pl.multinomial <- change.pl.multinomial.obj$minimum

        beta.pl.multinomial.after.optim <- beta.pl.multinomial.old +
          eta.pl.multinomial * upd.obj$change_mult_pl_bfgs

        func.val.after.pl.multinomial.optim <- log_lik_fun(
          beta.zc.multinomial.old,
          beta.pl.multinomial.after.optim,
          beta.nb.old,
          alpha.nb.old,
          beta.pl.old,
          c.pl,
          x.multinom.zc.extended,
          x.multinom.pl.extended,
          x.nb.extended,
          x.pl.extended,
          y
        )

        beta.pl.multinomial.new <- beta.pl.multinomial.after.optim
      }
    }
    ##################################
    ##############################################

    ########################
    ############################## Finished updating parameters

    par.end.em <- c(
      ini.val$Beta.multinom.ZC,
      beta.pl.multinomial.new,
      beta.nb.new,
      alpha.nb.new,
      beta.pl.new
    )

    par.diff <- par.end.em - par.start.em
    max.abs.par.diff <- max(abs(par.diff))

    func.val.after.bfgs <- log_lik_fun(
      ini.val$Beta.multinom.ZC,
      beta.pl.multinomial.new,
      beta.nb.new,
      alpha.nb.new,
      beta.pl.new,
      c.pl,
      x.multinom.zc.extended,
      x.multinom.pl.extended,
      x.nb.extended,
      x.pl.extended,
      y
    )

    par.all.tmp <- c(
      ini.val$Beta.multinom.ZC,
      beta.pl.multinomial.old,
      beta.nb.old,
      alpha.nb.old,
      beta.pl.old
    )

    beta.nb.old <- beta.nb.new
    alpha.nb.old <- alpha.nb.new
    beta.pl.old <- beta.pl.new
    beta.zc.multinomial.old <- ini.val$Beta.multinom.ZC
    beta.pl.multinomial.old <- beta.pl.multinomial.new

    # if(func.val.after.bfgs<upd.obj$func_val_before_bfgs){
    #  print("The function decreased")
    # par.diff.mod <- 0.01*par.diff
    # par.all.tmp <- par.start.em + par.diff.mod
    # beta.nb.multinomial.old <- par.all.tmp[1:n.beta.nb.multinomial]
    # par.all.tmp <- par.all.tmp[-c(1:n.beta.nb.multinomial)]
    # beta.pl.multinomial.old <- par.all.tmp[1:n.beta.pl.multinomial]
    # par.all.tmp <- par.all.tmp[-c(1:n.beta.pl.multinomial)]
    # beta.nb.old <- par.all.tmp[1:n.beta.nb]
    # par.all.tmp <- par.all.tmp[-c(1:n.beta.nb)]
    # alpha.nb.old <- par.all.tmp[1]
    # par.all.tmp <- par.all.tmp[-1]
    # beta.pl.old <- par.all.tmp
    # rm(par.all.tmp)
    # }

    #func.val.old <- log.lik.fun(beta.zc.multinomial.old,beta.pl.multinomial.old,theta.nb.old,beta.pl.old,c.pl)

    #func.val.old <- func.val

    cat(
      'Iteration ',
      i.em,
      ': The max abs diff in parameters is ',
      max.abs.par.diff,
      '. The function value is ',
      round(func.val.after.bfgs, 4),
      '\n',
      sep = ''
    )

    func.val.vec[i.em] <- func.val.after.bfgs
    i.em <- i.em + 1
    # print(beta.nb.multinomial.old)
    # print(beta.pl.multinomial.old)
    # print(beta.nb.old)
    # print(beta.pl.old)
    # print(alpha.nb.old)
  }

  #Gather results after finishing the estimation
  res <- list()
  par.mat <- list()
  par.mat$Props <- upd.obj$prop
  par.mat$Beta.multinom.ZC <- ini.val$Beta.multinom.ZC
  par.mat$Beta.multinom.PL <- beta.pl.multinomial.old
  par.mat$Beta.NB <- beta.nb.old
  par.mat$Alpha.NB <- alpha.nb.old
  par.mat$Beta.PL <- as.numeric(beta.pl.old)
  par.mat$C <- c.pl

  res$par.mat <- par.mat

  res$par.all <- par.end.em

  res$resp <- upd.obj$resp

  res$log.lik.vec <- func.val.vec
  res$log.lik <- func.val.after.bfgs

  if (i.em < control$max.no.em.steps) {
    res$converge <- TRUE
  } else {
    res$converge <- FALSE
  }
  return(res)
}
