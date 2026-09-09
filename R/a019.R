## Internal functions for the EVZINB and EVINB functions
zerinfl.nb.pl.reg.cond.c.fun <- function(y, x.obj, ini.val, control) {
  #Initialize parameters and data
  x.multinom.zc <- x.obj$X.multinom.ZC
  x.multinom.pl <- x.obj$X.multinom.PL

  x.nb <- x.obj$X.NB
  x.pl <- x.obj$X.PL

  n <- max(nrow(x.nb), nrow(x.pl), nrow(x.multinom.zc), nrow(x.multinom.pl))

  offset.nb <- x.obj$offset.nb
  if (is.null(offset.nb)) offset.nb <- rep(0, n)  # audit 4.4

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
    y,
    offset.nb
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
      control$no.m.bfgs.steps.nb,
      offset.nb
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
            y,
            offset.nb
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
            y,
            offset.nb
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
            y,
            offset.nb
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
          y,
          offset.nb
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
            y,
            offset.nb
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
          y,
          offset.nb
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
            y,
            offset.nb
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
          y,
          offset.nb
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
      y,
      offset.nb
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

  offset.nb <- x.obj$offset.nb
  if (is.null(offset.nb)) offset.nb <- rep(0, n)  # audit 4.4

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

  offset.nb <- x.obj$offset.nb
  if (is.null(offset.nb)) offset.nb <- rep(0, n)  # audit 4.4

  #Determine which values of c to be investigated
  c.range <- unique(sort(y))[
    unique(sort(y)) >= control$c.lim[1] & unique(sort(y)) <= control$c.lim[2]
  ]

  if (length(c.range) > 100 && control$prune.c.range == FALSE) {
    warning(
      "The c.range contains more than 100 values. If the estimation is slow consider reducing the c-range with the c.lim argument of the function, or pruning the c-range using the prune.c.range argument."
    )
  }
  if (is.numeric(control$prune.c.range)) {
    if (control$prune.c.range < 0 || control$prune.c.range > 1) {
      stop("The prune.c.range argument must be FALSE or be between 0 and 1")
    }
    gaps <- diff(c.range)
    sample.size <- ceiling(length(c.range) * (1 - control$prune.c.range)) - 2
    keep.indicies <- c(
      1,
      sample(
        seq(2, length(c.range) - 1),
        size = sample.size,
        prob = gaps[-c(length(c.range) - 1)] / sum(gaps)
      ),
      length(c.range)
    )
    c.range <- c.range[sort(keep.indicies)]
    if (length(c.range) > 100) {
      warning(
        "The pruned c.range still contains more than 100 values. If the estimation is slow consider reducing the c-range with the c.lim argument of the function, or increasing the prune.c.range argument."
      )
    }
  }

  ###When c is not well known we only run short versions of zerinfl.nb.pl.reg.cond.c.fun() in order to save time
  control.warmup <- control
  control.warmup$max.no.em.steps <- control$max.no.em.steps.warmup

  ###Initializing the warm-up phase when the update stops after control$max.no.em.steps.warmup EM steps
  #If the update of c is less than 1, stop the warmup
  c.abs.diff <- 100
  log.lik.vec.all <- NULL
  c_trace <- numeric(0)  # audit 4.5: sequence of chosen c values across iterations
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
        y,
        offset.nb
      )
    }

    c.pl.new <- c.range[which(log.lik.vec == max(log.lik.vec))]
    c_trace <- c(c_trace, c.pl.new)

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
        y,
        offset.nb
      )
    }

    c.pl.new <- c.range[which(log.lik.vec == max(log.lik.vec))]
    c_trace <- c(c_trace, c.pl.new)

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
  # audit 4.13: vectorised equivalent of the former per-observation loop.
  mu.nb.vec <- as.numeric(exp(x.nb.extended %*% prel.val$Beta.NB + offset.nb))
  alpha.pl.vec <- as.numeric(exp(x.pl.extended %*% prel.val$Beta.PL))
  exp.E.log.y <- c.pl.new * exp(1 / alpha.pl.vec)
  E.inv.y <- alpha.pl.vec / (c.pl.new * (alpha.pl.vec + 1))
  mean.pl.vec <- ifelse(
    alpha.pl.vec > 1,
    alpha.pl.vec * c.pl.new / (alpha.pl.vec - 1),
    NA_real_
  )
  median.pl.vec <- c.pl.new * 2^(1 / alpha.pl.vec)
  y.hat.plmedian <- props.old[, 2] * mu.nb.vec + props.old[, 3] * median.pl.vec
  y.hat.plmean <- props.old[, 2] * mu.nb.vec + props.old[, 3] * mean.pl.vec
  y.hat.plexpElogy <- props.old[, 2] * mu.nb.vec + props.old[, 3] * exp.E.log.y
  y.hat.pl.E.inv.y <- props.old[, 2] * mu.nb.vec + props.old[, 3] / E.inv.y

  par.all <- c(
    final.val$Beta.multinom.ZC,
    final.val$Beta.multinom.PL,
    final.val$Beta.NB,
    final.val$Alpha.NB,
    final.val$Beta.PL,
    final.val$C
  )
  n.par <- length(par.all)

  # audit 2.13: the trace maximum (func.val) is not guaranteed to be the
  # log-likelihood at the returned parameters. Recompute it there and, if they
  # disagree, use the value at the returned parameters for log.lik / AIC / BIC.
  ll.at.par <- log_lik_fun(
    final.val$Beta.multinom.ZC,
    final.val$Beta.multinom.PL,
    final.val$Beta.NB,
    final.val$Alpha.NB,
    final.val$Beta.PL,
    final.val$C,
    x.multinom.zc.extended,
    x.multinom.pl.extended,
    x.nb.extended,
    x.pl.extended,
    y,
    offset.nb
  )
  # audit R0.3: recompute silently and record whether the correction happened;
  # run_*() decides whether to tell the user (never inside a bootstrap).
  loglik_recomputed <- isTRUE(is.finite(ll.at.par) &&
                                abs(ll.at.par - func.val) > 1e-6)
  if (loglik_recomputed) {
    func.val <- ll.at.par
  }

  BIC <- log(n) * n.par - 2 * func.val
  AIC <- 2 * n.par - 2 * func.val

  #Gather results from estimation
  out <- list()

  out$control <- control
  out$par.mat <- final.val
  out$log.lik.vec.all <- log.lik.vec.all
  # audit 4.5: log-likelihood profile over the final candidate set, and the
  # sequence of chosen c values across EM iterations.
  out$c_profile <- data.frame(c = c.range, loglik = log.lik.vec)
  out$c_trace <- c_trace
  out$log.lik <- func.val
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
  out$loglik_recomputed <- loglik_recomputed

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

plinfl.nb.regression.fun <- function(y, x.obj, ini.val, control) {
  #Initialize data and parameters

  x.multinom.zc <- x.obj$X.multinom.ZC
  x.multinom.pl <- x.obj$X.multinom.PL
  x.nb <- x.obj$X.NB
  x.pl <- x.obj$X.PL

  n <- max(nrow(x.nb), nrow(x.pl), nrow(x.multinom.zc), nrow(x.multinom.pl))

  offset.nb <- x.obj$offset.nb
  if (is.null(offset.nb)) offset.nb <- rep(0, n)  # audit 4.4

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

  offset.nb <- x.obj$offset.nb
  if (is.null(offset.nb)) offset.nb <- rep(0, n)  # audit 4.4

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
  c_trace <- numeric(0)  # audit 4.5: sequence of chosen c values across iterations
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
        y,
        offset.nb
      )
    }

    c.pl.new <- c.range[which(log.lik.vec == max(log.lik.vec))]
    c_trace <- c(c_trace, c.pl.new)

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
        y,
        offset.nb
      )
    }

    c.pl.new <- c.range[which(log.lik.vec == max(log.lik.vec))]
    c_trace <- c(c_trace, c.pl.new)

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
  # audit 4.13: vectorised equivalent of the former per-observation loop.
  mu.nb.vec <- as.numeric(exp(x.nb.extended %*% prel.val$Beta.NB + offset.nb))
  alpha.pl.vec <- as.numeric(exp(x.pl.extended %*% prel.val$Beta.PL))
  exp.E.log.y <- c.pl.new * exp(1 / alpha.pl.vec)
  E.inv.y <- alpha.pl.vec / (c.pl.new * (alpha.pl.vec + 1))
  mean.pl.vec <- ifelse(
    alpha.pl.vec > 1,
    alpha.pl.vec * c.pl.new / (alpha.pl.vec - 1),
    NA_real_
  )
  median.pl.vec <- c.pl.new * 2^(1 / alpha.pl.vec)
  y.hat.plmedian <- props.old[, 2] * mu.nb.vec + props.old[, 3] * median.pl.vec
  y.hat.plmean <- props.old[, 2] * mu.nb.vec + props.old[, 3] * mean.pl.vec
  y.hat.plexpElogy <- props.old[, 2] * mu.nb.vec + props.old[, 3] * exp.E.log.y
  y.hat.pl.E.inv.y <- props.old[, 2] * mu.nb.vec + props.old[, 3] / E.inv.y

  par.all <- c(
    final.val$Beta.multinom.ZC,
    final.val$Beta.multinom.PL,
    final.val$Beta.NB,
    final.val$Alpha.NB,
    final.val$Beta.PL,
    final.val$C
  )
  n.par <- length(par.all)

  # audit 2.13: the trace maximum (func.val) is not guaranteed to be the
  # log-likelihood at the returned parameters. Recompute it there and, if they
  # disagree, use the value at the returned parameters for log.lik / AIC / BIC.
  ll.at.par <- log_lik_fun(
    final.val$Beta.multinom.ZC,
    final.val$Beta.multinom.PL,
    final.val$Beta.NB,
    final.val$Alpha.NB,
    final.val$Beta.PL,
    final.val$C,
    x.multinom.zc.extended,
    x.multinom.pl.extended,
    x.nb.extended,
    x.pl.extended,
    y,
    offset.nb
  )
  # audit R0.3: recompute silently and record whether the correction happened;
  # run_*() decides whether to tell the user (never inside a bootstrap).
  loglik_recomputed <- isTRUE(is.finite(ll.at.par) &&
                                abs(ll.at.par - func.val) > 1e-6)
  if (loglik_recomputed) {
    func.val <- ll.at.par
  }

  BIC <- log(n) * n.par - 2 * func.val
  AIC <- 2 * n.par - 2 * func.val

  #Gather results from estimation
  out <- list()

  out$control <- control
  out$par.mat <- final.val
  out$log.lik.vec.all <- log.lik.vec.all
  # audit 4.5: log-likelihood profile over the final candidate set, and the
  # sequence of chosen c values across EM iterations.
  out$c_profile <- data.frame(c = c.range, loglik = log.lik.vec)
  out$c_trace <- c_trace
  out$log.lik <- func.val
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
  out$loglik_recomputed <- loglik_recomputed

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

  offset.nb <- x.obj$offset.nb
  if (is.null(offset.nb)) offset.nb <- rep(0, n)  # audit 4.4

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
    y,
    offset.nb
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
      control$no.m.bfgs.steps.nb,
      offset.nb
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
            y,
            offset.nb
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
            y,
            offset.nb
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
            y,
            offset.nb
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
          y,
          offset.nb
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
            y,
            offset.nb
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
          y,
          offset.nb
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
      y,
      offset.nb
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
