#' Control settings for evzinb() / evinb()
#'
#' Bundles the ~20 tuning arguments of the EM fitting routine into a single
#' object, in the style of \code{stats::glm.control()} or
#' \code{pscl::zeroinfl.control()}. Pass the result as the \code{control}
#' argument of \code{\link{evzinb}()} / \code{\link{evinb}()}.
#'
#' @param max.diff.par EM convergence tolerance: the algorithm has converged when
#'   the maximum absolute change in the parameter estimates falls below this.
#' @param max.no.em.steps Maximum number of EM steps.
#' @param max.no.em.steps.warmup Number of EM steps in each warm-up round.
#' @param c.lim \code{NULL} or a numeric vector of length 2. The candidate set for
#'   \eqn{C_{EV}} is the unique observed response values within this range.
#'   \code{NULL} (the default) uses a data-driven range (see Details).
#' @param prune.c.range \code{FALSE}, or a number in \[0, 1]: thin the candidate
#'   set to about \code{length(c.lim) * (1 - prune.c.range)} values.
#' @param max.upd.par.zc.multinomial,max.upd.par.pl.multinomial,max.upd.par.nb,max.upd.par.pl
#'   Maximum parameter-change step sizes for the zero-inflation, extreme-value
#'   inflation, count and Pareto components.
#' @param no.m.bfgs.steps.multinomial,no.m.bfgs.steps.nb,no.m.bfgs.steps.pl
#'   Number of BFGS steps per M-step for the multinomial, count and Pareto blocks.
#' @param pdf.pl.type Pareto density approximation: \code{"approx"} or \code{"exact"}.
#' @param eta.int Interval for the eta line search, a numeric vector of length 2.
#' @param init.Beta.multinom.ZC,init.Beta.multinom.PL,init.Beta.NB,init.Beta.PL
#'   Optional starting values for the component coefficient vectors (\code{NULL}
#'   starts at zero).
#' @param init.Alpha.NB Starting value for the negative-binomial dispersion.
#' @param init.C \code{NULL} or a starting value for \eqn{C_{EV}} within
#'   \code{c.lim}. \code{NULL} uses the median of the candidate set.
#' @param alpha_floor,coef_limit Thresholds for flagging a bootstrap replicate
#'   as \emph{degenerate} (\code{$degenerate}, \code{$degenerate_reason}), so it
#'   is excluded from bootstrap summaries by default (see
#'   \code{exclude_degenerate}) and counted in \code{glance()} /
#'   \code{failed_bootstraps()}. A replicate is degenerate if (in this order):
#'   its EM did not converge; any fitted coefficient in a linear predictor
#'   (\code{Beta.*}) or the negative-binomial dispersion is non-finite or
#'   exceeds \code{coef_limit} in absolute value; or its smallest fitted Pareto
#'   shape is non-finite or below \code{alpha_floor}. \code{alpha_floor}
#'   (default \code{0.001}) is deliberately low --- it catches an outright tail
#'   collapse (shape near 0), not a merely heavy tail. \code{coef_limit}
#'   (default \code{50}) is on the linear-predictor scale, where \code{50} is
#'   already extreme.
#'
#' @details When \code{c.lim = NULL}, \code{evzinb()} / \code{evinb()} choose a
#'   range from the data: the lower bound is the 90th percentile of the positive
#'   response values (rounded down to an observed value), the upper bound is the
#'   third-largest unique response value; if fewer than 10 unique positive values
#'   lie in that range the lower bound drops to the 75th percentile. The chosen
#'   range is printed with a message and stored in \code{object$control$c.lim}.
#'
#' @return A list of class \code{"evinf_control"}.
#' @export
#'
#' @examples
#' evinf_control(max.no.em.steps = 300, c.lim = c(20, 500))
evinf_control <- function(
  max.diff.par = 1e-2,
  max.no.em.steps = 500,
  max.no.em.steps.warmup = 5,
  c.lim = NULL,
  prune.c.range = FALSE,
  max.upd.par.zc.multinomial = 0.5,
  max.upd.par.pl.multinomial = 0.5,
  max.upd.par.nb = 0.5,
  max.upd.par.pl = 0.5,
  no.m.bfgs.steps.multinomial = 3,
  no.m.bfgs.steps.nb = 3,
  no.m.bfgs.steps.pl = 3,
  pdf.pl.type = c("approx", "exact"),
  eta.int = c(-1, 1),
  init.Beta.multinom.ZC = NULL,
  init.Beta.multinom.PL = NULL,
  init.Beta.NB = NULL,
  init.Beta.PL = NULL,
  init.Alpha.NB = 0.01,
  init.C = NULL,
  alpha_floor = 0.001,
  coef_limit = 50
) {
  pdf.pl.type <- match.arg(pdf.pl.type, c("approx", "exact"))
  control <- list(
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
    alpha_floor = alpha_floor,
    coef_limit = coef_limit
  )
  validate_evinf_control(control)
}

# The names evinf_control() accepts, i.e. the tuning arguments that evzinb() /
# evinb() still take individually (deprecated) for backwards compatibility.
evinf_control_args <- function() {
  setdiff(names(formals(evinf_control)), "")
}

validate_evinf_control <- function(control) {
  # Default any element added after this control object was created (e.g. one
  # restored from a model fitted with an older evinf).
  if (is.null(control$alpha_floor)) {
    control$alpha_floor <- 0.001
  }
  if (is.null(control$coef_limit)) {
    control$coef_limit <- 50
  }
  pos_scalar <- c(
    "max.diff.par", "max.no.em.steps", "max.no.em.steps.warmup",
    "max.upd.par.zc.multinomial", "max.upd.par.pl.multinomial",
    "max.upd.par.nb", "max.upd.par.pl", "no.m.bfgs.steps.multinomial",
    "no.m.bfgs.steps.nb", "no.m.bfgs.steps.pl", "init.Alpha.NB", "alpha_floor",
    "coef_limit"
  )
  for (nm in pos_scalar) {
    v <- control[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop("evinf_control(): `", nm, "` must be a single positive number.",
           call. = FALSE)
    }
  }

  if (!is.null(control$c.lim)) {
    cl <- control$c.lim
    if (!is.numeric(cl) || length(cl) != 2L || anyNA(cl) || cl[1] <= 0 ||
        cl[2] <= cl[1]) {
      stop("evinf_control(): `c.lim` must be NULL or a length-2 numeric with ",
           "0 < c.lim[1] < c.lim[2].", call. = FALSE)
    }
  }

  if (!identical(control$prune.c.range, FALSE)) {
    p <- control$prune.c.range
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p > 1) {
      stop("evinf_control(): `prune.c.range` must be FALSE or a number in [0, 1].",
           call. = FALSE)
    }
  }

  if (!is.numeric(control$eta.int) || length(control$eta.int) != 2L ||
      control$eta.int[2] <= control$eta.int[1]) {
    stop("evinf_control(): `eta.int` must be an increasing length-2 numeric.",
         call. = FALSE)
  }

  if (!is.null(control$init.C) &&
      (!is.numeric(control$init.C) || length(control$init.C) != 1L ||
       !is.finite(control$init.C) || control$init.C <= 0)) {
    stop("evinf_control(): `init.C` must be NULL or a single positive number.",
         call. = FALSE)
  }

  for (nm in c("init.Beta.multinom.ZC", "init.Beta.multinom.PL",
               "init.Beta.NB", "init.Beta.PL")) {
    v <- control[[nm]]
    if (!is.null(v) && (!is.numeric(v) || anyNA(v))) {
      stop("evinf_control(): `", nm, "` must be NULL or a numeric vector.",
           call. = FALSE)
    }
  }

  if (!control$pdf.pl.type %in% c("approx", "exact")) {
    stop("evinf_control(): `pdf.pl.type` must be \"approx\" or \"exact\".",
         call. = FALSE)
  }

  class(control) <- "evinf_control"
  control
}

# Merge any deprecated individual arguments the user passed to evzinb()/evinb()
# into the control object, warning once, then re-validate.
resolve_evinf_control <- function(control, call, env, fn = "evzinb") {
  if (!inherits(control, "evinf_control")) {
    stop("`control` must be created with evinf_control().", call. = FALSE)
  }
  supplied <- intersect(evinf_control_args(), names(call))
  for (a in supplied) {
    control[[a]] <- get(a, envir = env, inherits = FALSE)
  }
  if (length(supplied)) {
    warning(
      "The following arguments of ", fn,
      "() are deprecated; pass them through `control = evinf_control()`: ",
      paste(supplied, collapse = ", "), ".",
      call. = FALSE
    )
  }
  validate_evinf_control(control)
}

# Data-driven candidate range for C_EV (audit 4.5), used when control$c.lim is NULL.
default_c_lim <- function(y) {
  y <- y[is.finite(y)]
  pos <- y[y > 0]
  uniq_pos <- sort(unique(pos))
  uniq_y <- sort(unique(y))

  pick_lower <- function(q) {
    thr <- stats::quantile(pos, q, names = FALSE)
    below <- uniq_pos[uniq_pos <= thr]
    if (length(below)) max(below) else min(uniq_pos)
  }

  lower <- pick_lower(0.9)
  upper <- if (length(uniq_y) >= 3L) uniq_y[length(uniq_y) - 2L] else max(uniq_y)
  if (upper <= lower) {
    upper <- max(uniq_y)
  }
  in_range <- uniq_pos[uniq_pos >= lower & uniq_pos <= upper]
  if (length(in_range) < 10L) {
    lower <- pick_lower(0.75)
  }
  c(lower, upper)
}

# Resolve NULL c.lim / init.C against the response, message once, flag the default.
evinf_resolve_c <- function(control, y) {
  if (is.null(control$c.lim)) {
    control$c.lim <- default_c_lim(y)
    control$c_lim_default <- TRUE
    message(
      "evinf: using a data-driven candidate range for C_EV: [",
      control$c.lim[1], ", ", control$c.lim[2],
      "]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override."
    )
  } else {
    control$c_lim_default <- FALSE
  }
  if (is.null(control$init.C)) {
    cand <- sort(unique(y[is.finite(y)]))
    cand <- cand[cand >= control$c.lim[1] & cand <= control$c.lim[2]]
    control$init.C <- if (length(cand)) {
      cand[which.min(abs(cand - stats::median(cand)))]
    } else {
      mean(control$c.lim)
    }
  }
  control
}

#' @export
print.evinf_control <- function(x, ...) {
  cat("<evinf_control>\n")
  show <- c("max.diff.par", "max.no.em.steps", "max.no.em.steps.warmup",
            "c.lim", "prune.c.range", "pdf.pl.type", "init.Alpha.NB", "init.C")
  for (nm in show) {
    v <- x[[nm]]
    cat(sprintf("  %-24s %s\n", nm,
                if (is.null(v)) "NULL (data-driven)" else paste(v, collapse = ", ")))
  }
  invisible(x)
}
