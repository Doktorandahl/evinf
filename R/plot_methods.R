# audit 4.6 - plot methods for evzinb / evinb models.

# round8 0.1: ggplot2 >= 3.5.0 renamed scale_*_continuous()'s `trans` argument
# to `transform`; `trans` still works but emits a deprecation message on the
# newer versions, and `transform` doesn't exist at all before 3.5.0. Keep the
# package usable across both by picking the argument name at call time.
# `ggplot2_version` is a parameter (not read internally) so tests can exercise
# both branches without needing two ggplot2 installations.
evinf_scale_log1p <- function(axis = c("x", "y"),
                              ggplot2_version = utils::packageVersion("ggplot2")) {
  axis <- match.arg(axis)
  scale_fun <- if (axis == "x") {
    ggplot2::scale_x_continuous
  } else {
    ggplot2::scale_y_continuous
  }
  if (ggplot2_version >= "3.5.0") {
    scale_fun(transform = "log1p")
  } else {
    scale_fun(trans = "log1p")
  }
}

#' Plots for evzinb / evinb models
#'
#' @param x A fitted \code{evzinb} / \code{evinb} model.
#' @param type
#'   \code{"states"} - prior state probabilities over \code{variable}, faceted by
#'     state;
#'   \code{"prediction"} - harmonic-mean prediction with a bootstrap ribbon plus
#'     the requested quantiles over \code{variable} (log1p y axis);
#'   \code{"coefficients"} - bootstrap densities of each coefficient, faceted by
#'     component and term, with the point estimate marked;
#'   \code{"ppc"} - posterior predictive check: a binned observed-vs-expected
#'     frequency panel;
#'   \code{"ppc_quantiles"} - posterior predictive check on the tails: the
#'     observed sample quantiles for probabilities
#'     \code{c(0.5, 0.75, 0.9, 0.95, 0.99, 0.999)} against the median simulated
#'     quantile (with a 5-95\% band across simulations) on \code{log1p} axes,
#'     with a 45-degree reference line;
#'   \code{"trace"} (round10 G.2, audit §5.10) - the log-likelihood
#'     (\code{$loglik_trace}) and \eqn{C_{EV}} (\code{$c_trace}) traces, in two
#'     stacked panels with the warm-up/convergence boundary marked (dashed
#'     line). With \code{evinf_control(n_starts > 1)}, every perturbed
#'     start's trace is overlaid in grey behind the winning start's, in black
#'     (traces are only stored per-start when \code{n_starts > 1}, to keep a
#'     single-start object the same size as before this type existed).
#' @param variable Covariate to vary (required for \code{"states"} and
#'   \code{"prediction"}).
#' @param quantiles Quantiles to draw for \code{type = "prediction"}.
#' @param ... Passed to \code{\link{predict_grid}}.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' plot(model, type = "coefficients")
#' plot(model, type = "prediction", variable = "x1")
#' }
#'
#' \dontrun{
#' data(hks)
#' hks_mod <- evzinb(osvAll ~ troopLag + lntpop + brv_AllLag_log,
#'                   data = hks, n_bootstraps = 5, multicore = FALSE)
#' plot(hks_mod, type = "prediction", variable = "troopLag")
#' }
plot.evzinb <- function(x, type = c("states", "prediction", "coefficients",
                                    "ppc", "ppc_quantiles", "trace"),
                        variable = NULL, quantiles = c(0.5, 0.95), ...) {
  rlang::check_installed("ggplot2", "for plot.evzinb() / plot.evinb()")
  type <- match.arg(type)

  switch(type,
    states = evinf_plot_states(x, variable, ...),
    prediction = evinf_plot_prediction(x, variable, quantiles, ...),
    coefficients = evinf_plot_coefficients(x),
    ppc = evinf_plot_ppc(x),
    ppc_quantiles = evinf_plot_ppc_quantiles(x),
    trace = evinf_plot_trace(x)
  )
}

#' @rdname plot.evzinb
#' @export
plot.evinb <- plot.evzinb


evinf_plot_states <- function(x, variable, ...) {
  if (is.null(variable)) {
    stop("`variable` is required for type = \"states\".", call. = FALSE)
  }
  g <- predict_grid(x, variable, type = "states", ...)
  g$state <- factor(sub("^pr_", "", g$type), levels = c("zero", "count", "evi"))
  ggplot2::ggplot(g, ggplot2::aes(.data$value, .data$estimate)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(.data$state)) +
    ggplot2::labs(x = variable, y = "prior state probability") +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal()
}

evinf_plot_prediction <- function(x, variable, quantiles, ...) {
  if (is.null(variable)) {
    stop("`variable` is required for type = \"prediction\".", call. = FALSE)
  }
  has_boot <- !is.null(x$bootstraps)
  harm <- predict_grid(x, variable, type = "harmonic", confint = has_boot, ...)
  harm$type <- "harmonic mean"
  qs <- lapply(quantiles, function(q) {
    d <- predict_grid(x, variable, type = "quantile", quantile = q, ...)
    d$type <- paste0("q", 100 * q)
    d
  })
  dat <- dplyr::bind_rows(c(list(harm), qs))

  p <- ggplot2::ggplot(dat, ggplot2::aes(.data$value, .data$estimate,
                                         colour = .data$type, fill = .data$type))
  if (has_boot) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high),
      alpha = 0.15, colour = NA
    )
  }
  p +
    ggplot2::geom_line() +
    evinf_scale_log1p("y") +
    ggplot2::labs(x = variable, y = "predicted y (log1p scale)",
                  colour = NULL, fill = NULL) +
    ggplot2::theme_minimal()
}

evinf_plot_coefficients <- function(x) {
  if (is.null(x$bootstraps)) {
    stop("type = \"coefficients\" needs a model fitted with bootstrap = TRUE.",
         call. = FALSE)
  }
  boot <- suppressWarnings(coefficient_extractor(x, component = "all"))
  long <- tidyr::pivot_longer(
    boot, cols = -dplyr::all_of(".component"),
    names_to = "term", values_to = "value"
  )
  long <- long[!is.na(long$value), ]
  long$component <- factor(long$.component, levels = c("count", "zero", "evi", "pareto"))

  est <- tibble::enframe(coef(x, "all"), name = "key", value = "estimate")
  est$component <- sub("_.*$", "", est$key)
  est$term <- sub("^[^_]*_", "", est$key)
  est <- est[est$component %in% c("count", "zero", "evi", "pareto"), ]
  est$component <- factor(est$component, levels = c("count", "zero", "evi", "pareto"))

  ggplot2::ggplot(long, ggplot2::aes(.data$value)) +
    ggplot2::geom_density(fill = "grey80", colour = NA) +
    ggplot2::geom_vline(data = est, ggplot2::aes(xintercept = .data$estimate),
                        linetype = "dashed") +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$component),
                        cols = ggplot2::vars(.data$term), scales = "free") +
    ggplot2::labs(x = "coefficient", y = "bootstrap density") +
    ggplot2::theme_minimal()
}

evinf_plot_ppc <- function(x) {
  y <- x$data$y
  sims <- simulate(x, nsim = 100)

  # observed vs expected frequencies in coarse bins
  brks <- c(-1, 0, 1, 4, 9, 49, 99, 999, Inf)
  labs <- c("0", "1", "2-4", "5-9", "10-49", "50-99", "100-999", "1000+")
  obs <- as.numeric(table(cut(y, brks, labels = labs)))
  exp_counts <- rowMeans(vapply(sims, function(s) {
    as.numeric(table(cut(s, brks, labels = labs)))
  }, numeric(length(labs))))
  freq <- tibble::tibble(
    bin = factor(rep(labs, 2), levels = labs),
    count = c(obs, exp_counts),
    source = rep(c("observed", "expected"), each = length(labs))
  )

  ggplot2::ggplot(freq, ggplot2::aes(.data$bin, .data$count, fill = .data$source)) +
    ggplot2::geom_col(position = "dodge") +
    evinf_scale_log1p("y") +
    ggplot2::labs(x = "y (binned)", y = "count (log1p scale)", fill = NULL,
                  title = "Posterior predictive check: observed vs expected frequencies") +
    ggplot2::theme_minimal()
}

evinf_plot_ppc_quantiles <- function(x) {
  y <- x$data$y
  probs <- c(0.5, 0.75, 0.9, 0.95, 0.99, 0.999)
  sims <- simulate(x, nsim = 100)

  obs_q <- stats::quantile(y, probs, names = FALSE)
  # sim_q: one column of quantiles per simulation
  sim_q <- vapply(sims, function(s) stats::quantile(s, probs, names = FALSE),
                  numeric(length(probs)))

  dat <- tibble::tibble(
    prob = probs,
    observed = obs_q,
    simulated = apply(sim_q, 1L, stats::median),
    lo = apply(sim_q, 1L, stats::quantile, 0.05, names = FALSE),
    hi = apply(sim_q, 1L, stats::quantile, 0.95, names = FALSE)
  )

  ggplot2::ggplot(dat, ggplot2::aes(.data$observed, .data$simulated)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                            colour = "grey50") +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = paste0(100 * .data$prob, "%")),
                       vjust = -0.8, size = 3) +
    evinf_scale_log1p("x") +
    evinf_scale_log1p("y") +
    ggplot2::labs(x = "observed quantile (log1p scale)",
                  y = "simulated quantile (log1p scale)",
                  title = paste("Posterior predictive check: tail quantiles",
                                "(5-95% band across simulations)")) +
    ggplot2::theme_minimal()
}

# round10 G.2 (audit §5.10): the log-likelihood ($loglik_trace) and C_EV
# ($c_trace) traces, stacked via facet_wrap() (one row per panel) rather
# than a second plotting package -- no new dependency. With n_starts > 1
# every perturbed start's own trace ($starts$loglik_trace/$c_trace, stored
# only in that case -- see evinf_run_starts()) is overlaid in grey behind
# the winning start's, in black.
evinf_plot_trace <- function(x) {
  if (is.null(x$loglik_trace) || is.null(x$c_trace)) {
    stop("type = \"trace\" needs a fit with $loglik_trace / $c_trace.",
         call. = FALSE)
  }
  trace_df <- function(ll, ct, label) {
    dplyr::bind_rows(
      tibble::tibble(step = seq_along(ll), value = ll,
                     panel = "log-likelihood", start = label),
      tibble::tibble(step = seq_along(ct), value = ct,
                     panel = "C_EV", start = label)
    )
  }
  dat <- trace_df(x$loglik_trace, x$c_trace, "best")

  if (!is.null(x$starts)) {
    others <- dplyr::bind_rows(lapply(seq_len(nrow(x$starts)), function(i) {
      trace_df(x$starts$loglik_trace[[i]], x$starts$c_trace[[i]],
              paste0("start ", i))
    }))
    dat <- dplyr::bind_rows(others, dat)
  }
  dat$panel <- factor(dat$panel, levels = c("log-likelihood", "C_EV"))

  vlines <- tibble::tibble(
    panel = factor(c("log-likelihood", "C_EV"), levels = c("log-likelihood", "C_EV")),
    xintercept = c(x$n_loglik_warmup %||% NA_integer_, x$n_c_iter_warmup %||% NA_integer_)
  )
  vlines <- vlines[is.finite(vlines$xintercept) & vlines$xintercept > 0, ]

  p <- ggplot2::ggplot(dat, ggplot2::aes(.data$step, .data$value, group = .data$start))
  if (!is.null(x$starts)) {
    p <- p + ggplot2::geom_line(
      data = dat[dat$start != "best", ], colour = "grey70", alpha = 0.6
    )
  }
  p <- p + ggplot2::geom_line(data = dat[dat$start == "best", ], colour = "black")
  if (nrow(vlines)) {
    p <- p + ggplot2::geom_vline(
      data = vlines, ggplot2::aes(xintercept = .data$xintercept),
      linetype = "dashed", colour = "grey40"
    )
  }
  p +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "EM step", y = NULL,
                  title = "EM trace (dashed: warm-up / convergence boundary)") +
    ggplot2::theme_minimal()
}
