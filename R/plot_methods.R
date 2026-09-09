# audit 4.6 - plot methods for evzinb / evinb models.

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
#'   \code{"ppc"} - posterior predictive check: observed vs simulated quantiles
#'     and a binned observed-vs-expected frequency panel.
#' @param variable Covariate to vary (required for \code{"states"} and
#'   \code{"prediction"}).
#' @param quantiles Quantiles to draw for \code{type = "prediction"}.
#' @param ... Passed to \code{\link{predict_grid}}.
#'
#' @return A \code{ggplot} object (or, for \code{"ppc"}, a patchwork-free
#'   two-panel plot built with \code{ggplot2}).
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' plot(model, type = "coefficients")
#' plot(model, type = "prediction", variable = "x1")
#' }
plot.evzinb <- function(x, type = c("states", "prediction", "coefficients", "ppc"),
                        variable = NULL, quantiles = c(0.5, 0.95), ...) {
  rlang::check_installed("ggplot2", "for plot.evzinb() / plot.evinb()")
  type <- match.arg(type)

  switch(type,
    states = evinf_plot_states(x, variable, ...),
    prediction = evinf_plot_prediction(x, variable, quantiles, ...),
    coefficients = evinf_plot_coefficients(x),
    ppc = evinf_plot_ppc(x)
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
    ggplot2::scale_y_continuous(transform = "log1p") +
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
    ggplot2::scale_y_continuous(transform = "log1p") +
    ggplot2::labs(x = "y (binned)", y = "count (log1p scale)", fill = NULL,
                  title = "Posterior predictive check: observed vs expected frequencies") +
    ggplot2::theme_minimal()
}
