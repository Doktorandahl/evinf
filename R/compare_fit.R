# audit 4.7 - paired-bootstrap fit comparison for evzinbcomp objects.

# Per-bootstrap fit statistics for the evinf model (constant npar, so AIC/BIC are
# recomputed from the trace log-likelihood rather than trusting a bootstrap
# object's own AIC, which for evinb still counts the fixed ZI coefficients).
evinf_boot_fit_stats <- function(comp) {
  ev <- comp$model
  npar <- length(ev$par.all)
  nobs <- nrow(ev$data$x.nb)
  lapply(ev$bootstraps, function(b) {
    if (inherits(b, "try-error")) {
      return(NULL)
    }
    ll <- b$log.lik
    c(logLik = ll, npar = npar,
      AIC = 2 * npar - 2 * ll, BIC = log(nobs) * npar - 2 * ll,
      rmse = NA_real_, rmsle = NA_real_)
  })
}

# Per-bootstrap fit statistics for a compared nbboot / zinbboot slot.
compared_boot_fit_stats <- function(slot) {
  lapply(slot$bootstraps, function(b) {
    if (inherits(b, "try-error") || is.null(b$fit_stats)) {
      return(NULL)
    }
    c(b$fit_stats,
      rmse = if (is.null(b$oob_rmse)) NA_real_ else b$oob_rmse,
      rmsle = if (is.null(b$oob_rmsle)) NA_real_ else b$oob_rmsle)
  })
}

#' Compare the extreme-value model against its plainer counterparts
#'
#' For each compared model (\code{nb}, \code{zinb}, and any winsorized /
#' razorized variants) and each metric, computes the paired bootstrap difference
#' \code{compared - evinf} (so a negative median favours the extreme-value
#' model), the proportion of bootstraps in which the evinf model is better, and
#' the number of bootstrap pairs where both fits succeeded.
#'
#' @param comp An \code{evzinbcomp} object from \code{\link{compare_models}()}.
#' @param metrics Which metrics: any of \code{"aic"}, \code{"bic"}, \code{"rmse"},
#'   \code{"rmsle"}. RMSE / RMSLE come from the out-of-bag predictions and need
#'   \code{\link{oob_evaluation}()} on the evinf model.
#' @param ... Unused.
#'
#' @return A tibble of class \code{evinf_compare_fit} with columns \code{model},
#'   \code{metric}, \code{median_difference}, \code{prop_evinf_better},
#'   \code{n_pairs}.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 10)
#' compare_fit(compare_models(model))
#' }
compare_fit <- function(comp, metrics = c("aic", "bic", "rmse", "rmsle"), ...) {
  if (!inherits(comp, "evzinbcomp")) {
    stop("`comp` must be an evzinbcomp object from compare_models().", call. = FALSE)
  }
  metrics <- match.arg(metrics, c("aic", "bic", "rmse", "rmsle"), several.ok = TRUE)
  metric_key <- c(aic = "AIC", bic = "BIC", rmse = "rmse", rmsle = "rmsle")

  ev_stats <- evinf_boot_fit_stats(comp)
  need_oob <- any(c("rmse", "rmsle") %in% metrics)
  if (need_oob) {
    ev_rmse <- suppressWarnings(oob_evaluation(comp$model, metric = "rmse"))
    ev_rmsle <- suppressWarnings(oob_evaluation(comp$model, metric = "rmsle"))
    for (k in seq_along(ev_stats)) {
      if (!is.null(ev_stats[[k]])) {
        ev_stats[[k]]["rmse"] <- ev_rmse[k]
        ev_stats[[k]]["rmsle"] <- ev_rmsle[k]
      }
    }
  }

  model_slots <- setdiff(names(comp), c("model", "evzinb"))
  rows <- list()
  for (slot in model_slots) {
    cmp_stats <- compared_boot_fit_stats(comp[[slot]])
    for (m in metrics) {
      key <- metric_key[[m]]
      diffs <- vapply(seq_along(ev_stats), function(k) {
        a <- ev_stats[[k]]; b <- cmp_stats[[k]]
        if (is.null(a) || is.null(b)) NA_real_ else unname(b[key] - a[key])
      }, numeric(1))
      diffs <- diffs[is.finite(diffs)]
      rows[[length(rows) + 1L]] <- tibble::tibble(
        model = slot, metric = m,
        median_difference = if (length(diffs)) stats::median(diffs) else NA_real_,
        prop_evinf_better = if (length(diffs)) mean(diffs > 0) else NA_real_,
        n_pairs = length(diffs)
      )
    }
  }
  out <- dplyr::bind_rows(rows)
  class(out) <- c("evinf_compare_fit", class(out))
  out
}

#' @export
print.evinf_compare_fit <- function(x, ...) {
  cat("Paired bootstrap fit comparison (compared - evinf)\n")
  cat("  negative median favours the extreme-value model\n\n")
  df <- as.data.frame(x)
  df$median_difference <- signif(df$median_difference, 4)
  df$prop_evinf_better <- round(df$prop_evinf_better, 3)
  print(df, row.names = FALSE)
  invisible(x)
}

# Long tibble of every per-bootstrap paired difference (compared - evinf).
evinf_compare_fit_long <- function(comp, metrics) {
  ev_stats <- evinf_boot_fit_stats(comp)
  if (any(c("rmse", "rmsle") %in% metrics)) {
    r <- suppressWarnings(oob_evaluation(comp$model, metric = "rmse"))
    rl <- suppressWarnings(oob_evaluation(comp$model, metric = "rmsle"))
    for (k in seq_along(ev_stats)) if (!is.null(ev_stats[[k]])) {
      ev_stats[[k]]["rmse"] <- r[k]; ev_stats[[k]]["rmsle"] <- rl[k]
    }
  }
  key <- c(aic = "AIC", bic = "BIC", rmse = "rmse", rmsle = "rmsle")
  slots <- setdiff(names(comp), c("model", "evzinb"))
  dplyr::bind_rows(lapply(slots, function(slot) {
    cmp <- compared_boot_fit_stats(comp[[slot]])
    dplyr::bind_rows(lapply(metrics, function(m) {
      d <- vapply(seq_along(ev_stats), function(k) {
        a <- ev_stats[[k]]; b <- cmp[[k]]
        if (is.null(a) || is.null(b)) NA_real_ else unname(b[key[[m]]] - a[key[[m]]])
      }, numeric(1))
      tibble::tibble(model = slot, metric = m, difference = d[is.finite(d)])
    }))
  }))
}

#' @rdname compare_fit
#' @param x An \code{evzinbcomp} object.
#' @param metrics Metrics to show.
#' @export
plot.evzinbcomp <- function(x, metrics = c("aic", "bic", "rmse", "rmsle"), ...) {
  rlang::check_installed("ggplot2", "for plot.evzinbcomp()")
  metrics <- match.arg(metrics, c("aic", "bic", "rmse", "rmsle"), several.ok = TRUE)
  long <- evinf_compare_fit_long(x, metrics)
  ggplot2::ggplot(long, ggplot2::aes(.data$difference)) +
    ggplot2::geom_density(fill = "grey80", colour = NA) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$model),
                        cols = ggplot2::vars(.data$metric), scales = "free") +
    ggplot2::labs(x = "compared - evinf  (negative favours evinf)",
                  y = "bootstrap density") +
    ggplot2::theme_minimal()
}

#' Stacked tidy / glance for a model comparison
#'
#' @param x An \code{evzinbcomp} object.
#' @param ... Passed to the per-model \code{tidy()} / \code{glance()}.
#' @return A tibble with a leading \code{model} column.
#' @name evzinbcomp-tidiers
#' @export
tidy.evzinbcomp <- function(x, ...) {
  slots <- setdiff(names(x), "evzinb")
  slots <- c("model", setdiff(slots, "model"))
  dplyr::bind_rows(lapply(slots, function(s) {
    td <- tryCatch(suppressWarnings(generics::tidy(x[[s]], ...)), error = function(e) NULL)
    if (is.null(td)) NULL else dplyr::mutate(td, model = s, .before = 1)
  }))
}

#' @rdname evzinbcomp-tidiers
#' @export
glance.evzinbcomp <- function(x, ...) {
  slots <- setdiff(names(x), "evzinb")
  slots <- c("model", setdiff(slots, "model"))
  dplyr::bind_rows(lapply(slots, function(s) {
    gl <- tryCatch(suppressWarnings(generics::glance(x[[s]], ...)), error = function(e) NULL)
    if (is.null(gl)) NULL else dplyr::mutate(gl, model = s, .before = 1)
  }))
}
