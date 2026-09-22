# Diagnostics: a single convergence/health report (round10 G.4, audit §5.10).

#' Diagnostic report for a fitted evzinb / evinb model
#'
#' A single report gathering the convergence and health checks otherwise
#' spread across \code{print()}, \code{glance()} and \code{failed_bootstraps()}:
#' EM and \eqn{C_{EV}}-profile convergence, whether \eqn{C_{EV}} landed on a
#' candidate-grid boundary (full sample and across bootstrap replicates), the
#' smallest fitted Pareto shape against \code{alpha_pl_floor}, failed and
#' degenerate bootstrap replicates, start agreement (see \code{n_starts} in
#' \code{\link{evinf_control}}), and, for the block bootstrap schemes, the
#' out-of-bag row fraction.
#'
#' @param object A fitted \code{evzinb} / \code{evinb} model.
#' @return A classed tibble (\code{"evinf_check"}) with columns \code{check},
#'   \code{status} (\code{"ok"}, \code{"note"} or \code{"warning"}) and
#'   \code{detail} (a one-line description, with a suggestion for anything
#'   not \code{"ok"}).
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' check_evinf(model)
#' }
check_evinf <- function(object) {
  if (!inherits(object, c("evzinb", "evinb"))) {
    stop("`object` must be a fitted evzinb / evinb model.", call. = FALSE)
  }
  rows <- list()
  add <- function(check, status, detail) {
    rows[[length(rows) + 1L]] <<- tibble::tibble(check = check, status = status, detail = detail)
  }

  # --- EM / C_EV-profile convergence ----------------------------------------
  add("EM converged", if (isTRUE(object$converge)) "ok" else "warning",
      if (isTRUE(object$converge)) {
        "The EM algorithm converged."
      } else {
        "The EM algorithm did not converge -- see $converge, $c_converged and raise max.diff.par/max.no.em.steps/max.c.iter in evinf_control() if needed."
      })

  c_converged <- if (is.null(object$c_converged)) NA else isTRUE(object$c_converged)
  add("C_EV profile converged (convergence phase)",
      if (is.na(c_converged)) "note" else if (c_converged) "ok" else "warning",
      if (is.na(c_converged)) {
        "Not recorded on this fit (fitted before this field existed)."
      } else if (c_converged) {
        "The C_EV profile settled within max.c.iter."
      } else {
        "The C_EV profile did not settle within max.c.iter (kept oscillating between candidates) -- consider raising max.c.iter in evinf_control()."
      })

  c_warmup_capped <- isTRUE(object$c_warmup_capped)
  add("C_EV profile converged (warm-up phase)",
      if (c_warmup_capped) "note" else "ok",
      if (c_warmup_capped) {
        "The warm-up phase hit max.c.iter without settling -- usually harmless (a short exploratory phase) if the convergence phase above settled."
      } else {
        "The warm-up phase settled within max.c.iter."
      })

  # --- C_EV on a candidate-grid boundary ------------------------------------
  on_boundary <- !is.null(object$c_profile) &&
    (object$coef$C <= min(object$c_profile$c) || object$coef$C >= max(object$c_profile$c))
  n_c_bnd <- object$n_c_on_boundary
  n_boot_total <- if (is.null(object$bootstraps)) NA_integer_ else length(object$bootstraps)
  boundary_detail <- if (on_boundary) {
    paste0("The fitted C_EV (", object$coef$C, ") lies on the edge of the candidate range -- ",
          "consider widening c.lim in evinf_control().")
  } else {
    "The fitted C_EV lies inside the candidate range."
  }
  if (isTRUE(n_c_bnd > 0)) {
    boundary_detail <- paste0(boundary_detail, " ", n_c_bnd, " of ", n_boot_total,
                              " bootstrap replicates also landed on a boundary.")
  }
  add("C_EV on candidate-grid boundary", if (on_boundary) "warning" else "ok", boundary_detail)

  # --- Fitted Pareto shape vs the floor -------------------------------------
  min_alpha_pl <- if (is.null(object$fitted$alpha.pl)) NA_real_ else min(object$fitted$alpha.pl)
  alpha_pl_floor <- object$control$alpha_pl_floor %||% 0.01
  collapsed <- isTRUE(is.finite(min_alpha_pl) && min_alpha_pl < alpha_pl_floor)
  add("Fitted Pareto shape (alpha_pl) not collapsed",
      if (is.na(min_alpha_pl)) "note" else if (collapsed) "warning" else "ok",
      if (is.na(min_alpha_pl)) {
        "Not available on this fit."
      } else if (collapsed) {
        paste0("min(alpha_pl) = ", signif(min_alpha_pl, 3), " is below alpha_pl_floor = ",
              alpha_pl_floor, "; predictions involving it are floored. See glance()$min_alpha_pl and predict(type = 'explog')'s warning.")
      } else {
        paste0("min(alpha_pl) = ", signif(min_alpha_pl, 3), ", above the floor.")
      })

  # --- Bootstrap replicates -------------------------------------------------
  boot <- evinf_boot_counts(object)
  if (is.na(boot$n_bootstraps)) {
    add("Bootstrap replicates usable", "note", "This model was fitted without bootstrapping.")
  } else {
    n_bad <- boot$n_failed_bootstraps + boot$n_degenerate_bootstraps
    add("Bootstrap replicates usable",
        if (n_bad == 0) "ok" else "note",
        paste0(boot$n_bootstraps, " usable of ", boot$n_bootstraps + n_bad, " (",
              boot$n_failed_bootstraps, " failed, ", boot$n_degenerate_bootstraps,
              " degenerate)",
              if (n_bad > 0) " -- see failed_bootstraps() for the details." else "."))
  }

  # --- Multiple starts -------------------------------------------------------
  if (!is.null(object$starts)) {
    st <- evinf_n_starts_at_best(object)
    frac_at_best <- st$n_starts_at_best / st$n_starts
    add("Multiple starts agree",
        if (frac_at_best >= 0.5) "ok" else "note",
        paste0(st$n_starts_at_best, " of ", st$n_starts,
              " starts reached the best log-likelihood (within 1e-4).",
              if (frac_at_best < 0.5) {
                " Most starts landed elsewhere -- the likelihood surface may have multiple local optima; consider more starts (n_starts in evinf_control())."
              } else {
                ""
              }))
  }

  # --- Out-of-bag fraction, for the block schemes ---------------------------
  if (isTRUE(object$bootstrap_scheme %in% c("moving_block", "stationary"))) {
    oob <- evinf_oob_fraction_summary(object)
    add("Out-of-bag fraction (block scheme)", "note",
        if (is.na(oob$oob_fraction_mean)) {
          "Not available (no usable bootstrap replicates)."
        } else {
          paste0("mean ", signif(oob$oob_fraction_mean, 3), " (range ",
                signif(oob$oob_fraction_min, 3), "-", signif(oob$oob_fraction_max, 3),
                ") of rows left out per replicate under bootstrap_scheme = \"",
                object$bootstrap_scheme, "\" -- overlapping blocks leave out less than i.i.d. resampling, so out-of-bag error from it is optimistic.")
        })
  }

  out <- dplyr::bind_rows(rows)
  class(out) <- c("evinf_check", class(out))
  out
}

#' @export
print.evinf_check <- function(x, ...) {
  cat("evinf diagnostic report\n")
  n_warn <- sum(x$status == "warning")
  n_note <- sum(x$status == "note")
  cat(sprintf("  %d check%s, %d warning%s, %d note%s\n\n",
             nrow(x), if (nrow(x) != 1) "s" else "",
             n_warn, if (n_warn != 1) "s" else "",
             n_note, if (n_note != 1) "s" else ""))
  icon <- c(ok = "[ok]", note = "[note]", warning = "[!!]")
  for (i in seq_len(nrow(x))) {
    cat(sprintf("%-8s %-50s\n", icon[[x$status[i]]], x$check[i]))
    cat("         ", x$detail[i], "\n", sep = "")
  }
  invisible(x)
}
