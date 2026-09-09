#' Log-likelihood profile over the candidate values of C_EV
#'
#' Returns the log-likelihood evaluated at each candidate value of \eqn{C_{EV}}
#' from the final EM iteration. Use it to judge whether \eqn{C_{EV}} is well
#' identified and whether the candidate range (\code{c.lim}) is wide enough.
#'
#' @param object A fitted \code{evzinb} / \code{evinb} object.
#' @return A tibble with columns \code{c} and \code{loglik}.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' c_profile(model)
#' }
c_profile <- function(object) {
  if (!inherits(object, c("evzinb", "evinb"))) {
    stop("`object` must be a fitted evzinb / evinb model.", call. = FALSE)
  }
  if (is.null(object$c_profile)) {
    stop("This model does not carry a C_EV profile (refit with evinf >= 0.10.0).",
         call. = FALSE)
  }
  tibble::as_tibble(object$c_profile)
}

#' Plot the log-likelihood profile over the candidate values of C_EV
#'
#' @param object A fitted \code{evzinb} / \code{evinb} object.
#' @return A \code{ggplot} object.
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#' plot_c_profile(model)
#' }
plot_c_profile <- function(object) {
  rlang::check_installed("ggplot2", "for plot_c_profile()")
  prof <- c_profile(object)
  c_ev <- object$coef$C
  on_boundary <- c_ev <= min(prof$c) || c_ev >= max(prof$c)

  p <- ggplot2::ggplot(prof, ggplot2::aes(x = .data$c, y = .data$loglik)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::geom_vline(xintercept = c_ev, linetype = "dashed") +
    ggplot2::labs(
      x = expression(C[EV]), y = "log-likelihood",
      title = "Log-likelihood profile over the candidate range for C_EV",
      subtitle = if (on_boundary) {
        "Warning: the estimate lies on the boundary of c.lim - widen it."
      } else {
        NULL
      }
    ) +
    ggplot2::theme_minimal()
  p
}
