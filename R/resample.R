# Resolve evzinb()/evinb()'s bootstrap_scheme = argument: NULL defaults to
# "cluster" when block is given (today's behaviour, unchanged) else "iid";
# an explicit scheme is validated and, for "moving_block"/"stationary",
# checked against `time` being supplied (evinf_resample_ids() would error
# anyway, but here at evzinb()/evinb() call time names the right argument).
evinf_resolve_bootstrap_scheme <- function(scheme, block, time) {
  if (is.null(scheme)) {
    return(if (is.null(block)) "iid" else "cluster")
  }
  scheme <- match.arg(scheme, c("iid", "cluster", "moving_block", "stationary"))
  if (scheme == "cluster" && is.null(block)) {
    stop("bootstrap_scheme = \"cluster\" requires `block`.", call. = FALSE)
  }
  if (scheme %in% c("moving_block", "stationary") && is.null(time)) {
    stop("bootstrap_scheme = \"", scheme, "\" requires `time`.", call. = FALSE)
  }
  scheme
}

# round9 F (audit §5.7): panel / time-series bootstrap resampling schemes.
# evinf_resample_ids() is the single generator behind bootstrap_scheme =
# "iid" / "cluster" / "moving_block" / "stationary", replacing the two
# inline branches (plain sample() vs. block resampling) that used to be
# duplicated in bootrun_evzinb() and bootrun_evinb().

# Split 1:n by block_vec, preserving each unit's rows in their original
# (row-index) order -- never sorted by this function. Returns a named list
# of integer row-index vectors, one per unit, in first-appearance order.
evinf_split_by_unit <- function(n, block_vec) {
  if (is.null(block_vec)) {
    return(list(`1` = seq_len(n)))
  }
  split(seq_len(n), block_vec, drop = TRUE)
}

# Validate that time_vec is strictly increasing (no ties, no drops) within
# a unit's rows, in the order those rows already appear in the data -- a
# moving/stationary block bootstrap needs a well-defined "next" observation,
# and silently sorting would resample a different time order than the one
# the model was actually fit on (round9 F, review: "error, don't sort").
evinf_check_unit_time_order <- function(idx, time_vec, unit_label) {
  t <- time_vec[idx]
  if (anyNA(t)) {
    stop("evinf_resample_ids(): `time` has missing values in unit '", unit_label,
         "'.", call. = FALSE)
  }
  if (is.unsorted(t, strictly = TRUE)) {
    stop("evinf_resample_ids(): `time` is not strictly increasing within unit '",
         unit_label, "' (duplicated or out-of-order timestamps). Sort the data ",
         "by (block, time) before fitting -- this function never reorders it ",
         "for you.", call. = FALSE)
  }
}

# round10 0.4 (review §3): the check above only ran inside each bootstrap
# replicate, so a duplicated or out-of-order `time` under moving_block /
# stationary let the full-sample fit complete and only surfaced as every
# bootstrap replicate being a try-error. Called once, up front, from
# run_evzinb()/run_evinb()/add_bootstraps() -- before any EM fitting -- so it
# errors immediately instead. evinf_resample_ids()'s own per-replicate check
# is unchanged, still guarding any direct/standalone caller.
#
# When `block` is NULL, evinf_split_by_unit() treats the whole data set as
# one unit, which is the usual mistake this guards against: forgetting
# `block =` on panel data makes every unit's genuinely-increasing timestamps
# look like duplicates once concatenated, so the error names that as the
# likely cause.
evinf_validate_time <- function(n, scheme, block_vec, time_vec) {
  if (!scheme %in% c("moving_block", "stationary")) {
    return(invisible(NULL))
  }
  units <- evinf_split_by_unit(n, block_vec)
  for (u in names(units)) {
    tryCatch(
      evinf_check_unit_time_order(units[[u]], time_vec, u),
      error = function(e) {
        msg <- conditionMessage(e)
        if (is.null(block_vec) && grepl("not strictly increasing", msg, fixed = TRUE)) {
          msg <- paste0(
            msg, " This usually means `block =` is missing: without it, the ",
            "whole data set is treated as a single time series, so repeated ",
            "timestamps across what are really separate units look like ",
            "duplicates within it."
          )
        }
        stop(msg, call. = FALSE)
      }
    )
  }
  invisible(NULL)
}

# Default block length for one unit of length T (round9 F): ceiling(T^(1/3)),
# the usual block-bootstrap rule of thumb. Always <= T for T >= 1, so it
# never needs clamping the way a user-supplied block_length might.
evinf_default_block_length <- function(len) {
  max(1L, ceiling(len^(1 / 3)))
}

# Resolve block_length against block_vec's units into a named per-unit
# vector: NULL -> ceiling(T^(1/3)) per unit (messaged once here); a single
# number -> applied to every unit; an already-resolved named vector (what
# run_evzinb()/run_evinb() store on the fitted object once, per moving_block/
# stationary fit) -> used as-is, so a bootstrap's many replicates each call
# evinf_resample_ids() without ever re-hitting the NULL branch or its
# message.
evinf_resolve_block_length <- function(n, block_vec, block_length) {
  units <- evinf_split_by_unit(n, block_vec)
  lens <- vapply(units, length, integer(1))
  if (is.null(block_length)) {
    Ls <- vapply(lens, evinf_default_block_length, numeric(1))
    message("evinf_resample_ids(): block_length not supplied; using ceiling(T^(1/3)) ",
            "per unit (block length", if (length(unique(Ls)) > 1) "s used: " else " used: ",
            paste(sort(unique(Ls)), collapse = ", "), "). Pass `block_length` to override.")
  } else if (!is.null(names(block_length)) && setequal(names(block_length), names(units))) {
    Ls <- block_length[names(units)]
  } else {
    if (length(block_length) != 1L) {
      stop("evinf_resample_ids(): block_length must be NULL, a single number, or a ",
           "named vector matching the units in block_vec.", call. = FALSE)
    }
    Ls <- stats::setNames(rep(block_length, length(units)), names(units))
  }
  if (any(Ls > lens)) {
    bad <- names(units)[Ls > lens]
    stop("evinf_resample_ids(): block_length exceeds the length of unit",
         if (length(bad) > 1) "s " else " ", paste(bad, collapse = ", "),
         " (T = ", paste(lens[bad], collapse = ", "), "). Pass a smaller ",
         "block_length, or NULL to use ceiling(T^(1/3)) per unit.", call. = FALSE)
  }
  Ls
}

# Moving-block bootstrap indices for one unit already in time order (`idx`,
# length T): draw ceiling(T/L) overlapping length-L blocks (starts in
# 1:(T-L+1)) with replacement, concatenate, and trim to the first T -- the
# standard moving-block bootstrap (Kunsch 1989). Returns row indices (values
# from `idx`), not positions.
evinf_moving_block_unit <- function(idx, L) {
  T_len <- length(idx)
  L <- as.integer(round(L))
  if (L >= T_len) {
    return(idx)
  }
  n_blocks <- ceiling(T_len / L)
  starts <- sample.int(T_len - L + 1L, n_blocks, replace = TRUE)
  pos <- unlist(lapply(starts, function(s) s:(s + L - 1L)))
  idx[pos[seq_len(T_len)]]
}

# Stationary bootstrap indices for one unit already in time order (`idx`,
# length T): blocks of Geometric(1/L)-distributed length (mean L, support
# {1, 2, ...}) starting at a uniformly-drawn position, wrapping circularly
# within the unit once a block runs past the end (Politis & Romano 1994).
# Draws blocks until the concatenated length reaches T, then trims.
evinf_stationary_unit <- function(idx, L) {
  T_len <- length(idx)
  L <- as.integer(round(L))
  if (T_len == 1L) {
    return(idx)
  }
  pos <- integer(0)
  while (length(pos) < T_len) {
    start <- sample.int(T_len, 1L)
    len <- stats::rgeom(1L, prob = 1 / L) + 1L
    block <- ((start - 1L + seq_len(len) - 1L) %% T_len) + 1L
    pos <- c(pos, block)
  }
  idx[pos[seq_len(T_len)]]
}

#' Bootstrap resample indices under a chosen scheme
#'
#' The single generator behind \code{bootstrap_scheme = } in \code{evzinb()} /
#' \code{evinb()} / \code{add_bootstraps()} (round9 F, audit §5.7).
#'
#' @param n Number of rows in the estimation data.
#' @param scheme One of \code{"iid"} (plain row resampling), \code{"cluster"}
#'   (resample whole units named by \code{block_vec}, keeping every row of a
#'   drawn unit -- what \code{block = } has always done), \code{"moving_block"}
#'   or \code{"stationary"} (block-resample each unit's own time series, see
#'   Kunsch 1989 / Politis and Romano 1994).
#' @param block_vec Optional unit identifier, one value per row. Required for
#'   \code{"cluster"}. For \code{"moving_block"} / \code{"stationary"}, treated
#'   as a single unit (the whole data) when \code{NULL}.
#' @param time_vec Required for \code{"moving_block"} / \code{"stationary"}: a
#'   numeric/orderable time index, one value per row, that must already be
#'   strictly increasing within every unit's rows as they appear in the data
#'   -- this function never sorts it for you (see Details).
#' @param block_length Block length \eqn{L} for \code{"moving_block"} /
#'   \code{"stationary"}. \code{NULL} (the default) uses
#'   \eqn{\lceil T^{1/3} \rceil} for each unit's own length \eqn{T} (a message
#'   names the value(s) used); a user-supplied \code{block_length} that
#'   exceeds some unit's \eqn{T} is an error rather than a silent clamp.
#'
#' @details With overlapping blocks (\code{"moving_block"} /
#'   \code{"stationary"}) the set of rows never drawn in a given replicate is
#'   both smaller and more temporally correlated than under i.i.d. resampling,
#'   so an out-of-bag error estimated from it (see \code{\link{oob_evaluation}})
#'   is optimistic relative to genuine forecasting performance.
#'
#' @return An integer vector of length \code{n}: row indices into the
#'   estimation data, with replacement, resampled under \code{scheme}.
#'
#' @keywords internal
evinf_resample_ids <- function(n, scheme = c("iid", "cluster", "moving_block", "stationary"),
                               block_vec = NULL, time_vec = NULL, block_length = NULL) {
  scheme <- match.arg(scheme)

  if (scheme == "iid") {
    return(sample.int(n, n, replace = TRUE))
  }

  if (scheme == "cluster") {
    if (is.null(block_vec)) {
      stop("evinf_resample_ids(): scheme = \"cluster\" requires block_vec ",
           "(evzinb()/evinb()'s `block =`).", call. = FALSE)
    }
    uniques <- unique(block_vec)
    boot_units <- sample(uniques, length(uniques), replace = TRUE)
    return(unlist(lapply(boot_units, function(u) which(block_vec == u))))
  }

  # scheme %in% c("moving_block", "stationary")
  if (is.null(time_vec)) {
    stop("evinf_resample_ids(): scheme = \"", scheme, "\" requires `time`.",
         call. = FALSE)
  }
  units <- evinf_split_by_unit(n, block_vec)
  Ls <- evinf_resolve_block_length(n, block_vec, block_length)

  unit_fun <- if (scheme == "moving_block") evinf_moving_block_unit else evinf_stationary_unit
  out <- lapply(names(units), function(u) {
    idx <- units[[u]]
    evinf_check_unit_time_order(idx, time_vec, u)
    unit_fun(idx, Ls[[u]])
  })
  unlist(out)
}
