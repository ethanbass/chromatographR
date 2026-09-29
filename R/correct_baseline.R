#' Correct baseline along the time axis
#'
#' Estimates and subtracts a baseline from each chromatogram in a `chrom_list`,
#' along the time axis, using any of the available baseline estimators (see
#' [`baseline_estimators`]).
#'
#' @details Unlike [`preprocess`], this function only subtracts the baseline; it
#' does not interpolate, smooth, or rescale.
#'
#' @param x A numeric matrix (time × wavelength), a list of such matrices, or a
#' plain numeric vector for a single trace. Row names (or, for a vector,
#' names) must correspond to time points and column names to wavelengths.
#' @param baseline Baseline estimator to use along the time axis. Either
#' asymmetric least squares (`"als"`, the default, via
#' [`asysm`][ptw::asysm]), asymmetrically reweighted penalized least squares
#' (`"arpls"`), adaptive smoothness parameter penalized least squares
#' (`"aspls"`), or statistics-sensitive non-linear iterative peak-clipping
#' (`"snip"`). See [`baseline_estimators`] for a description of each.
#' @param attach_baseline Logical. Whether to retain the estimated baseline
#' alongside the corrected data, as a `"baseline"` attribute on each returned
#' matrix, retrievable with [`get_baseline`]. The baseline is subtracted either
#' way; this argument only controls whether the estimate is kept for
#' inspection. Defaults to `FALSE`: the baseline is estimated separately for
#' each wavelength, so it is the same size as the data and retaining it
#' doubles memory use.
#' @param zero_floor Logical. Whether to clip negative values in the corrected
#' data to zero. Defaults to `FALSE`. Clipping rectifies the noise, destroying
#' the symmetric below-baseline residual that diagnostics of over-fitting rely
#' on, so it is off by default. It can be useful with `baseline = "als"`, which
#' can overshoot the signal; `"snip"` cannot overshoot by construction, so
#' clipping is never needed there. Note that [`subtract_blanks`] also has a
#' `zero_floor` argument, but it defaults to `TRUE`, since a negative value
#' there means the blank exceeded the sample rather than that an estimator
#' overshot.
#' @param cl Either an integer specifying the number of cores to use for
#' parallel processing or a cluster object created by
#' [`makeCluster`][parallel::makeCluster]. Defaults to `2`. On Windows integer
#' values are ignored.
#' @param show_progress Logical. Whether to show progress bar. Defaults to
#' `TRUE` if [`pbapply`][pbapply::pbapply] is installed.
#' @param ... Additional arguments to the selected `baseline` estimator: e.g.
#' `stiffness` and `p` for `"als"` ([`asysm`][ptw::asysm]), `stiffness` for
#' `"arpls"`, `stiffness` and `k` for `"aspls"`, or `m` for `"snip"`.
#' @return The baseline-corrected data, in the same shape as `x`, with row and
#' column names (or, for a vector, names) unchanged.
#' @author Ethan Bass
#' @family preprocessing functions
#' @seealso [`baseline_estimators`], [`preprocess`], [`get_baseline`],
#' [`subtract_blanks`]
#' @inherit baseline_estimators references
#' @examples
#' data(Sa_pr)
#' bc <- correct_baseline(Sa_pr[1], cl = 1, show_progress = FALSE)
#' @md
#' @export correct_baseline

correct_baseline <- function(x,
                             baseline = c("als", "arpls", "aspls", "snip"),
                             attach_baseline = FALSE,
                             zero_floor = FALSE,
                             cl = 2, show_progress = NULL, ...){
  baseline <- match.arg(baseline, c("als", "arpls", "aspls", "snip"))
  vector_in <- is.numeric(x) && is.null(dim(x))
  if (vector_in) x <- as.matrix(x)
  laplee <- choose_apply_fnc(show_progress = show_progress, cl = cl)
  prep <- as_chrom_list(x, arg = "x")
  out <- laplee(prep$chrom_list, FUN = correct_baseline_matrix,
                baseline = baseline,
                attach_baseline = attach_baseline,
                zero_floor = zero_floor,
                ...)
  if (prep$return_matrix){
    if (vector_in) undim_trace(out[[1]]) else out[[1]]
  } else{
    class(out) <- prep$class
    out
  }
}

#' @noRd
correct_baseline_matrix <- function(x, baseline,
                                    attach_baseline = FALSE,
                                    zero_floor = FALSE, ...){
  if (!is.matrix(x))
    stop("x should be a matrix!")
  corrected <- subtract_baseline(x, baseline = baseline,
                                        dt = axis_dt(rownames(x)),
                                 attach_baseline = attach_baseline, ...)
  out <- corrected$x
  if (zero_floor && min(out, na.rm = TRUE) < 0){
    out[out < 0] <- 0
  }
  if (attach_baseline){
    bl <- corrected$baseline
    dimnames(bl) <- dimnames(out)
    attr(out, "baseline") <- bl
  }
  if (!is.null(corrected$param)) attr(out, "baseline_param") <- corrected$param
  out
}

#' Convert a single-column corrected matrix back to a named vector
#' @noRd
undim_trace <- function(x){
  nm <- rownames(x)
  bl <- attr(x, "baseline", exact = TRUE)
  dim(x) <- NULL
  names(x) <- nm
  if (!is.null(bl)) attr(x, "baseline") <- setNames(as.vector(bl), nm)
  x
}
