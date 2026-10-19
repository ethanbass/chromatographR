#' Crop chromatograms
#'
#' Restricts chromatograms to a range of retention times and/or wavelengths by
#' subsetting the existing grid.
#'
#' @details Cropping is pure subsetting: the returned data are a submatrix of
#' the input, so values and dimnames pass through verbatim. This is the
#' difference from cropping via [`preprocess`], which reaches a narrower range
#' by interpolating onto a new axis — resampling the signal unless the requested
#' grid happens to coincide with the existing time points, and reformatting the
#' retention times either way.
#'
#' Bounds are inclusive, and either end of a range may be `NA` to leave it
#' open. Ranges are expressed in the units of the data itself: retention times
#' in whatever unit the chromatogram carries, as returned by [`get_times`], and
#' wavelengths as labeled in the column names, as returned by [`get_lambdas`].
#'
#' An `"estimated_baseline"` attribute attached by [`preprocess`] or [`correct_baseline`]
#' is cropped alongside the data, so it stays aligned and can still be
#' retrieved with [`get_baseline`].
#'
#' A range that does not overlap the data at all is an error. A range that
#' merely extends past one end is honored as far as the data allow, with a
#' warning reporting what was actually kept, since this is usually a typo or a
#' unit mismatch. Cropping chromatograms whose time axes are offset from one
#' another — as unaligned chromatograms generally are — can leave them with
#' different numbers of rows, which also warns: functions that read the time
#' axis from the first chromatogram would then mislabel the rest, so align
#' before cropping if the result needs a common time axis.
#'
#' @param x A numeric matrix (time × wavelength) or a list of such matrices. Row
#' names must correspond to time points and column names to wavelengths.
#' @param rt_range Numeric vector of length 2 giving the inclusive range of
#' retention times to keep, in the time unit of the data. Either element may be
#' `NA` for an open bound. Defaults to `NULL`, which keeps every time point.
#' @param lambda_range Numeric vector of length 2 giving the inclusive range of
#' wavelengths to keep. Either element may be `NA` for an open bound. Defaults
#' to `NULL`, which keeps every wavelength.
#' @return The cropped data, in the same shape as `x`.
#' @author Ethan Bass
#' @family preprocessing functions
#' @seealso [`preprocess`], [`correct_baseline`], [`get_times`], [`get_lambdas`]
#' @examples
#' data(Sa_pr)
#' cropped <- crop_chroms(Sa_pr, rt_range = c(12, 16))
#' trimmed <- crop_chroms(Sa_pr, lambda_range = c(220, NA))
#' @md
#' @export crop_chroms

crop_chroms <- function(x, rt_range = NULL, lambda_range = NULL){
  if (is.null(rt_range) && is.null(lambda_range))
    stop("Supply at least one of `rt_range` or `lambda_range`.")
  rt_range <- check_crop_range(rt_range, "rt_range")
  lambda_range <- check_crop_range(lambda_range, "lambda_range")
  prep <- as_chrom_list(x, arg = "x")
  # the axes are parsed once per sample and reused for both the crop and the
  # out-of-range warnings, and an axis with no range is not parsed at all
  axes <- lapply(prep$chrom_list, function(xx){
    list(rt = if (is.null(rt_range)) NULL
              else suppressWarnings(as.numeric(rownames(xx))),
         lambda = if (is.null(lambda_range)) NULL
                  else suppressWarnings(as.numeric(colnames(xx))))
  })
  out <- lapply(seq_along(prep$chrom_list), function(i)
    crop_matrix(prep$chrom_list[[i]], rt_range = rt_range,
                lambda_range = lambda_range, axis = axes[[i]]))
  names(out) <- names(prep$chrom_list)
  warn_crop_clamped(prep$chrom_list, out, rt_range, "rt_range", "rt", axes,
                    rownames, prep$return_matrix)
  warn_crop_clamped(prep$chrom_list, out, lambda_range, "lambda_range",
                    "lambda", axes, colnames, prep$return_matrix)
  if (prep$return_matrix){
    out[[1]]
  } else{
    warn_crop_ragged(out)
    class(out) <- prep$class
    out
  }
}

#' @noRd
sample_labels <- function(x){
  nm <- names(x)
  if (is.null(nm)) paste("sample", seq_along(x)) else sQuote(nm)
}

#' @noRd
warn_crop_clamped <- function(before, after, range, arg, which_axis, axes,
                              axis_fn, unlabeled = FALSE){
  if (is.null(range)) return(invisible(NULL))
  lo <- range[1]
  hi <- range[2]
  vals <- lapply(axes, `[[`, which_axis)
  hit <- which(vapply(vals, function(v){
    (is.finite(lo) && lo < min(v)) || (is.finite(hi) && hi > max(v))
  }, logical(1)))
  if (length(hit) == 0) return(invisible(NULL))
  nm <- sample_labels(before)
  detail <- vapply(hit, function(i){
    v <- suppressWarnings(as.numeric(axis_fn(after[[i]])))
    span <- sprintf("[%s, %s]", format(min(v)), format(max(v)))
    if (unlabeled) sprintf("kept %s", span)
    else sprintf("%s kept %s", nm[i], span)
  }, character(1))
  warning(sprintf(paste("`%s` extends beyond the available data, so the crop",
                        "was limited to what is present.\n\t%s."),
                  arg, format_list(detail, quote = FALSE)), immediate. = TRUE)
  invisible(NULL)
}

#' @noRd
warn_crop_ragged <- function(out){
  dims <- vapply(out, function(x) paste(dim(x), collapse = " x "), character(1))
  if (length(unique(dims)) < 2) return(invisible(NULL))
  warning(sprintf(paste("The cropped chromatograms have unequal dimensions:",
                        "%s.\n\tFunctions that read the time axis from the",
                        "first chromatogram will mislabel the others. Align the",
                        "chromatograms before cropping, or choose a range that",
                        "every sample covers identically."),
                  format_list(sprintf("%s (%s)", sample_labels(out), dims),
                              quote = FALSE)), immediate. = TRUE)
  invisible(NULL)
}

#' @noRd
check_crop_range <- function(range, arg){
  if (is.null(range)) return(NULL)
  if (length(range) == 2 && all(is.na(range)))
    stop(sprintf("`%s` must specify at least one bound.", arg))
  if (!is.numeric(range) || length(range) != 2)
    stop(sprintf("`%s` must be a numeric vector of length 2.", arg))
  if (!anyNA(range) && range[1] > range[2])
    stop(sprintf("`%s` must be increasing, but %s is greater than %s.", arg,
                 format(range[1]), format(range[2])))
  c(if (is.na(range[1])) -Inf else range[1],
    if (is.na(range[2])) Inf else range[2])
}

#' @noRd
crop_axis_idx <- function(vals, range, arg, what, n){
  if (is.null(range)) return(seq_len(n))
  if (length(vals) == 0 || anyNA(vals))
    stop(sprintf("Cannot crop by `%s` because the %s could not be read from ",
                 arg, what), "the dimnames.")
  keep <- which(vals >= range[1] & vals <= range[2])
  if (length(keep) == 0)
    stop(sprintf("`%s` does not overlap the data, which spans [%s, %s].", arg,
                 format(min(vals)), format(max(vals))))
  keep
}

#' @noRd
crop_matrix <- function(x, rt_range = NULL, lambda_range = NULL,
                        axis = list()){
  ri <- crop_axis_idx(axis$rt, rt_range, "rt_range", "retention times",
                      n = nrow(x))
  ci <- crop_axis_idx(axis$lambda, lambda_range, "lambda_range", "wavelengths",
                      n = ncol(x))
  bl <- attr(x, "estimated_baseline", exact = TRUE)
  out <- transfer_metadata(x[ri, ci, drop = FALSE], x, transfer_class = FALSE)
  if (!is.null(bl))
    attr(out, "estimated_baseline") <- bl[ri, ci, drop = FALSE]
  out
}
