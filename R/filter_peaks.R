#' Filter peak lists
#'
#' Filters peaks from a `peak_list` object by height, area, width, retention
#' time, or fit diagnostics (`converged` and `window_limited`).
#' 
#' @param peak_list A `peak_list` object, consisting of a nested list of peak
#' tables, where the first level is the sample, and the second level is the
#' spectral wavelength. Every sample × wavelength component is described by a
#' `data.frame` where every row represents a peak, and the columns contain
#' information on retention time, peak shape, height, and area.
#' @param min_height Minimum peak height. Peaks whose `height` is `NA`, such as
#' failed fits kept by [`get_peaks`], are dropped if `min_height` is supplied
#' and retained otherwise.
#' @param min_area Minimum peak area. Peaks whose `area` is `NA` are treated
#' the same way.
#' @param min_sd The `min_sd` argument is deprecated. Please use `min_fwhm`
#' instead.
#' @param max_sd The `max_sd` argument is deprecated. Please use `max_fwhm`
#' instead.
#' @param min_rt Minimum retention time.
#' @param max_rt Maximum retention time.
#' @param min_fwhm Minimum peak width, as a full width at half maximum, in the
#' time units of the peak list. Peaks whose `FWHM` could not be measured
#' (`NA`) are retained.
#' @param max_fwhm Maximum peak width, as a full width at half maximum, in the
#' time units of the peak list. Expressed in the same units as `max_fwhm` in
#' [`get_peaks`], so the same value transfers between the two. Peaks whose
#' `FWHM` could not be measured (`NA`) are retained.
#' @param converged Logical. If `TRUE`, retains only peaks whose fit converged;
#' if `FALSE`, retains only those that did not. A peak has `converged == 0`
#' either because the optimizer stopped short of convergence or because no
#' optimizer ran at all, the latter occurring when the integration window held
#' fewer points than the model has parameters. Requires a `converged` column,
#' which peak lists made with `fit = "raw"` or before chromatographR 0.9.0 lack.
#' Peaks whose `converged` is `NA` are retained.
#' @param window_limited Logical. If `FALSE`, drops peaks whose fitted shape
#' hit a limit set by the integration window, meaning the window rather than
#' the data determined the reported shape; if `TRUE`, retains only those.
#' Requires a `window_limited` column, which peak lists made with
#' `fit = "raw"` or before chromatographR 0.9.0 lack. Peaks whose
#' `window_limited` is `NA` are retained.
#' @return A filtered `peak_list` object containing only peaks that satisfy the
#' specified criteria.
#' @author Ron Wehrens, Ethan Bass
#' @note The structure of `filter_peaks` is adapted from the `filterPeaks`
#' function in Ron Wehrens' [alsace](https://github.com/rwehrens/alsace)
#' package.
#' @seealso [`get_peaks`], [`filter_peaktable`]
#' @examples
#' data(Sa_warp)
#' pks <- get_peaks(Sa_warp[1], lambda = "210")
#' filter_peaks(pks, min_height = 100)
#' filter_peaks(pks, converged = TRUE, window_limited = FALSE)
#' @export filter_peaks

filter_peaks <- function(peak_list, min_height, min_area,
                         min_sd, max_sd, min_rt, max_rt,
                         min_fwhm, max_fwhm, converged, window_limited){
  if (missing(min_height) & missing(min_area) &
      missing(min_sd) & missing(max_sd) &
      missing(min_fwhm) & missing(max_fwhm) &
      missing(min_rt) & missing(max_rt) &
      missing(converged) & missing(window_limited)) {
    warning("Nothing to filter...")
    return(peak_list)
  }
  if (!missing(min_sd) || !missing(max_sd)){
    warning("The `min_sd` and `max_sd` arguments are deprecated. They are ",
            "expressed as a standard deviation, which is not a peak width for ",
            "any model except `\"gaussian\"`. Please use `min_fwhm` and ",
            "`max_fwhm` instead, which are expressed as a full width at half ",
            "maximum.", immediate. = TRUE)
    if (!missing(min_sd)){
      if (!missing(min_fwhm))
        stop("Supply either `min_fwhm` or `min_sd`, not both.")
      min_fwhm <- FWHM_CONST * min_sd
    }
    if (!missing(max_sd)){
      if (!missing(max_fwhm))
        stop("Supply either `max_fwhm` or `max_sd`, not both.")
      max_fwhm <- FWHM_CONST * max_sd
    }
  }
  x <- do.call(rbind, do.call(rbind, peak_list))
  keep_na_height <- missing(min_height)
  keep_na_area <- missing(min_area)
  if (!missing(min_height) && !anyNA(x$height)){
    if (min_height < min(x$height))
      warning("The 'min_height' is less than minimum peak height.")}
  if (!missing(min_area) && !anyNA(x$area)){
    if (min_area < min(x$area))
      warning("The 'min_area' is less than minimum peak area.")}
  measured <- x$FWHM[!is.na(x$FWHM)]
  if (!missing(min_fwhm) && length(measured) > 0){
    if (min_fwhm < min(measured))
      warning("The 'min_fwhm' is less than the minimum peak width.")}
  if (!missing(max_fwhm) && length(measured) > 0){
    if (max_fwhm > max(measured))
      warning("The 'max_fwhm' is greater than the maximum peak width.")}

  if (missing(max_fwhm)) max_fwhm <- Inf
  if (missing(min_fwhm)) min_fwhm <- 0
  if (missing(max_rt)) max_rt <- Inf
  if (missing(min_rt)) min_rt <- 0
  if (missing(min_height)) min_height <- 0
  if (missing(min_area)) min_area <- 0
  converged <- check_flag_filter(converged, "converged", x)
  window_limited <- check_flag_filter(window_limited, "window_limited", x)

  result <- lapply(peak_list,
         function(smpl)
           transfer_metadata(lapply(smpl,
                  function(comp)
                    comp[which((is.na(comp[, "FWHM"]) |
                              (comp[, "FWHM"] >= min_fwhm &
                                 comp[, "FWHM"] <= max_fwhm)) &
                           (comp[, "height"] >= min_height |
                              (keep_na_height & is.na(comp[, "height"]))) &
                           (comp[, "area"] >= min_area |
                              (keep_na_area & is.na(comp[, "area"]))) &
                           comp[, "rt"] >= as.numeric(min_rt) &
                           comp[, "rt"] <= as.numeric(max_rt) &
                           keep_by_flag(comp, "converged", converged) &
                           keep_by_flag(comp, "window_limited",
                                        window_limited)), , drop = FALSE]), smpl))
  transfer_metadata(result, peak_list, transfer_class = TRUE)
}

#' Validate a fit-diagnostic filter
#'
#' @param want The value supplied by the user, or missing.
#' @param col The column it filters on.
#' @param x The peak list, flattened, used to check the column is present.
#' @return `want` as a single logical, or `NULL` if it was not supplied.
#' @author Ethan Bass
#' @noRd
check_flag_filter <- function(want, col, x){
  if (missing(want) || is.null(want)) return(NULL)
  if (!is.logical(want) || length(want) != 1L || is.na(want)){
    stop(sprintf("`%s` must be either TRUE or FALSE.", col), call. = FALSE)
  }
  if (!col %in% colnames(x)){
    stop(sprintf(
      "This peak list has no `%s` column, so it cannot be filtered on. It is reported for the peak models but not for `fit = \"raw\"`, and was added in chromatographR 0.9.0.",
      col), call. = FALSE)
  }
  want
}

#' Peaks retained by a fit-diagnostic filter
#' @noRd
keep_by_flag <- function(comp, col, want){
  if (is.null(want)) return(rep(TRUE, nrow(comp)))
  is.na(comp[, col]) | comp[, col] == as.numeric(want)
}

#' Filter peak table
#' 
#' Filters features (columns) in a `peak_table`. Filtering is
#' applied consistently across all components of the peak table, including the 
#' intensity matrix (`tab`), feature metadata (`pk_meta`), and reference spectra
#' (if present). Filtering can be based on retention time, wavelength, feature
#' intensity summarized across samples (mean, median, or max), feature
#' sparsity (i.e., proportion of zero values across samples), and/or fit
#' convergence.
#' 
#' @param peak_table A peak_table object from [`get_peaktable`].
#' @param rts Vector of retention times used to select features. Values are
#' mapped to the closest retention times in the peak table. If no 
#' match is found within `tol`, the value is ignored and a warning is issued.
#' Cannot be combined with `min_rt` or `max_rt`.
#' @param min_rt Minimum retention time for features to be retained.
#' @param max_rt Maximum retention time for features to be retained.
#' @param min_value Minimum threshold for feature intensity summarized across
#' samples (using the method specified by `what`).
#' @param max_zeros Maximum allowed feature sparsity (proportion of zero values 
#' across samples).
#' @param min_converged Minimum value of the `converged` row of `pk_meta`,
#' between 0 and 1. That row summarizes `converged` over the peaks in each
#' feature by the `summarize_by` method given to [`get_peaktable`]: with
#' `"mean"`, it is the proportion of those peaks whose fit converged, so `0.8`
#' retains features that converged in at least 80% of them; with the default
#' `"weighted.mean"`, that proportion is weighted by peak intensity. Requires a
#' `converged` row in `pk_meta`, so it is unavailable for peak tables built from
#' `fit = "raw"` peak lists.
#' @param what Method used to summarize feature intensity across samples for
#' filtering. One of `"mean"`, `"median"` (default), or `"max"`.
#' @param lambda Wavelength(s) used to select features in the peak table.
#' Only features matching the specified wavelength(s) are retained.
#' @param tol Tolerance for matching of retention times to `rts`. A
#' feature is considered a match if the absolute difference is at most `tol`.
#' Defaults to `0`, which requires an exact match.
#' @return A filtered `peak_table` containing only features (columns) that
#' satisfy the specified criteria. The same feature indices are applied to all
#' `peak_table` components.
#' @author Ethan Bass
#' @seealso [`get_peaktable`], [`filter_peaks`]
#' @examples
#' data(pk_tab)
#' pk_tab <- filter_peaktable(pk_tab, min_rt = 10, max_rt = 16)
#' @export filter_peaktable

filter_peaktable <- function(peak_table, rts, min_rt, max_rt, min_value, 
                             max_zeros, lambda, min_converged,
                             what = c("median", "mean", "max"), tol = 0){
  check_peaktable(peak_table)
  if (!missing(rts) && (!missing(min_rt) || !missing(max_rt))){
    stop("Supply either `rts` or `min_rt`/`max_rt`, not both.", call. = FALSE)
  }
  if (missing(rts) & missing(min_rt) &
      missing(max_rt) & missing(min_value) & 
      missing(lambda) & missing(max_zeros) & missing(min_converged)) {
    warning("Nothing to filter...")
    return(peak_table)
  }
  peak_table$pk_meta["rt",] <- as.numeric(peak_table$pk_meta["rt",])
  what <- match.arg(what, c("median", "mean", "max"))
  if (!missing(rts)){
    rts <- as.numeric(rts)
    if (!inherits(rts, c("numeric"))){
      stop("`rts` should be a vector of retention times.")
    }
    idx.rt <- as.numeric(sapply(rts, function(rt){
      diffs <- abs(rt - peak_table$pk_meta["rt", ])
      ifelse(min(diffs) <= tol, which.min(diffs), NA)
    }))
    nas <- is.na(idx.rt)
    if (any(nas)){
      warning(paste0("The following retention times were not identified in the peak table: ",
                     paste(rts[nas], collapse = ', '),
              ". \n", "     You can try increasing the tolerance (`tol`) to permit fuzzier matching."),
              immediate. = TRUE)
      idx.rt <- idx.rt[-which(nas)]
    }
  } else if (!missing(min_rt) | !missing(max_rt)){
    if (missing(min_rt))
      min_rt <- 0
    if (missing(max_rt))
      max_rt <- Inf
    idx.rt <- which(peak_table$pk_meta["rt",] >= min_rt &
                      peak_table$pk_meta["rt",] <= max_rt)
  } else{idx.rt <- seq_along(peak_table$pk_meta)}
  if (!missing(min_value)){
    val <- apply(peak_table$tab, 2, eval(what), na.rm = TRUE)
    idx.val <- which(val >= min_value)
  } else (idx.val <- seq_along(peak_table$tab))
  if (!missing(max_zeros)){
    if (!is.numeric(max_zeros) || max_zeros < 0 || max_zeros > 1)
      stop("`max_zeros` must be a numeric value between 0 and 1.")
    prop_zero <- colMeans(peak_table$tab == 0, na.rm = TRUE)
    idx.zeros <- which(prop_zero <= max_zeros)
  } else (idx.zeros <- seq_along(peak_table$tab))
  if (!missing(lambda)){
    idx.lambda <- which(peak_table$pk_meta["lambda",] %in% lambda)
  } else (idx.lambda <- seq_along(peak_table$tab))
  if (!missing(min_converged)){
    if (!is.numeric(min_converged) || length(min_converged) != 1L ||
        is.na(min_converged) || min_converged < 0 || min_converged > 1){
      stop("`min_converged` must be a numeric value between 0 and 1.",
           call. = FALSE)
    }
    if (!"converged" %in% rownames(peak_table$pk_meta)){
      stop("This peak table has no `converged` row in `pk_meta`, so it cannot be filtered on. It is carried over from the peak list, which reports it for the peak models but not for `fit = \"raw\"`.",
           call. = FALSE)
    }
    idx.conv <- which(as.numeric(peak_table$pk_meta["converged", ]) >=
                        min_converged)
  } else (idx.conv <- seq_along(peak_table$tab))
  idx <- Reduce(intersect, list(idx.rt, idx.val, idx.lambda, idx.zeros,
                                idx.conv))
  peak_table[, idx]
}
