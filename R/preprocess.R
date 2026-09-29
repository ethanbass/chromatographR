#' Preprocess time/wavelength data
#' 
#' Standard preprocessing of time × wavelength response matrices
#' (e.g. HPLC-DAD/UV data), consisting of: (i) optional interpolation in either
#' dimension for dimensionality reduction, (ii) optional smoothing along the
#' spectral axis, and (iii) optional baseline correction along the time axis.
#' For densely sampled data (e.g., UV-VIS spectra), the size of the matrix can
#' be reduced by interpolation. By default, the data are smoothed along the
#' spectral axis using cubic smoothing splines ([`smooth.spline`]), then
#' baseline-corrected along the time axis by asymmetric least squares
#' ([`asysm`][ptw::asysm]). Three other baseline estimators are available
#' through `baseline`, and `baseline = "none"` skips the step; see
#' [`baseline_estimators`].
#' 
#' @details
#' Note that the operations are applied in the order given above: interpolation
#' first, then spectral smoothing, then baseline correction. Because the
#' baseline is estimated after the time axis is interpolated, `stiffness` is
#' tied to whichever `dim1` grid was chosen: the penalty applies on the index
#' axis, so the same value gives a different baseline on a differently sampled
#' grid. Choose `dim1` first, then tune `stiffness` against that grid.
#'
#' @importFrom ptw asysm
#' @importFrom stats approx smooth.spline median
#' @param X A numeric matrix (time × wavelength) or a list of such matrices.
#' Row names must correspond to time points and column names to wavelengths.
#' @param dim1 Numeric vector specifying the target time grid, usually built
#' with [`seq`]. The data will be interpolated onto these time points. The
#' range of the new values should not exceed the range of the original time
#' points, and the values must be strictly increasing. Note that writing
#' `c(from, to, by = x)` in place of `seq(from, to, by = x)` yields a
#' three-element vector rather than a sequence; this is detected and raises an
#' error.
#' @param dim2 Numeric vector specifying the target wavelength grid, usually
#' built with [`seq`]. The data will be interpolated onto these wavelengths.
#' The range of the new values should not exceed the range of the original
#' wavelengths, and the values must be strictly increasing (unless
#' `interpolate_cols = FALSE`, in which case `dim2` is ignored).
#' @param remove.time.baseline The `remove.time.baseline` argument is
#' deprecated. Please use `baseline = "none"` instead.
#' @param baseline Baseline estimator to use along the time axis. Either
#' asymmetric least squares (`"als"`, the default, via
#' [`asysm`][ptw::asysm]), asymmetrically reweighted penalized least squares
#' (`"arpls"`), adaptive smoothness parameter penalized least squares
#' (`"aspls"`), statistics-sensitive non-linear iterative peak-clipping
#' (`"snip"`), or `"none"` to skip baseline correction along the time axis
#' altogether. See [`baseline_estimators`] for a description of each.
#' @param attach_baseline Logical. Whether to retain the estimated baseline
#' alongside the corrected data as a `"baseline"` attribute. The smoothing
#' parameter that produced it is always recorded, as a `"baseline_param"`
#' attribute, whether or not the baseline itself is kept. The baseline is 
#' subtracted either way; this argument only controls whether the estimate is
#' kept for inspection. Defaults to `FALSE`: the baseline is estimated
#' separately for each wavelength, so it is the same size as the data and
#' retaining it roughly doubles memory use. Adding the attribute back to the returned data recovers the
#' uncorrected matrix only when `zero_floor = FALSE` and `maxI` is unset:
#' clipping alters the data without altering the stored estimate, and `maxI`
#' rescales both by a factor derived from the clipped data.
#' @param zero_floor Logical. Whether to clip negative values in the
#' corrected data to zero. Defaults to `FALSE`. Clipping rectifies the noise,
#' destroying the symmetric below-baseline residual that diagnostics of
#' over-fitting rely on, so it is off by default.
#' @param smooth_spectra Logical, indicating whether smoothing should be done along
#' the spectral axis, according to [`smooth.spline`]. Defaults to `TRUE`.
#' @param spec.smooth The `spec.smooth` argument is deprecated. Please use
#' `smooth_spectra` instead.
#' @param maxI If supplied, all values are rescaled so that the maximum
#' intensity equals `maxI`.
#' @param interpolate_rows Logical. Whether to interpolate along the time axis 
#' (`dim1`). Defaults to `TRUE`.
#' @param interpolate_cols Logical. Whether to interpolate along the spectral 
#' axis (`dim2`). Defaults to `TRUE`.
#' @param cl Either an integer specifying the number of cores to use for
#' parallel processing or a cluster object created by
#' [`makeCluster`][parallel::makeCluster]. Defaults to `2`. On Windows integer
#' values will be ignored.
#' @param show_progress Logical. Whether to show progress bar. Defaults to 
#' `TRUE` if [`pbapply`][pbapply::pbapply] is installed.
#' @param outlier_cutoff Threshold (in minutes) for excluding chromatograms
#' that end earlier than expected. Samples ending more than this value before
#' the median end time are removed. The value is converted into the time unit
#' of the data, so it is always expressed in minutes regardless of how the
#' chromatograms are stored. Defaults to `5/60`, i.e. 5 seconds. Only applies
#' when `dim1` is not specified.
#' @param ... Additional arguments to the selected `baseline` estimator: e.g.
#' `stiffness` and `p` for `"als"` ([`asysm`][ptw::asysm]), `stiffness` for
#' `"arpls"`, `stiffness` and `k` for `"aspls"`, or `m` for `"snip"`. Note that
#' `maxI` follows `...` in the signature and so must be named in
#' full: otherwise R's partial matching would bind `m` 
#' to `maxI` instead of passing them through to the estimator.
#' @return The function returns the preprocessed data matrix (or list of 
#' matrices), with row names and column names indicating the time points and 
#' wavelengths, respectively.
#' @author Ethan Bass
#' @family preprocessing functions
#' @seealso [`baseline_estimators`], [`correct_baseline`], [`crop_chroms`],
#' [`get_baseline`]
#' @note Adapted from the `preprocess` function in the alsace package by 
#' Ron Wehrens: <https://github.com/rwehrens/alsace/blob/master/R/preprocess.R>.
#' @inherit baseline_estimators references
#' @examples
#' data(Sa)
#' new.ts <- seq(10,18.66,by=.01) # choose time-points
#' new.lambdas <- seq(200, 318, by = 2) # choose wavelengths
#' Sa_pr <- preprocess(Sa[[1]], dim1 = new.ts, dim2 = new.lambdas)
#' @md
#' @export preprocess

preprocess <- function(X, dim1, dim2,
                       remove.time.baseline = NULL,
                       baseline = c("als", "arpls", "aspls", "snip", "none"),
                       attach_baseline = FALSE,
                       zero_floor = FALSE,
                       smooth_spectra = TRUE,
                       interpolate_rows = TRUE,
                       interpolate_cols = TRUE,
                       cl = 2, show_progress = NULL,
                       outlier_cutoff = 5/60, ...,
                       maxI = NULL,
                       spec.smooth = NULL){
  smooth_spectra <- resolve_deprecated(spec.smooth, smooth_spectra)
  laplee <- choose_apply_fnc(show_progress = show_progress, cl = cl)
  baseline_supplied <- !missing(baseline)
  baseline <- match.arg(baseline, c("als", "arpls", "aspls", "snip", "none"))
  baseline <- resolve_remove_time_baseline(remove.time.baseline, baseline,
                                           baseline_supplied)
  prep <- as_chrom_list(X)
  X <- prep$chrom_list
  time_unit <- get_time_unit(X, na_value = "min")
  tfac <- time_factor(time_unit)
  if (ncol(X[[1]]) == 1){
    interpolate_cols <- FALSE
  }
  # the order check applies only to an axis that is actually interpolated: a
  # supplied `dim2` is ignored when `interpolate_cols = FALSE`, and a scalar is
  # legitimate there
  if (!missing(dim1))
    check_interp_grid(dim1, "dim1", check_order = interpolate_rows)
  if (!missing(dim2))
    check_interp_grid(dim2, "dim2", check_order = interpolate_cols)
  if (interpolate_rows){
    limits <- sapply(X, function(x){
      ts <- as.numeric(rownames(x))
      c(head(ts, 1), tail(ts, 1))
    })
    if (missing(dim1)){
      message("...Times not provided. Extrapolating from matrix dimensions for interpolation.")
      start <- ceiling(max(limits[1,])*100)/100
      outliers <- which(limits[2,] < median(limits[2,]) - outlier_cutoff*tfac)
      if (length(outliers) > 0){
        limits <- limits[,-outliers]
        warning(sprintf("Excluding short chromatograms: %s.", 
                        format_list(names(outliers))),
                immediate. = TRUE)
        X <- X[-outliers]
      }
      end <- floor((min(limits[2,]))*100)/100
      dim1 <- seq(start, end, by = .01)
    } else{
      too_early <- which(limits[1,] > head(dim1, 1))
      too_late <- which(limits[2,] < tail(dim1, 1))
      outliers <- union(too_early, too_late)
      
      if (length(outliers) == length(X)) {
        stop(sprintf(
          "No extrapolation allowed. New time axis (dim1) range [%.2f, %.2f] is incompatible with actual data.",
          head(dim1, 1), tail(dim1, 1)
        ))
      }
      
      if (length(too_early) > 0) {
        labeled <- format_list(paste0(too_early, " (", sQuote(names(too_early)), ")"),
                               quote = FALSE)
        plural <- plural(length(too_early))
        warning(sprintf(paste(
          "No extrapolation allowed. New time axis (dim1) begins at %.2f, before actual data begins.",
          "\n\tExcluding chromatogram%s: %s."),
          head(dim1, 1), plural, labeled
        ), immediate. = TRUE)
      }
      
      if (length(too_late) > 0) {
        labeled <- format_list(paste0(too_late, " (", sQuote(names(too_late)), ")"),
                               quote = FALSE)
        plural <- plural(length(too_late))
        warning(sprintf(paste(
          "No extrapolation allowed. New time axis (dim1) ends at %.2f, after actual data ends.",
          "\n\tExcluding chromatogram%s: %s."),
          tail(dim1, 1), plural, labeled
        ), immediate. = TRUE)
      }
      
      # Remove all outliers
      if (length(outliers) > 0) {
        X <- X[-outliers]
      }
    }
  }
  if (interpolate_cols && missing(dim2)){
    message("...Wavelengths not provided. Extrapolating from matrix dimensions for interpolation.")
    dim2 <- as.numeric(colnames(X[[1]]))
  }
  if (attach_baseline && baseline == "none"){
    warning("`attach_baseline` is ignored when `baseline = \"none\"`, ",
            "since no baseline is estimated.", immediate. = TRUE)
    attach_baseline <- FALSE
  }
  warn_na_gaps(X, bridge = interpolate_rows || interpolate_cols,
               baseline = baseline)
  X <- laplee(X, FUN = preprocess_matrix,
                          dim1 = dim1,
                          dim2 = dim2,
                          dim1_names = if (interpolate_rows) format_axis(dim1),
                          dim2_names = if (interpolate_cols) format_axis(dim2),
                          dt = if (interpolate_rows) axis_dt(dim1),
                          baseline = baseline,
                          attach_baseline = attach_baseline,
                          zero_floor = zero_floor,
                          smooth_spectra = smooth_spectra, maxI = maxI,
                          interpolate_rows = interpolate_rows,
                          interpolate_cols = interpolate_cols,
                          ...)
  if (prep$return_matrix){
    X[[1]]
  } else{
    class(X) <- prep$class
    X
  }
}

#' @noRd
warn_na_gaps <- function(X, bridge, baseline){
  if (!bridge) return(invisible(NULL))
  gaps <- lapply(X, function(x){
    if (!anyNA(x)) return(NULL)
    runs <- apply(x, 2, function(y){
      r <- rle(is.na(y))
      if (any(r$values)) max(r$lengths[r$values]) else 0L
    })
    c(n = sum(is.na(x)), run = max(runs), dt = axis_dt(rownames(x)))
  })
  hit <- which(!vapply(gaps, is.null, logical(1)))
  if (length(hit) == 0) return(invisible(NULL))
  nm <- names(X)[hit]
  if (is.null(nm)) nm <- paste("sample", hit)
  detail <- vapply(seq_along(hit), function(i){
    g <- gaps[[hit[i]]]
    span <- if (is.finite(g[["dt"]])){
      sprintf(" = %.3g time units", g[["run"]] * g[["dt"]])
    } else ""
    sprintf("%s (%d missing, longest run %d timepoints%s)", nm[i],
            g[["n"]], g[["run"]], span)
  }, character(1))
  consequence <- if (baseline == "none") "" else
    ", so the baseline will be fit partly through interpolated data"
  warning(sprintf("Missing values will be bridged by linear interpolation%s.\n\t%s.",
                  consequence, format_list(detail, quote = FALSE)),
          immediate. = TRUE)
  invisible(NULL)
}

#' @noRd
subtract_baseline <- function(X, baseline,
                              dt = NA_real_, attach_baseline = TRUE,
                              ...){
  bl_args <- list(...)
  bl <- apply(X, 2, function(y)
    do.call(estimate_baseline, c(list(y, baseline = baseline, dt = dt), bl_args)))
  dim(bl) <- dim(X)
  # the baseline is the same size as the data, so it is not worth keeping alive
  # for the rest of the call when nobody will read it
  list(x = structure(X - bl, baseline = NULL),
       baseline = if (attach_baseline) bl else NULL,
       param = effective_baseline_param(baseline, bl_args, dt = dt))
}

#' @noRd
preprocess_matrix <- function(X,
                              dim1, ## time axis
                              dim2, ## spectral axis
                              baseline,
                              attach_baseline = FALSE,
                              zero_floor = FALSE,
                              smooth_spectra = TRUE,
                              interpolate_rows = TRUE,
                              interpolate_cols = TRUE,
                              ...,
                              maxI = NULL,
                              dim1_names = NULL,
                              dim2_names = NULL,
                              dt = NULL) {
  if (!is.matrix(X))
    stop("X should be a matrix!")
  metadata <- attributes(X)
  metadata[c("dimnames", "names", "row.names", "dim", "class", "levels",
             "baseline", "baseline_param")] <- NULL
  if (interpolate_rows){
    if (length(tpoints <- as.numeric(rownames(X))) == 0)
      tpoints <- seq_nrow(X)
    if (min(dim1) < min(tpoints) | max(dim1) > max(tpoints))
      stop("No extrapolation allowed - check dim1 argument")
    X <- apply(X, 2, function(xx) approx(tpoints, xx, dim1)$y)
  } else dim1 <- rownames(X)
  if (interpolate_cols){
    if (length(lambdas <- as.numeric(colnames(X))) == 0)
      lambdas <- seq_ncol(X)
    if (min(dim2) < min(lambdas) | max(dim2) > max(lambdas))
      stop("No extrapolation allowed - check dim2 argument")
    X <- t(apply(X, 1, function(xx) approx(lambdas, xx, dim2)$y)) 
  } else dim2 <- colnames(X)
  if (ncol(X) < 4){
    smooth_spectra <- FALSE
  }
  if (smooth_spectra){
    sm_x <- suppressWarnings(as.numeric(dim2))
    if (length(sm_x) != ncol(X) || !all(is.finite(sm_x)))
      sm_x <- seq_ncol(X)
    X <- t(apply(X, 1, function(xxx) smooth.spline(sm_x, xxx)$y))
  }
  if (baseline != "none"){
    # `ptw::baseline.corr(y, ...)` is `y - asysm(y, ...)` for a vector, so
    # subtracting the estimate explicitly is bit-identical to the previous
    # implementation while making the baseline available to keep.
    corrected <- subtract_baseline(X, baseline = baseline,
                                            dt = dt %||% axis_dt(dim1),
                                   attach_baseline = attach_baseline, ...)
    X <- corrected$x
    bl <- corrected$baseline
    bl_param <- corrected$param
  } else {
    bl_param <- NULL
  }
  if (zero_floor && min(X, na.rm = TRUE) < 0){
    X[X < 0] <- 0
  }
  if (!is.null(maxI)){
    scale_factor <- maxI / max(X, na.rm = TRUE)
    X <- scale_factor * X
    if (attach_baseline) bl <- scale_factor * bl
  }
  dimnames(X) <- list(
    dim1_names %||% (if (interpolate_rows) format_axis(dim1) else dim1),
    dim2_names %||% (if (interpolate_cols) format_axis(dim2) else dim2))
  attributes(X) <- c(attributes(X), metadata)
  if (attach_baseline){
    dimnames(bl) <- dimnames(X)
    attr(X, "baseline") <- bl
  }
  if (!is.null(bl_param)) attr(X, "baseline_param") <- bl_param
  X
}
