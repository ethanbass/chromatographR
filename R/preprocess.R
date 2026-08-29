#' Preprocess time/wavelength data
#' 
#' Standard preprocessing of time × wavelength response matrices
#' (e.g. HPLC-DAD/UV data), consisting of: (i) baseline correction along the 
#' time axis, (ii) smoothing along the spectral axis, and (iii) optional 
#' interpolation in either dimension for dimensionality reduction. For densely 
#' sampled data, (e.g.,  UV-VIS spectra, the size of the matrix can be reduced 
#' by interpolation. By default, the data are baseline-corrected along the time
#' axis using [`baseline.corr`][ptw::baseline.corr], and smoothed along the
#' spectral axis using cubic smoothing splines ([`smooth.spline`]).
#' 
#' @importFrom ptw baseline.corr
#' @importFrom stats approx smooth.spline
#' @param X A numeric matrix (time × wavelength) or a list of such matrices.
#' Row names must correspond to time points and column names to wavelengths.
#' @param dim1 Numeric vector specifying the target time grid. The data will be
#' interpolated onto these time points. The range of the new values should not 
#' exceed the range of the original time points.
#' @param dim2 Numeric vector specifying the target wavelength grid.
#' The data will be interpolated onto these wavelengths. The range of the new 
#' values should not exceed the range of the original wavelengths.
#' @param remove.time.baseline Logical, indicating whether baseline correction
#' should be done along the time axis, according to
#' [`baseline.corr`][ptw::baseline.corr]. Defaults to `TRUE`.
#' @param spec.smooth Logical, indicating whether smoothing should be done along
#' the spectral axis, according to [`smooth.spline`]. Defaults to `TRUE`.
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
#' @param outlier_cutoff Threshold (in seconds) for excluding chromatograms
#' that end earlier than expected. Samples ending more than this value before
#' the median end time are removed. Defaults to `5` seconds. Only applies when 
#' `dim1` is not specified.
#' @param ... Additional arguments to [`baseline.corr`][ptw::baseline.corr].
#' @return The function returns the preprocessed data matrix (or list of 
#' matrices), with row names and column names indicating the time points and 
#' wavelengths, respectively.
#' @author Ethan Bass
#' @note Adapted from the `preprocess` function in the alsace package by 
#' Ron Wehrens: <https://github.com/rwehrens/alsace/blob/master/R/preprocess.R>.
#' @references 
#' * Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast parametric 
#' time warping of peak lists. *Bioinformatics* 31:3063-3065. 
#' \doi{10.1093/bioinformatics/btv299}.
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. 
#' *Metabolomics* 11:1:143-154. \doi{10.1007/s11306-014-0683-5}.
#' @examplesIf interactive()
#' data(Sa)
#' new.ts <- seq(10,18.66,by=.01) # choose time-points
#' new.lambdas <- seq(200, 318, by = 2) # choose wavelengths
#' Sa_pr <- preprocess(Sa[[1]], dim1 = new.ts, dim2 = new.lambdas)
#' @md
#' @export preprocess

preprocess <- function(X, dim1, dim2,
                       remove.time.baseline = TRUE,
                       spec.smooth = TRUE, maxI = NULL,
                       interpolate_rows = TRUE,
                       interpolate_cols = TRUE,
                       cl = 2, show_progress = NULL,
                       outlier_cutoff = 5/60, ...){
  laplee <- choose_apply_fnc(show_progress = show_progress, cl = cl)
  time_unit <- get_time_unit(X, na_value = "min")
  tfac <- switch(time_unit, "min" = 1, "s" = 60, "ms" = 60*1000)
  X_class <- class(X)
  if (is.matrix(X)){
    X <- list(X)
    return_matrix <- TRUE
  } else return_matrix <- FALSE
  if (base::inherits(X, "list") &&
      base::all(base::sapply(X, function(x) {
        base::inherits(x, "list") && "dad" %in% base::names(x) && base::is.matrix(x$dad)
      }))) {
    X <- base::lapply(X, function(x) x$dad)
  }
  if (!(inherits(X, "matrix") || (inherits(X, "list") && all(sapply(X, is.matrix)))))
    stop("X should be a matrix or a list of matrices")
  if (ncol(X[[1]]) == 1){
    interpolate_cols <- FALSE
  }
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
                        paste(sQuote(names(outliers)), collapse = ", ")),
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
        labeled <- paste0(too_early, " (", sQuote(names(too_early)), ")", collapse = ", ")
        plural <- ifelse(length(too_early) > 1, "s", "")
        warning(sprintf(paste(
          "No extrapolation allowed. New time axis (dim1) begins at %.2f, before actual data begins.",
          "\n\tExcluding chromatogram%s: %s."),
          head(dim1, 1), plural, labeled
        ), immediate. = TRUE)
      }
      
      if (length(too_late) > 0) {
        labeled <- paste0(too_late, " (", sQuote(names(too_late)), ")", collapse = ", ")
        plural <- ifelse(length(too_late) > 1, "s", "")
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
  X <- laplee(X, FUN = preprocess_matrix,
                          dim1 = dim1,
                          dim2 = dim2,
                          remove.time.baseline = remove.time.baseline,
                          spec.smooth = spec.smooth, maxI = maxI,
                          interpolate_rows = interpolate_rows,
                          interpolate_cols = interpolate_cols,
                          ...)
    if (return_matrix){
      X[[1]]
    } else{
      class(X) <- X_class
      X
    }
}

#' @noRd
preprocess_matrix <- function(X,
                              dim1, ## time axis
                              dim2, ## spectral axis
                              remove.time.baseline = TRUE,
                              spec.smooth = TRUE,
                              maxI = NULL,
                              interpolate_rows = TRUE,
                              interpolate_cols = TRUE,
                              ...) {
  if (!is.matrix(X))
    stop("X should be a matrix!")
  metadata <- attributes(X)
  metadata[c("dimnames", "names", "row.names", "dim", "class", "levels")] <- NULL
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
  if (ncol(X) == 1){
    spec.smooth <- FALSE
  }
  if (spec.smooth){
    X <- t(apply(X, 1, function(xxx) smooth.spline(xxx)$y))
  }
  if (remove.time.baseline){
    X <- apply(X, 2, ptw::baseline.corr, ...)
  }
  if (min(X, na.rm = TRUE) < 0){
    X[X < 0] <- 0
  }
  if (!is.null(maxI)){
    X <- maxI * X / max(X)
  }
  dimnames(X) <- list(dim1, dim2)
  attributes(X) <- c(attributes(X), metadata)
  X
}
