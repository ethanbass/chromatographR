#' Preprocess time/wavelength data
#' 
#' Standard pre-processing of response matrices, consisting of a time axis and
#' a spectral axis (e.g. HPLC-DAD/UV data). For smooth data, like UV-VIS data,
#' the size of the matrix can be reduced by interpolation. By default,
#' the data are baseline-corrected in the time direction
#' (\code{\link[ptw:baseline.corr]{baseline.corr}}) and smoothed in the 
#' spectral dimension using cubic smoothing splines
#' (\code{\link[stats:smooth.spline]{smooth.spline}}).
#' 
#' @import ptw
#' @importFrom stats approx smooth.spline
#' @param X A numerical data matrix, or list of data matrices. Missing values
#' are not allowed. If rownames or colnames attributes are used, they should be
#' numerical and signify time points and wavelengths, respectively.
#' @param dim1 A new, usually shorter, set of time points (numerical). The
#' range of these should not exceed the range of the original time points.
#' @param dim2 A new, usually shorter, set of wavelengths (numerical). The
#' range of these should not exceed the range of the original wavelengths.
#' @param remove.time.baseline Logical, indicating whether baseline correction
#' should be done in the time direction, according to
#' \code{\link[ptw:baseline.corr]{baseline.corr}}. Default is \code{TRUE}.
#' @param spec.smooth Logical, indicating whether smoothing should be done in
#' the spectral direction, according to
#' \code{\link[stats:smooth.spline]{smooth.spline}}. Default is \code{TRUE}.
#' @param maxI if given, the maximum intensity in the matrix is set to this
#' value.
#' @param interpolate_rows Logical. Whether to interpolate along the time axis 
#' (\code{dim1}). Defaults to \code{TRUE}.
#' @param interpolate_cols Logical. Whether to interpolate along the spectral 
#' axis (\code{dim2}). Defaults to \code{TRUE}.
#' @param cl Argument to \code{\link[pbapply]{pblapply}} or \code{\link[parallel]{mclapply}}.
#' Either an integer specifying the number of clusters to use for parallel
#' processing or a cluster object created by \code{\link[parallel]{makeCluster}}.
#' Defaults to 2. On Windows integer values will be ignored.
#' @param show_progress Logical. Whether to show progress bar. Defaults to 
#' \code{TRUE} if \code{\link[pbapply]{pbapply}} is installed.
#' @param \dots Further optional arguments to
#' \code{\link[ptw:baseline.corr]{baseline.corr}}.
#' @return The function returns the preprocessed data matrix (or list of 
#' matrices), with row names and column names indicating the time points and 
#' wavelengths, respectively.
#' @author Ethan Bass
#' @note Adapted from the
#' \href{https://github.com/rwehrens/alsace/blob/master/R/preprocess.R}{preprocess}
#' function in the alsace package by Ron Wehrens.
#' @references 
#' * Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast
#' parametric time warping of peak lists. \emph{Bioinformatics}
#' \bold{31}:3063-3065. \doi{10.1093/bioinformatics/btv299}.
#' 
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. \emph{
#' Metabolomics} \bold{11:1}:143-154. \doi{10.1007/s11306-014-0683-5}.
#' @examplesIf interactive()
#' data(Sa)
#' new.ts <- seq(10,18.66,by=.01) # choose time-points
#' new.lambdas <- seq(200, 318, by = 2) # choose wavelengths
#' Sa_pr <- preprocess(Sa[[1]], dim1 = new.ts, dim2 = new.lambdas)
#' @md
#' @export preprocess

preprocess <- function(X, dim1, ## time axis
                          dim2, ## spectral axis
                          remove.time.baseline = TRUE,
                          spec.smooth = TRUE, maxI = NULL,
                          interpolate_rows = TRUE,
                          interpolate_cols = TRUE,
                          cl = 2, show_progress = NULL, ...){
  laplee <- choose_apply_fnc(show_progress = show_progress, cl = cl)
  if (is.matrix(X)){
    X <- list(X)
    return_matrix <- TRUE
  } else return_matrix <- FALSE
  if (!is.list(X) | !all(sapply(X, is.matrix)))
    stop("X should be a matrix or a list of matrices")
  if (ncol(X[[1]]) == 1){
    interpolate_cols <- FALSE
  }
  if (interpolate_rows && missing(dim1)){
    message("...Times not provided. Extrapolating from matrix dimensions for interpolation.")
    limits <- sapply(X,function(x){
      ts <- rownames(x)
      c(head(ts,1), tail(ts,1))})
    start <- ceiling(max(as.numeric(limits[1,]))*100)/100
    end <- floor((min(as.numeric(limits[2,])))*100)/100
    dim1 <- seq(start, end, by = .01)
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
    } else X
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
  
  ## possibly resize matrix to a lower dimension - faster, noise averaging
  if (interpolate_rows){
    if (length(tpoints <- as.numeric(rownames(X))) == 0)
      tpoints <- seq_len(nrow(X))
    if (min(dim1) < min(tpoints) |
        max(dim1) > max(tpoints))
      stop("No extrapolation allowed - check dim1 argument")
    X <- apply(X, 2, function(xx) approx(tpoints, xx, dim1)$y)
  } else dim1 <- rownames(X)
  if (interpolate_cols){
    if (length(lambdas <- as.numeric(colnames(X))) == 0)
      lambdas <- seq_len(ncol(X))
    if (min(dim2) < min(lambdas) |
        max(dim2) > max(lambdas))
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
    X <- apply(X, 2, baseline.corr, ...)
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
