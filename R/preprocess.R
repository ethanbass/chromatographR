#' Preprocessing smooth time-wavelength data
#' 
#' Standard preprocessing of response matrices where the first axis is a time
#' axis, and the second a spectral axis. An example is HPLC-DAD data. For
#' smooth data, like UV-VIS data, there is the option to decrease the size of
#' the matrix by interpolation. By default, the data are baseline-corrected in
#' the time direction and smoothed in the spectral dimension.
#' 
#' @import ptw
#' @importFrom stats approx smooth.spline
#' @param X A numerical data matrix, or list of data matrices. Missing values
#' are not allowed. If rownames or colnames attributes are used, they should be
#' numerical and signify time points and wavelengths, respectively.
#' @param dim1 A new, usually shorter, set of time points (numerical). The
#' range of these should not be outside the range of the original time points,
#' otherwise the function stops with an error message.
#' @param dim2 A new, usually shorter, set of wavelengths (numerical). The
#' range of these should not be outside the range of the original wavelengths,
#' otherwise the function stops with an error message.
#' @param remove.time.baseline Logical, indicating whether baseline correction
#' should be done in the time direction, according to
#' \code{\link[ptw:baseline.corr]{baseline.corr}}. Default is TRUE.
#' @param spec.smooth Logical, indicating whether smoothing should be done in
#' the spectral direction, according to
#' \code{\link[stats:smooth.spline]{smooth.spline}}. Default is TRUE.
#' @param maxI if given, the maximum intensity in the matrix is set to this
#' value.
#' @param parallel Logical, indicating whether to use parallel processing.
#' Defaults to TRUE.
#' @param mc.cores How many cores to use for parallel processing. Defaults to 2.
#' @param \dots Further optional arguments to
#' \code{\link[ptw:baseline.corr]{baseline.corr}}.
#' @return The function returns the preprocessed data matrix, with rownames and
#' colnames indicating the time points and wavelengths, respectively.
#' @author Ethan Bass
#' @note Adapted from the preprocess function authored by Ron Wehrens in the
#' alsace package.
#' @examples
#' \dontrun{
#' data(Sa)
#' new.ts <- seq(1,38,by=.01) # choose time-points
#' new.lambdas <- seq(200, 400, by = 2) # choose wavelengths
#' Sa.pr <-preprocess(Sa[[1]], dim1 = new.ts, dim2 = new.lambdas)
#' }
#' @export preprocess

preprocess <- function(X,
                              dim1, ## time axis
                              dim2, ## spectral axis
                              remove.time.baseline = TRUE,
                              spec.smooth = TRUE,
                              maxI, parallel=TRUE, mc.cores=2, ...){
  if (!is.list(X) & !is.matrix(X)){
    stop("X should be a matrix or a list of matrices")
  }
  if (is.list(X)){
    if (mean(sapply(X,is.matrix)) != 1){
      stop("X should be a matrix or a list of matrices")
    }
    if (parallel){
      if (length(find.package('parallel', quiet=TRUE))==0){
        stop("Parallel must be installed to enable parallel processing.")
      }
      parallel::mclapply(X, FUN=preprocess_matrix,
               dim1=dim1,
               dim2=dim2,
               remove.time.baseline = remove.time.baseline,
               spec.smooth = spec.smooth,
               maxI=maxI, mc.cores=mc.cores,
               ...)
    } else{
      lapply(X, FUN=preprocess_matrix,
             dim1=dim1,
             dim2=dim2,
             remove.time.baseline = remove.time.baseline,
             spec.smooth = spec.smooth,
             maxI=maxI,
             ...)
    }
    
  } else{
    preprocess_matrix(X, dim1=dim1,
                      dim2=dim2,
                      remove.time.baseline = remove.time.baseline,
                      spec.smooth = spec.smooth,
                      maxI=maxI,
                      ...)
  }
}

#' @noRd
preprocess_matrix <- function(X,
                              dim1, ## time axis
                              dim2, ## spectral axis
                              remove.time.baseline = TRUE,
                              spec.smooth = TRUE,
                              maxI, ...) {
  if (!is.matrix(X))
    stop("X should be a matrix!")
  
  ## possibly resize matrix to a lower dimension - faster, noise averaging
  if (length(tpoints <- as.numeric(rownames(X))) == 0)
    tpoints <- seq_len(nrow(X))
  if (length(lambdas <- as.numeric(colnames(X))) == 0)
    lambdas <- seq_len(ncol(X))
  if (min(dim1) < min(tpoints) |
      max(dim1) > max(tpoints))
    stop("No extrapolation allowed - check dim1 argument")
  
  X <- apply(X, 2, function(xx) approx(tpoints, xx, dim1)$y)
  
  if (min(dim2) < min(lambdas) |
      max(dim2) > max(lambdas))
    stop("No extrapolation allowed - check dim2 argument")
  
  X <- t(apply(X, 1, function(xx) approx(lambdas, xx, dim2)$y))
  
  if (spec.smooth)
    X <- t(apply(X, 1, function(xxx) smooth.spline(xxx)$y))
  
  if (remove.time.baseline)
    X <- apply(X, 2, baseline.corr, ...)
  if (min(X) < 0)
    # X <- X - min(X)
    X[X<0] <- 0
  if (!missing(maxI))
    X <- maxI * X / max(X)
  
  dimnames(X) <- list(dim1, dim2)
  X
}
