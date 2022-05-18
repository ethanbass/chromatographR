#' Preprocess time/wavelength data
#' 
#' Standard pre-processing of response matrices, consisting of a time axis and
#' a spectral axis (e.g. HPLC-DAD/UV data). For smooth data, like UV-VIS data,
#' the size of the matrix can be reduced by interpolation. By default,
#' the data are baseline-corrected in the time direction and smoothed in the 
#' spectral dimension.
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
#' Defaults to TRUE (unless you're on Windows).
#' @param interpolate_rows Logical. Whether to interpolate along dim1. Defaults
#' to TRUE.
#' @param interpolate_cols Logical. Whether to interpolate along dim2. Defaults
#' to TRUE.
#' @param mc.cores How many cores to use for parallel processing. Defaults to 2.
#' @param \dots Further optional arguments to
#' \code{\link[ptw:baseline.corr]{baseline.corr}}.
#' @return The function returns the preprocessed data matrix, with rownames and
#' colnames indicating the time points and wavelengths, respectively.
#' @author Ethan Bass
#' @note Adapted from
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
                          spec.smooth = TRUE,
                          maxI, parallel = TRUE, 
                          interpolate_rows = TRUE,
                          interpolate_cols = TRUE,
                          mc.cores=2, ...){
  if (parallel & .Platform$OS.type == "windows"){
    parallel <- FALSE
    warning("Parallel processing is not currently available on Windows.")
  }
  if (is.matrix(X)){
    X <- list(X)
    return_matrix <- TRUE
  } else return_matrix <- FALSE
  if (!is.list(X) | mean(sapply(X,is.matrix)) != 1)
    stop("X should be a matrix or a list of matrices")
  # if (missing(dim1) | missing(dim2)){
  #   warning("...Times or wavelengths not provided.
  #           Extrapolating matrix dimensions for interpolation.",
  #           immediate. = TRUE)
  # }
  if (missing(dim1) & interpolate_rows){
    warning("...Times not provided. Extrapolating from matrix dimensions for interpolation.",
            immediate. = TRUE)
    limits <- sapply(X,function(x){
      ts <- rownames(x)
      c(head(ts,1), tail(ts,1))})
    start <- round(max(as.numeric(limits[1,])),2)
    end <- floor(min(as.numeric(limits[2,])))
    dim1 <- seq(start,end,by=.01)
  }
  if (missing(dim2) & interpolate_cols){
    warning("...Wavelengths not provided. Extrapolating from matrix dimensions for interpolation.",
            immediate. = TRUE)
    dim2 <- as.numeric(colnames(X[[1]]))
  }
    if (parallel){
      if (length(find.package('parallel', quiet=TRUE))==0){
        stop("Parallel must be installed to enable parallel processing.")
    }
      X <- parallel::mclapply(X, FUN=preprocess_matrix,
               dim1=dim1,
               dim2=dim2,
               remove.time.baseline = remove.time.baseline,
               spec.smooth = spec.smooth, maxI = maxI,
               interpolate_rows = interpolate_rows,
               interpolate_cols = interpolate_cols,
               mc.cores=mc.cores,
               ...)
    } else{
      X <- lapply(X, FUN=preprocess_matrix,
             dim1=dim1,
             dim2=dim2,
             remove.time.baseline = remove.time.baseline,
             spec.smooth = spec.smooth, maxI = maxI,
             interpolate_rows = interpolate_rows,
             interpolate_cols = interpolate_cols,
             ...)
    }
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
                              maxI,
                              interpolate_rows = TRUE,
                              interpolate_cols = TRUE,
                              ...) {
  if (!is.matrix(X))
    stop("X should be a matrix!")
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
  
  if (spec.smooth)
    X <- t(apply(X, 1, function(xxx) smooth.spline(xxx)$y))
  
  if (remove.time.baseline)
    X <- apply(X, 2, baseline.corr, ...)
  if (min(X) < 0)
    # X <- X - min(X)
    X[X < 0] <- 0
  if (!missing(maxI))
    X <- maxI * X / max(X)
  
  dimnames(X) <- list(dim1, dim2)
  X
}
