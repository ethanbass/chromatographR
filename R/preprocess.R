preprocess <- function(X,
                       dim1 = tpoints, ## time axis
                       dim2 = lambdas, ## spectral axis
                       remove.time.baseline = TRUE,
                       spec.smooth = TRUE,
                       maxI, ...) {
  if (!is.matrix(X))
      stop("X should be a matrix!")
  
  ## possibly resize matrix to a lower dimension - faster, noise averaging
  if (is.null(tpoints <- as.numeric(rownames(X)))) tpoints <- 1:nrow(X)
  if (is.null(lambdas <- as.numeric(colnames(X)))) lambdas <- 1:ncol(X)

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
      X <- X - min(X)
  if (!missing(maxI))
      X <- maxI * X / max(X)

  dimnames(X) <- list(dim1, dim2)
  X
}
