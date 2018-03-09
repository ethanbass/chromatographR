opa <- function(x, ncomp, initXref = NULL)
{
  if (is.list(x))
    x <- do.call("rbind", x)

  lambdas <- colnames(x) ## may be NULL, no problem...
    
  x <- t(apply(x, 1, function(xx) xx / rep(sqrt(crossprod(xx)), length(xx))))
  
  Xref <- matrix(0, ncomp, ncol(x))

  if (is.null(initXref)) {
    huhn <- colMeans(x)
    Xref[1,] <- huhn / rep(sqrt(crossprod(huhn)), length(huhn))
    ncomp.done <- 0
  } else {
    if (ncol(initXref) != ncol(x)) 
      stop("Incompatible x and initXref arguments (different ncol)")
    
    ncomp.done <- nrow(initXref)
    Xref[1:ncomp.done,] <-             ## normalize to unit length
      t(apply(initXref, 1, function(x) x / rep(sqrt(crossprod(x)), length(x))))
  }

  for (i in (ncomp.done + 1):ncomp) {
    Xs <- lapply(1:nrow(x),
                 function(ii, xx, xref) rbind(xref, xx[ii,]),
                 x, Xref[1:(i-1),])
    dissims <- sapply(Xs, function(xx) det(tcrossprod(xx)))

    Xref[i,] <- x[which.max(dissims),]
  }
  colnames(Xref) <- lambdas

  t(Xref)
}
