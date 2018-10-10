## return those components that never reach a certain
## concentraction value - these might be taken out
smallComps <- function(obj, Ithresh) {
  maxCvalues <- apply(sapply(obj$CList, function(x) apply(x, 2, max)),
                      1, max)
  list(smallComps = which(maxCvalues < Ithresh),
       maxCvalues = maxCvalues)
}

removeComps <- function(obj, toRemove, ...) {
  if (!all(toRemove %in% 1:ncol(obj$S)))
      stop("Invalid component vector: should be a subset of", 1:ncol(obj$S))

  ## regenerate the data, should be quick
  PsiList <- lapply(1:length(obj$CList),
                    function(i)
                    tcrossprod(obj$CList[[i]], obj$S) + obj$resid[[i]])
  names(PsiList) <- names(obj$CList)

  doALS(PsiList, obj$S[, -toRemove], ...)
}

## function combineComps allows the user to combine components.
## Suppose that in a set of 8 components the _real_
## component number 7 is a combination of comps 2 and 7. compList is
## then given like this: compList <- list(1, 2, 3, 4, 5, 6, c(2, 7),
## 8). Weights are given by the max intensities in the CList elements,
## or explicitly provided by the user.

combineComps <- function(obj, compList, weights, ...) {
  lambdas <- getWavelength(obj)
  ## regenerate the data, should be quick
  PsiList <- lapply(1:length(obj$CList),
                    function(i)
                    tcrossprod(obj$CList[[i]], obj$S) + obj$resid[[i]])
  names(PsiList) <- names(obj$CList)
  
  CLlength <- sapply(compList, length)
  if (max(CLlength) == 1)
      stop("Nothing to do: no combination of components given")

  if (missing(weights)) {
    maxCvalues <- apply(sapply(obj$CList,
                               function(x) apply(x, 2, max)),
                        1, max)
    scaledS <- sweep(obj$S, 2, maxCvalues, FUN = "*")
    newS <- sapply(compList,
                   function(x)
                     rowMeans(scaledS[, x, drop = FALSE]))
  } else {
    combineS <- function(origS, idx, wght) {
      rowSums(sweep(origS[,idx, drop = FALSE], 2, wght, FUN = "*"))
    }
    ## weights should be as long as compList, individual elements
    ## should have the same length, too
    newS <- sapply(seq(along = compList),
                   function(ii)
                     combineS(obj$S, compList[[ii]], weights[[ii]]))
  }
  
  ## scale to unit length
  newS <- apply(newS, 2, function(xx) xx/rep(sqrt(crossprod(xx)), length(xx)))
  dimnames(newS) <- list(lambdas, paste("Component", 1:ncol(newS)))

  doALS(PsiList, newS, ...)
}

## Experimental function, based on analysis of TEA data.
suggestCompCombis <- function(obj, indices, Ithresh = 0,
                              corthresh = 0.9,
                              clusterHeight = 0.6) {
  ## auxiliary function works on a sublist of the input argument
  cC.aux <- function(XCL) {
    mat1 <- matrix(0, ncol(XCL[[1]]), ncol(XCL[[1]]))
    for (s in 1:length(XCL)) {
      ## consider only those components with a max I above the threshold
      maxI <- which(apply(XCL[[s]], 2, max) > Ithresh)
      if (length(maxI) > 0) {
        ## find all components with a minimum correlation level
        corX <- which(cor(XCL[[s]]) > corthresh, arr.ind = TRUE)
        corX <- corX[apply(corX, 1, diff) < 0,,
                     drop = FALSE] ## no double entries
        corX <- corX[apply(corX, 1,
                           function(x) any(match(x, maxI, nomatch = 0))),,
                     drop = FALSE]
        
        ## if successful on both criteria, add one to the count
        mat1[corX] <- mat1[corX] + 1
      }      
    }
    
    mat1 / length(XCL)
  }

  huhn <- lapply(indices, function(idx) cC.aux(obj$CList[idx]))
  huhn.dist <- as.dist(1 - do.call("+", huhn))
  huhn.hcl <- hclust(huhn.dist, method = "single")
  huhn.clusters <- cutree(huhn.hcl, h = clusterHeight)

  ## reformat the cluster vector into something that can be directly
  ## imported into combineComps...
  list(clusters = lapply(1:max(huhn.clusters),
           function(x) which(huhn.clusters == x)),
       agreements = huhn)
}
 
