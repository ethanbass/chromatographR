## Jan 14, 2014: call doALS for one iteration after having identified which
## components are the same across time windows:
mergeTimeWindows <- function(obj, simSThreshold = .9, simCThreshold = .9,
                             verbose = FALSE)
{
  if (class(obj[[1]]) == "ALS") {
    overlap <- sum(rownames(obj[[1]]$CList[[1]]) %in%
                   rownames(obj[[2]]$CList[[1]])) / 2
    
    ## merge ALS objects: all components that are representing the
    ## same spectrum in different time windows need to be merged into
    ## one. Criteria: spectra similarity, and if overlap is present,
    ## also overlap in chromatographic profiles. In the latter case,
    ## if the similarity is big enough, the chromatographic profiles
    ## will be a weighted average. This should lead to smooth
    ## transitions across the window boarders.
    nwindows <- length(obj)
    window.times <- sapply(obj, function(x) rownames(x$CList[[1]]))
    timepoints <- sort(as.numeric(unique(unlist(window.times)))) ## overall
    nsamples <- length(obj[[1]]$CList)
    lambdas <- getWavelength(obj[[1]])
    nlambdas <- length(lambdas)
    ncomps <- sapply(obj, function(x) ncol(x$CList[[1]]))
    
    ## additional check on spectral similarity threshold: spectra from
    ## the same window should _not_ be put in the same cluster. So we
    ## use the max of the within-window correlations as a lower bound.
    wwcors <- sapply(obj,
                     function(x) {
                       huhn <- cor(x$S)
                       max(huhn[row(huhn)>col(huhn)])
                     })
    if (simSThreshold < max(wwcors)) {
      simSThreshold <- max(wwcors)+.001
      warning("Increased simSThreshold to avoid merging components within one time window")
    }
    
    ## first check spectral similarity. Step 1: gather all spectra.
    allS <- do.call(cbind, lapply(obj, function(x) x$S))
    Swin <- rep(1:length(obj), sapply(obj, function(x) ncol(x$S)))
    Snr <- unlist(sapply(obj, function(x) 1:ncol(x$S)))
    
    ## Step 2: crude clustering. First create a binary distance matrix
    ## based on the similarity threshold for spectral correlations,
    ## and then build a complete linkage tree that is cut at height
    ## 0.6. 
    corS <- cor(allS)
    identS <- 1 - (corS > simSThreshold)
    iS.cl <- cutree(hclust(as.dist(identS), method = "complete"), h = .6)
    
    ## if overlap is present check for the clusters in iS.cl the
    ## similarity in elution profiles. Things stay in the same cluster
    ## unless there are real differences.
    if (overlap > 0) {
      doubles <- as.numeric(names(which(table(iS.cl) > 1)))
      for (i in seq(along = doubles)) {
        dbl.idx <- which(iS.cl == doubles[i])
        iC.cl <- matchCprofiles(obj, Swin[dbl.idx], Snr[dbl.idx],
                                overlap, simCThreshold)
        
        if (length(unique(iC.cl)) > 1)
            iS.cl[dbl.idx] <-   ## simply adjust the label!
                paste(doubles[i], letters[iC.cl], sep = "")
      }
    }
    
    ## rename all clusters so that they again correspond to numbers
    iS.cl <- factor(as.integer(factor(iS.cl)))
    
    ## we have to count doubles again!
    doubles <- names(which(table(iS.cl) > 1))
    nS <- nlevels(iS.cl)
    S <- t(aggregate(t(allS), list(iS.cl), mean)[,-1])
    dimnames(S) <- list(lambdas, paste("Component", 1:nS))
    
    if (verbose){
      if (max(table(iS.cl)) > 1) {
        cat("\nMerging the following components:")
        iS.cl <- as.numeric(iS.cl)
        clusters <- lapply(1:max(iS.cl),
                           function(x) which(iS.cl == x))
        merges <- which(sapply(clusters, length) > 1)
        for (ii in merges)
            cat("\n\tComponents",
                paste(Snr[ clusters[[ii]] ], " (window ",
                      Swin[ clusters[[ii]] ], ")", sep = "", collapse = ", "))
        cat("\n")
      } ## else: nothing is merged!
    }
    
    ## simply do one iteration to estimate concentration profiles over
    ## the whole time range with the new spectra
    ## First re-generate the data for all windows
    PsiListList <- lapply(obj,
                          function(wobj)
                          lapply(1:length(wobj$CList),
                                 function(i)
                                 tcrossprod(wobj$CList[[i]], wobj$S) +
                                 wobj$resid[[i]]))
    PsiList <- merge.mats(PsiListList, overlap)
    names(PsiList) <- names(obj[[1]]$CList)
    
    newObj <- doALS(PsiList, PureS = S, maxiter = 1)
    attr(newObj, "overlap") <- overlap
    attr(newObj, "splitpoints") <-
        as.numeric(sapply(obj[-1],
                          function(x)
                          rownames(x$CList[[1]])[overlap]))
    
    newObj
  } else {   ## merge data lists
    if (is.null(overlap <- attr(obj, "overlap")))
        overlap <-
            sum(rownames(obj[[1]]) %in% rownames(obj[[1]])) / 2
        
    merge.mats(obj, overlap)
  }
}


## Addition Oct 21, 2013: overlapping time windows
splitTimeWindow <- function(datalist, splitpoints, overlap = 0) {
  tpoints <- as.numeric(rownames(datalist[[1]]))
  if (length(tpoints) == 0) {## no rownames available!
    tpoints <- 1:nrow(datalist[[1]])
    for (i in 1:length(datalist))
        rownames(datalist[[i]]) <- tpoints
  }
  if (min(splitpoints) < min(tpoints) | max(splitpoints) > max(tpoints))
      stop("Splitpoints outside time range")
  
  split.idx <- sapply(splitpoints, function(x) which(tpoints > x)[1])
  start.idx <- c(1, split.idx - overlap)
  end.idx <- c(split.idx + overlap - 1, length(tpoints))
  ## if (any(diff(start.idx) < 0 | diff(end.idx) < 0))
  ##     stop("Overlap larger than window size")
  if (any(diff(start.idx) < overlap | diff(end.idx) < overlap))
      stop("Non-overlapping area too small: should be at least as large as the overlap")
  
  splitData <- lapply(1:length(start.idx),
                      function(ii) {
                        lapply(datalist,
                               function(xx) xx[start.idx[ii]:end.idx[ii],])
                      })
  names(splitData) <- paste("Window", 1:(length(start.idx)))
  
  attr(splitData, "overlap") <- overlap
  
  splitData
}



## auxiliary function to merge a list of matrices.
## Time is in the rows, wavelength in the columns...
## we eliminate the overlapping bits at the beginning of all but
## the first window, and update the last rows of the
## preceding window accordingly. NOT EXPORTED.
unify.mats <- function(matlist, overlap, wghts) {
  if (overlap == 0) { ## for completeness
    do.call("rbind", matlist)
  } else {
    nrows <- sapply(matlist, nrow)
    for (i in 2:length(matlist)) {
      idx <- (nrows[i-1] - 2*overlap + 1):nrows[i-1]
      matlist[[i-1]][idx, ] <-
          sweep(matlist[[i-1]][idx, ,drop = FALSE], 1, rev(wghts),
                FUN = "*") +
                    sweep(matlist[[i]][1:(2*overlap), ,drop = FALSE], 1, wghts,
                          FUN = "*")
    }
    c(matlist[1],
      lapply(matlist[2:length(matlist)],
             function(x) x[-(1:(2*overlap)),,drop = FALSE]))
  }
}
  
  ## second auxiliary function to get from a list with time windows as
  ## the top level, and sample-based matrices as the second level, to
## one list that is sample-based. NOT EXPORTED.
  merge.mats <- function(windowlist, overlap) {
    if (overlap > 0) wghts <- (1:(2*overlap)) / (2*overlap + 1)

    mats <- lapply(1:length(windowlist[[1]]),
                   function(y)
                 lapply(windowlist, function(samp) samp[[y]]))
  names(mats) <- names(windowlist[[1]])
  
  if (overlap > 0)
      mats <- lapply(mats, unify.mats, overlap, wghts)
  
  lapply(mats, function(x) do.call("rbind", x))
}
  

## function matchCprofiles calculates correlation coefficients
## between the relevant parts of the elution profiles. If a
## coefficient is higher than the threshold, a match is recorded
## in a similarity matrix, which later again is used as input
## for a complete-linkage clustering. NOT EXPORTED.
matchCprofiles <- function(alsObj, winNrs, compNrs, overlap, simCThreshold) {
  ## we might try to calculate a correlation for cases where
  ## there is no signal, i.e. a sd of zero, leading to a missing
  ## value. This will be replaced by zero, but the user should
  ## not be alarmed by a warning...
  owarn <- options(warn = -1)
  on.exit(options(owarn))
  
  result <- matrix(0, length(winNrs), length(winNrs)) ## all equal
  for (i in seq(along = winNrs[-1])+1) {
    if (diff(winNrs[(i-1):i]) == 1) { ## neighbours
      Cprofiles1 <- lapply(alsObj[[winNrs[i-1]]]$CList,
                           function(x) x[,compNrs[i-1] ])
      Cprofiles2 <- lapply(alsObj[[winNrs[i]]]$CList,
                           function(x) x[,compNrs[i] ])
      n1 <- length(Cprofiles1[[1]])
      corC <- sapply(seq(along = Cprofiles1),
                     function(ii)
                     cor(Cprofiles1[[ii]][(n1 - 2*overlap + 1):n1],
                         Cprofiles2[[ii]][1:(2*overlap)]))
      corC[is.na(corC)] <- 0 ## e.g. when comparing empty profiles
      if (median(corC) < simCThreshold) ## mean? min? max?
          result[i, i-1] <- result[i-1, i] <- 1
    }
  }
  if (max(result) == 0) { ## confirmed: all equal!
    rep(1, length(winNrs))
  } else {
    cutree(hclust(as.dist(result), method = "complete"), h = .6)
  }
}


