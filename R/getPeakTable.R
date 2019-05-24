getPeakTable <- function(peakList, response = c("area", "height"),
                         use.cor = TRUE, maxdiff = .2, plotIt = FALSE,
                         ask = plotIt)
{
  response <- match.arg(response)

  rt <- ifelse(use.cor, "rt.cor", "rt")
  
  ncomp <- length(peakList[[1]]) ## all elements should have the same length
  if (plotIt) {
    if(!require(lattice)) {
      plotIt <- FALSE
      warning("Install the lattice package for the plots...")
    }

    opar <- par(ask = ask, no.readonly = TRUE)
    on.exit(par(opar))
    
    myPalette <- colorRampPalette(c("green", "blue", "purple", "red", "orange"))
  }

  clusterPeaks <- function(comp, pkLst) {
    ## if a peak is lost in the alignment, there will
    ## be an NA... We disregard these cases.
    pkLst <- lapply(pkLst,
                    function(x)
                    lapply(x,
                           function(y)
                           if (nrow(y) > 0) {
                             y[!is.na(y[,rt]), , drop = FALSE]
                           } else {
                             y
                           }))
    
    file.idx <- rep(names(pkLst),
                    sapply(pkLst, function(samp) nrow(samp[[comp]])))
    pkcenters <- unlist(lapply(pkLst,
                               function(samp) samp[[comp]][,rt]))
    names(pkcenters) <- NULL

    if (length(pkcenters) < 2) return(NULL)
    
    pkcenters.hcl <- hclust(dist(pkcenters), method = "complete")
    pkcenters.cl <- cutree(pkcenters.hcl, h = maxdiff)
    cl.centers <-
        aggregate(pkcenters, list(pkcenters.cl), "mean")[,2]
    ncl <- length(cl.centers)
    
    ## reorder clusters from small to large rt
    pkcenters.cl <- order(order(cl.centers))[pkcenters.cl]
    cl.centers <- sort(cl.centers)
    
    metaInfo <- cbind("Component" = rep(comp, ncl),
                      "Peak" = 1:ncl,
                      "RT" = cl.centers)
    
    if (plotIt) {
      mycols <- myPalette(length(cl.centers))
      cl.df <- data.frame(peaks = pkcenters,
                          files = factor(file.idx),
                          cluster = pkcenters.cl)
      ##      dev.new()
      print(stripplot(files ~ peaks, data = cl.df,
                      col = mycols[pkcenters.cl],
                      pch = pkcenters.cl %% 14, ## max pch
                      xlab = "Retention time", ylab = "",
                      main = paste("Component", comp),
                      panel = function(...) {
                        panel.stripplot(...)
                        panel.abline(v = cl.centers, col = mycols)
                      }))
    }

    if (max(clusCount <- table(file.idx, pkcenters.cl)) > 1) 
      warning(paste("More than one peak of one injection in the same cluster",
                    paste("for component ", comp, ".", sep = ""),
                    "Keeping only the most intense one.",
                    "", sep = "\n"))
    
    allIs <- unlist(lapply(pkLst,
                           function(samp) samp[[comp]][,response]))
    Iinfo <- matrix(0, ncl, length(pkLst),
                    dimnames = list(NULL, names(pkLst)))
    for (i in seq(along = allIs))
        Iinfo[pkcenters.cl[i], file.idx[i] ] <-
            max(allIs[i], Iinfo[pkcenters.cl[i], file.idx[i] ])
    
    cbind(metaInfo, Iinfo)
  }

  result <- lapply(1:ncomp, clusterPeaks, peakList)

  do.call("rbind", result)
}
