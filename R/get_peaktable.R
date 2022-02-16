get_peaktable <- function(peak_list, chrom_list = NULL, response = c("area", "height"),
                          use.cor = FALSE, hmax = 0.2, plotIt = FALSE,
                          ask = plotIt, clust = c("rt","sp.rt"),
                          sigma.t = 2, sigma.r = 0.5,
                          deepSplit = FALSE, out = c('data.frame', 'matrix')){
  
  response <- match.arg(response, c("area", "height"))
  clust <- match.arg(clust, c("rt","sp.rt"))
  out <- match.arg(out, c('data.frame', 'matrix'))
  rt <- ifelse(use.cor, "rt.cor", "rt")
  ncomp <- length(peak_list[[1]]) ## all elements should have the same length
  if (plotIt) {
    opar <- par(ask = ask, no.readonly = TRUE)
    on.exit(par(opar))
    myPalette <- colorRampPalette(c("green", "blue", "purple", "red", "orange"))
  }
  clusterPeaks <- function(comp, pkLst){
    pkLst <- lapply(pkLst, function(x) lapply(x, function(y) if (nrow(y) > 0){
      y[!is.na(y[, rt]), , drop = FALSE]
    }
    else {
      y
    }))
    file.idx <- rep(names(pkLst), sapply(pkLst, function(samp) nrow(samp[[comp]])))
    pkcenters <- unlist(lapply(pkLst, function(samp) samp[[comp]][, 
                                                                  rt]))
    names(pkcenters) <- NULL
    if (length(pkcenters) < 2) 
      return(NULL)
    if (clust == 'rt'){
    pkcenters.hcl <- hclust(dist(pkcenters), method = "complete")
    pkcenters.cl <- cutree(pkcenters.hcl, h = hmax)
    }
    if (clust == 'sp.rt'){
      if (is.null(chrom_list)){
        stop('Must provide list of chromatograms for spectral clustering.')
      }
      sp <- sapply(1:length(pkcenters), function(i){
        rescale(t(chrom_list[[file.idx[i]]][pkcenters[i],]))
      }, simplify=T)
      c <- cor(sp,sp,method = "pearson")
      mint <- abs(outer(unlist(pkcenters),unlist(pkcenters), FUN="-"))
      S<-exp((-(1-abs(c))^2)/(2*sigma.r^2))*exp(-(mint^2)/2*sigma.t^2)
      D<-1-S
      linkage = "average"
    
      pkcenters.hcl <- hclust(as.dist(D), method = linkage)
      pkcenters.cl <- cutreeDynamicTree(pkcenters.hcl, maxTreeHeight = hmax, deepSplit = deepSplit, minModuleSize = 2)
      sing <- which(pkcenters.cl == 0)
      pkcenters.cl[sing] <- max(pkcenters.cl) + 1:length(sing)
    }
    cl.centers <- aggregate(pkcenters, list(pkcenters.cl), 
                            "mean")[, 2]
    ncl <- length(cl.centers)
    ## reorder clusters from small to large rt
    
    pkcenters.cl <- order(order(cl.centers))[pkcenters.cl]
    cl.centers <- sort(cl.centers)
    metaInfo <- cbind(Component = rep(comp, ncl),
                      Peak = 1:ncl, 
                      RT = cl.centers)
    if (plotIt){
      mycols <- myPalette(length(cl.centers))
      cl.df <- data.frame(peaks = pkcenters, files = factor(file.idx), 
                          cluster = pkcenters.cl)
      message(stripplot(files ~ peaks, data = cl.df, col = mycols[pkcenters.cl], 
                      pch = pkcenters.cl%%14, xlab = "Retention time", 
                      ylab = "", main = paste("Component", comp), panel = function(...) {
                        panel.stripplot(...)
                        panel.abline(v = cl.centers, col = mycols)
                      }))
    }
    if (max(clusCount <- table(file.idx, pkcenters.cl)) > 
        1) 
      Warning(paste("More than one peak of one injection in the same cluster", 
                paste("for component ", comp, ".", sep = ""), 
                "Keeping only the most intense one.", "", sep = "\n"))
    allIs <- unlist(lapply(pkLst, function(samp) samp[[comp]][, 
                                                              response]))
    Iinfo <- matrix(0, ncl, length(pkLst), dimnames = list(NULL, 
                                                           names(pkLst)))
    for (i in seq(along = allIs)) Iinfo[pkcenters.cl[i], 
                                        file.idx[i]] <- max(allIs[i], Iinfo[pkcenters.cl[i], 
                                                                            file.idx[i]])
    cbind(metaInfo, Iinfo)
  }
  result <- lapply(1:ncomp, clusterPeaks, peak_list)
  result <- t(do.call("rbind", result))
  if (out == "data.frame"){
    return(data.frame(result))
  } else return(result)
}

getPeakTable <- function(peak_list, chrom_list = NULL, response = c("area", "height"),
                          use.cor = FALSE, hmax = 0.2, plotIt = FALSE,
                          ask = plotIt, clust = c("rt","sp.rt"),
                          sigma.t = 2, sigma.r = 0.5,
                          deepSplit = FALSE, out = c('data.frame', 'matrix')){
  msg<-"The getPeakTable function is deprecated. Please use get_peaktable instead"
  .Deprecated(get_peaktable, package="chromatographR", msg,
              old = as.character(sys.call(sys.parent()))[1L])
  response <- match.arg(response)
  clust <- match.arg(clust, c('rt','sp.rt'))
  out <- match.arg(out, c('data.frame', 'matrix'))
  rt <- ifelse(use.cor, "rt.cor", "rt")
  ncomp <- length(peak_list[[1]]) ## all elements should have the same length
  if (plotIt) {
    opar <- par(ask = ask, no.readonly = TRUE)
    on.exit(par(opar))
    myPalette <- colorRampPalette(c("green", "blue", "purple", "red", "orange"))
  }
  clusterPeaks <- function(comp, pkLst){
    pkLst <- lapply(pkLst, function(x) lapply(x, function(y) if (nrow(y) > 0){
      y[!is.na(y[, rt]), , drop = FALSE]
    }
    else {
      y
    }))
    file.idx <- rep(names(pkLst), sapply(pkLst, function(samp) nrow(samp[[comp]])))
    pkcenters <- unlist(lapply(pkLst, function(samp) samp[[comp]][, 
                                                                  rt]))
    names(pkcenters) <- NULL
    if (length(pkcenters) < 2) 
      return(NULL)
    if (clust == 'rt'){
      pkcenters.hcl <- hclust(dist(pkcenters), method = "complete")
      pkcenters.cl <- cutree(pkcenters.hcl, h = hmax)
    }
    if (clust == 'sp.rt'){
      if (is.null(chrom_list)){
        stop('Must provide list of chromatograms for spectral clustering.')
      }
      sp <- sapply(1:length(pkcenters), function(i){
        scales::rescale(t(chrom_list[[file.idx[i]]][pkcenters[i],]))
      }, simplify=T)
      c <- cor(sp,sp,method = "pearson")
      mint <- abs(outer(unlist(pkcenters),unlist(pkcenters), FUN="-"))
      S<-exp((-(1-abs(c))^2)/(2*sigma.r^2))*exp(-(mint^2)/2*sigma.t^2)
      D<-1-S
      linkage = "average"
      
      pkcenters.hcl <- hclust(as.dist(D), method = linkage)
      pkcenters.cl <- cutreeDynamicTree(pkcenters.hcl, maxTreeHeight = hmax, deepSplit = deepSplit, minModuleSize = 2)
      sing <- which(pkcenters.cl == 0)
      pkcenters.cl[sing] <- max(pkcenters.cl) + 1:length(sing)
    }
    cl.centers <- aggregate(pkcenters, list(pkcenters.cl), 
                            "mean")[, 2]
    ncl <- length(cl.centers)
    ## reorder clusters from small to large rt
    
    pkcenters.cl <- order(order(cl.centers))[pkcenters.cl]
    cl.centers <- sort(cl.centers)
    metaInfo <- cbind(Component = rep(comp, ncl),
                      Peak = 1:ncl, 
                      RT = cl.centers)
    if (plotIt){
      mycols <- myPalette(length(cl.centers))
      cl.df <- data.frame(peaks = pkcenters, files = factor(file.idx), 
                          cluster = pkcenters.cl)
      message(stripplot(files ~ peaks, data = cl.df, col = mycols[pkcenters.cl], 
                      pch = pkcenters.cl%%14, xlab = "Retention time", 
                      ylab = "", main = paste("Component", comp), panel = function(...) {
                        panel.stripplot(...)
                        panel.abline(v = cl.centers, col = mycols)
                      }))
    }
    if (max(clusCount <- table(file.idx, pkcenters.cl)) > 
        1) 
      warning(paste("More than one peak of one injection in the same cluster", 
                paste("for component ", comp, ".", sep = ""), 
                "Keeping only the most intense one.", "", sep = "\n"))
    allIs <- unlist(lapply(pkLst, function(samp) samp[[comp]][, 
                                                              response]))
    Iinfo <- matrix(0, ncl, length(pkLst), dimnames = list(NULL, 
                                                           names(pkLst)))
    for (i in seq(along = allIs)) Iinfo[pkcenters.cl[i], 
                                        file.idx[i]] <- max(allIs[i], Iinfo[pkcenters.cl[i], 
                                                                            file.idx[i]])
    cbind(metaInfo, Iinfo)
  }
  result <- lapply(1:ncomp, clusterPeaks, peak_list)
  result <- t(do.call("rbind", result))
  if (out == "data.frame"){
    return(data.frame(result))
  } else return(result)
}