# peakList <- pks
 response = "area"
 use.cor = F
 maxdiff = 0.2
# plotIt = FALSE
# ask = plotIt
 rt<-"rt"
 sigma.t = 2
 sigma.r = 0.5

require(fastcluster)
require(dynamicTreeCut)

# group spectra by similarity
# cosine distance versus pearson correlation?
# molecular networking versus chemical dendrogram.

peakList <- pkTab
pkLst <- pks.lvs
result <- lapply(1:3, clusterPeaks, peakList)
pkLst<-pks.lvs

cluster_spectra <- function(peakList, hmax = 0.2,
                         deepSplit = FALSE, out = c('data.frame', 'matrix')){
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

      sp <- sapply(1:length(pkcenters), function(i){
        scales::rescale(t(dat.pr[[file.idx[i]]][pkcenters[i],]))
      }, simplify=T)
      
      c <- cor(sp,sp,method = "pearson")
      mint <- abs(outer(unlist(pkcenters),unlist(pkcenters), FUN="-"))
      S <- exp((-(1-abs(c))^2)/(2*sigma.r^2))*exp(-(mint^2)/2*sigma.t^2)
      D<-1-S
      linkage = "average"
      
      pkcenters.hcl <- fastcluster::hclust(as.dist(D), method = linkage)
      # or maybe choose a representative chromatogram from each peak instead
      pkcenters.cl <- dynamicTreeCut::cutreeDynamicTree(pkcenters.hcl, maxTreeHeight = hmax, deepSplit = deepSplit, minModuleSize = 2)
      sing <- which(pkcenters.cl == 0)
      pkcenters.cl[sing] <- max(pkcenters.cl) + 1:length(sing)
      
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
      print(stripplot(files ~ peaks, data = cl.df, col = mycols[pkcenters.cl], 
                      pch = pkcenters.cl%%14, xlab = "Retention time", 
                      ylab = "", main = paste("Component", comp), panel = function(...) {
                        panel.stripplot(...)
                        panel.abline(v = cl.centers, col = mycols)
                      }))
    }
    if (max(clusCount <- table(file.idx, pkcenters.cl)) > 
        1) 
      cat(paste("Warning!", "More than one peak of one injection in the same cluster", 
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
  result <- lapply(1:ncomp, clusterPeaks, peakList)
  result <- t(do.call("rbind", result))
  if (out == "data.frame"){
    return(data.frame(result))
  } else return(result)
}

# tab <- getPeakTable(pks, response = "area",
#                           use.cor = F, hmax = 0.2, plotIt = T,
#                           ask = plotIt, clust="sp.rt",
#                           sigma.t=2, sigma.r=0.5)
# 
# tab <- getPeakTable(pks, response = "area",
#                     use.cor = F, hmax = 0.2, plotIt = T,
#                     ask = plotIt, clust="rt",
#                     sigma.r=1, sigma.t=1)

#?alsace_dev:getPeakTable

###########

# peakList <- pks
response = "area"
use.cor = F
maxdiff = 0.2
# plotIt = FALSE
# ask = plotIt
rt<-"rt"
sigma.t = 2
sigma.r = 0.5

require(fastcluster)
require(dynamicTreeCut)

# group spectra by similarity
# cosine distance versus pearson correlation?
# molecular networking versus chemical dendrogram.

peakList <- pkTab
pkLst <- pks.lvs
result <- lapply(1:3, clusterPeaks, peakList)
pkLst<-pks.lvs
par(mfrow=c(2,1))
plot_spectrum(peak=20, peak_table=pkTab,chrom_list = dat.pr)

cluster_spectra2 <- function(peakList, hmax = 0.2,
                            deepSplit = FALSE, out = c('data.frame', 'matrix')){
  rep <- sapply(1:ncol(dat.pr[,1]), function(j){
    sp <- plot_all_spectra(peak=j, peak_table=pkTab, chrom_list = dat.pr,
                           scale_spectrum=F)
    #sp.s <- apply(sp,2,scales::rescale)
    sp.m <- apply(sp,2,max)
    # which.max(apply(abs(cor(sp.s))*outer(sp.m,sp.m,FUN='/'),2,mean))
    # which.max(apply(cor(sp),2,mean))
    # which.max(sp.m)
    # plot(sp[,38],type='l')
    # plot(sp[,26],type='l')
    ch<-which.max(sp.m)
    #scales::rescale(outer(sp.m,sp.m,FUN='/'))
    #apply(outer(sp.m,sp.m,FUN='/'),2,scales::rescale)
      #height should be factored in here
    # c <- cor(sp,sp,method = "pearson")
    # mint <- abs(outer(unlist(pkcenters),unlist(pkcenters), FUN="-"))
    # S <- exp((-(1-abs(c))^2)/(2*sigma.r^2))*exp(-(mint^2)/2*sigma.t^2)
    # D<-1-S
    # linkage = "average"
    #ch<-which.max(apply(cor(sp),2,mean))
    df <- data.frame(sp[,ch])
    df
})

  rep <- data.frame(do.call(cbind,rep))
  names(rep)<-paste0('V',1:ncol(rep))
  rep.s <- apply(rep,2,scales::rescale)
  d<-1-abs(cor(rep.s,method="pearson"))

  cl<-fastcluster::hclust(as.dist(d), method = linkage)

  plot(cl)
  matplot(rep.s[,c(7,3,9,2)],type='l')
  matplot(rep.s[,c(4,6)],type='l')
  matplot(rep.s[,c(1,5)],type='l')
  matplot(rep.s[,c(10)],type='l',add=T)
}


########
require(pvclust)
require(fastcluster)
#require(dynamicTreeCut)

cluster_spectra <- function(pkTab, chrom_list, 
                            deepSplit = FALSE, mn=5, mx=100,
                            alpha=0.95, nboot=1000, plot_dend=T, plot_spectra=T,
                            verbose=T, parallel=T){
  if (verbose==T) print('...collecting representative spectra')
  rep <- sapply(1:ncol(pkTab), function(j){
    sp <- plot_spectrum(peak=j, peak_table=pkTab, chrom_list = dat.pr,
                        scale_spectrum=T, plot_trace=F, export_spectrum = T, plot_spectrum=F,verbose=F)
  })
  rep <- data.frame(do.call(cbind,rep))
  names(rep)<-paste0('V',1:ncol(rep))
  d<-1-abs(cor(rep,method="pearson"))
  
  #cl<-fastcluster::hclust(as.dist(d), method = linkage)
  #plot(cl)
  if (verbose==T) print('...clustering spectra')
  result <- pvclust(rep, method.dist="cor", method.hclust="average", nboot=nboot, parallel=parallel)
  if (plot_dend==T){
  pvrect(result, alpha=alpha, max.only = F)
  }
  #saveRDS(result, 'pvclust.RDS')
  p <- pvpick(result, alpha=alpha, max.only=F)
  l <- sapply(p$clusters, length)
  sub <- p$clusters[which(l > mn & l < mx)]
  pval<-result$edges[p$edges[which(l > mn & l < mx)],'au']
  sub <- lapply(1:length(sub), function(i) new("cluster", peaks=sub[[i]], pval=pval[i]))
  pval=format(round(result$edges[p$edges[which(l > mn & l < mx)],'au'],2), nsmall=2)
  names(sub)<-paste0('c',1:length(sub))
  if (plot_spectra==T){
    if (verbose==T) print('...plotting clustered spectra')
    sapply(1:length(sub), function(i){ 
      matplot(new.lambdas,rep[,as.numeric(gsub('V','',sub[[i]]@peaks))],
              type='l', ylab='', yaxt='n', xlab=expression(lambda),
              main=paste0('cluster ', i, '; p = ', format(round(sub[[i]]@pval,2),nsmall=2))
                          )})
  }
  return(sub)
}

cluster_spectra(pkTab=pkTab, chrom_list=dat.pr)

setClass("cluster", representation(peaks = "character", pval = "numeric"))
New("cluster", membership=sub[[i]], pval=p)
pvclust_show_signif(dend = as.dendrogram(result$hclust), pvclust_obj = result)
pvclust_edges(result)