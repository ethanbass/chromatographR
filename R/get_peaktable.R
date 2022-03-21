#' Convert peak list into an ordered peak table
#' 
#' Function returns a matrix of intensities, where rows correspond to samples
#' and columns correspond to aligned features. The function performs a complete
#' linkage clustering of retention times across all samples, and cuts at a
#' height given by the user (which can be interpreted as the maximal
#' inter-cluster retention time difference) in the simple case based on
#' retention times. If two peaks from the same sample are assigned to the same
#' cluster, a warning message is given.
#' 
#' The clustering based on spectral similarity used a distance function adapted
#' from Broeckling et al., 2014:
#' \deqn{\exp({-\frac{(1-c_{ij})^2}{2\sigma_r^2}})*\exp({-\frac{(1-(t_i -
#' t_j)^2}{2\sigma_t^2}})} If one sees warnings about peaks from the same
#' sample sharing a cluster label, one option is to reduce the \code{maxdiff}
#' variable - this, however, will increase the number of clusters. Another
#' option is to filter the peaks on intensity: perhaps one of the two peaks in
#' the cluster is only a very small feature.
#' 
#' @aliases get_peaktable getPeakTable
#' @importFrom dynamicTreeCut cutreeDynamicTree
#' @importFrom fastcluster hclust
#' @importFrom stats dist cutree as.dist aggregate
#' @importFrom lattice panel.stripplot panel.abline stripplot
#' @importFrom grDevices colorRampPalette 
#' @importFrom scales rescale
#' @param peak_list A nested list of peak tables: the first level is the
#' sample, and the second level is the component. Every component is described
#' by a matrix where every row is one peak, and the columns contain information
#' on retention time, full width at half maximum (FWHM), peak width, height,
#' and area.
#' @param chrom_list A list of chromatographic matrices.
#' @param response An indicator whether peak area or peak height is to be used
#' as intensity measure. Default is peak area.
#' @param use.cor Logical, indicating whether to use corrected retention times
#' (by default) or raw retention times (not advised!).
#' @param hmax Height at which the complete linkage dendrogram will be cut. Can
#' be interpreted as the maximal inter-cluster retention time difference.
#' @param plot_it Logical. If TRUE, for every component a stripplot will be
#' shown indicating the clustering.
#' @param ask Logical. Ask before showing new plot?
#' @param clust Specify whether to perform hierarchical clustering based on
#' spectral similarity and retention time ("sp.rt") or retention time alone
#' ("rt").
#' @param sigma.t Width of gaussian in retention time distance function.
#' Controls weight of retention time.
#' @param sigma.r Width of gaussian in spectral similarity function. Controls
#' weight of spectral correlation.
#' @param deepSplit Logical. Controls sensitivity to cluster splitting. If
#' TRUE, will return more smaller clusters. See documentation for
#' \code{\link{cutreeDynamic}}.
#' @param out Specify "data.frame" or "matrix" as output. Defaults to
#' `data.frame`.
#' @return The function returns a data frame where the first couple of columns
#' contain meta-information on the features (component, peak, retention time)
#' and the other columns contain the intensities of the features in the
#' individual injections.
#' @author Ethan Bass
#' @note Adapted from getPeakTable function in alsace package by Ron Wehrens.
#' @references Broeckling, C. D., F. A. Afsar, S. Neumann, A. Ben-Hur, and J.
#' E. Prenni. 2014. RAMClust: A Novel Feature Clustering Method Enables
#' Spectral-Matching-Based Annotation for Metabolomics Data. \emph{Anal. Chem.}
#' \bold{86}:6812-6817.
#' @examples
#' 
#' data(Sa)
#' new.ts <- seq(1,38,by=.01) # choose time-points
#' new.lambdas <- seq(200, 400, by = 2) # choose wavelengths
#' dat.pr <- lapply(X=Sa,FUN=preprocess,
#'                  dim1=new.ts,
#'                  dim2=new.lambdas)
#' warping.models <- correct_rt(dat.pr, what = "models", lambdas=c('210','260','360'))
#' warp <- correct_rt(chrom_list=dat.pr, models=warping.models)
#' pks <- get_peaks(warp, lambdas="210")
#' get_peaktable(pks, response = "area")
#' 
#' @export get_peaktable
get_peaktable <- function(peak_list, chrom_list = NULL, response = c("area", "height"),
                          use.cor = FALSE, hmax = 0.2, plot_it = FALSE,
                          ask = plot_it, clust = c("rt","sp.rt"),
                          sigma.t = NULL, sigma.r = 0.5,
                          deepSplit = FALSE, out = c('data.frame', 'matrix')){
  response <- match.arg(response, c("area", "height"))
  clust <- match.arg(clust, c("rt","sp.rt"))
  out <- match.arg(out, c('data.frame', 'matrix'))
  rt <- ifelse(use.cor, "rt.cor", "rt")
  if (class(peak_list) != "peak_list")
    stop("Peak list must be of the associated class.")
  if (is.null(chrom_list)){
    chrom_list <- try(get(attr(peak_list, "chrom_list")))
    if (class(chrom_list)=="try-error") stop("Chromatograms not found")
  }
  ncomp <- length(peak_list[[1]]) ## all elements should have the same length
  if (plot_it) {
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
    pkcenters <- unlist(lapply(pkLst, function(samp) samp[[comp]][,rt]))
    names(pkcenters) <- NULL
    if (length(pkcenters) < 2) 
      return(NULL)
    if (clust == 'rt'){
    pkcenters.hcl <- hclust(dist(pkcenters), method = "complete")
    pkcenters.cl <- cutree(pkcenters.hcl, h = hmax)
    }
    if (clust == 'sp.rt'){
      if (is.null(sigma.t)){
        sigma.t <- .5*mean(do.call(rbind,unlist(pkLst,recursive = F))$end - do.call(rbind,unlist(pkLst,recursive = F))$start)
      }
      ts<- as.numeric(rownames(chrom_list[[1]]))
      sp <- sapply(seq_along(pkcenters), function(i){
        rescale(t(chrom_list[[file.idx[i]]][which(elementwise.all.equal(ts, pkcenters[i])),]))
      }, simplify=T)
      cor.matrix <- cor(sp, method = "pearson")
      mint <- abs(outer(unlist(pkcenters), unlist(pkcenters), FUN="-"))
      S <- (exp((-(1-abs(cor.matrix))^2)/(2*sigma.r^2)))*exp(-(mint^2)/(2*sigma.t^2))
      D <- 1-S
      linkage <- "average"
      pkcenters.hcl <- hclust(as.dist(D), method = linkage)
      pkcenters.cl <- cutreeDynamicTree(pkcenters.hcl, maxTreeHeight = hmax, deepSplit = deepSplit, minModuleSize = 2)
      sing <- which(pkcenters.cl == 0)
      pkcenters.cl[sing] <- max(pkcenters.cl) + seq_along(sing)
    }
    cl.centers <- aggregate(pkcenters, list(pkcenters.cl), 
                            "mean")[, 2]
    ncl <- length(cl.centers)
    ## reorder clusters from small to large rt
    pkcenters.cl <- order(order(cl.centers))[pkcenters.cl]
    cl.centers <- sort(cl.centers)
    metaInfo <- cbind(Component = rep(comp, ncl),
                      Peak = 1:ncl, 
                      RT = round(cl.centers,2)
                      )
    if (plot_it){
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
    #cbind(metaInfo, Iinfo)
    return(list(Iinfo, metaInfo))
  }
  result <- lapply(1:ncomp, clusterPeaks, peak_list)
  result <- t(do.call("rbind", result))
  result <- list(tab=data.frame(t(result[[1]])),
                 pk_meta=data.frame(t(result[[2]])),
                 sample_meta=NA,
                 args=c(peak_list = deparse(substitute(peak_list)),
                        chrom_list = attr(peak_list,"chrom_list"),
                        response = response,
                        use.cor = use.cor,
                        hmax = hmax,
                        clust = clust,
                        sigma.t = sigma.t,
                        sigma.r = sigma.r,
                        deepSplit = deepSplit
                        ))
  class(result) <- "peak_table"
  return(result)
}
#' @importFrom utils head
#' @noRd
#' @rdname head.peak_table
#' @export
head.peak_table <- function(x,...){
  head(x$tab)
}

#' @noRd
#' @export
tail.peak_table <- function(x,...){
  tail(x$tab)
}

#' @noRd
#' @export
print.peak_table <- function(x, ...){
  print(x$tab)
}

#' @noRd
#' @export
dim.peak_table <- function(x){
  dim(x$tab)
}

