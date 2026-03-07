#' Converts peak list into an ordered peak table.
#' 
#' The function performs a complete linkage clustering of retention times across
#' all samples, and cuts at a height given by the user (which can be understood
#' as the maximal inter-cluster retention time difference) in the simple case
#' based on retention times. Clustering can also incorporate information about
#' spectral similarity using a distance function adapted from Broeckling et al.,
#' 2014:
#' \deqn{e^{-\frac{(1-c_{ij})^2}{2\sigma_r^2}} \cdot e^{-\frac{(1-(t_i-t_j)^2)}{2\sigma_t^2}}}
#' If two peaks from the same sample are assigned to the same cluster, a warning
#' message is printed to the console. These warnings can usually be ignored, but
#' one could also consider reducing the \code{hmax} variable. However, this may 
#' lead to splitting of peaks across multiple clusters. Another option is to
#' filter the peaks by intensity to remove small features.
#' 
#' @name get_peaktable
#' @aliases get_peaktable
#' @importFrom dynamicTreeCut cutreeDynamicTree
#' @importFrom fastcluster hclust
#' @importFrom stats dist cutree as.dist aggregate
#' @importFrom lattice panel.stripplot panel.abline stripplot
#' @importFrom grDevices colorRampPalette 
#' @importFrom scales rescale
#' @importFrom graphics par
#' @param peak_list A `peak_list` object created by \code{\link{get_peaks}},
#' containing a nested list of peak tables: the first level is the
#' sample, and the second level is the spectral wavelength. Every component is
#' described by a \code{data.frame} with a row for each peak and columns
#' containing information on various peak parameters.
#' @param chrom_list A list of chromatographic matrices.
#' @param response Indicates whether peak area or peak height is to be used
#' as intensity measure. Defaults to `area` setting.
#' @param use.cor Logical. Indicates whether to use corrected retention times
#' (\code{rt.cor} column) or raw retention times (\code{rt} column). Unless
#' otherwise specified, the \code{rt.cor} column will be used by default if it 
#' exists in the provided \code{peak_list}.
#' @param hmax Height at which the complete linkage dendrogram will be cut. Can
#' be interpreted as the maximal intercluster retention time difference.
#' @param plot_it Logical. If \code{TRUE}, for every component a strip plot will be
#' shown indicating the clustering.
#' @param ask Logical. Ask before showing new plot? Defaults to \code{TRUE}.
#' @param clust Specify whether to perform hierarchical clustering based on
#' spectral similarity and retention time (\code{sp.rt}) or retention time alone
#' (\code{rt}). Defaults to \code{rt}. The \code{sp.rt} option is experimental
#' and should be used with caution.
#' @param sigma.t Width of gaussian in retention time distance function.
#' Controls weight given to retention time if \code{sp.rt} is selected.
#' @param sigma.r Width of gaussian in spectral similarity function. Controls
#' weight given to spectral correlation if \code{sp.rt} is selected.
#' @param deepSplit Logical. Controls sensitivity to cluster splitting. If
#' \code{TRUE}, function will return more smaller clusters. See documentation for
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for additional information.
#' @param verbose Logical. Whether to print warning when combining peaks into 
#' single time window. Defaults to \code{FALSE}.
#' @param out Specify \code{data.frame} or \code{matrix} as output. Defaults to
#' \code{data.frame}.
#' @md
#' @return The function returns an S3 \code{\link{peak_table}} object, 
#' containing the following elements:
#' * `tab`: The peak table itself -- a \code{data.frame} of intensities in a
#' sample x peak configuration.
#' * `pk_meta`: A \code{data.frame} containing peak meta-data (e.g., the 
#' spectral component, peak number, and average retention time).
#' * `sample_meta`: A \code{data.frame} of sample meta-data. Must be added using
#' \code{\link{attach_metadata}}.
#' * `ref_spectra`: A \code{data.frame} of reference spectra (in a wavelength × 
#' peak configuration). Must be added using \code{\link{attach_ref_spectra}}.
#' * `args`: A vector of arguments given to \code{\link{get_peaktable}} to 
#' generate the peak table.
#' @author Ethan Bass
#' @note This function is adapted from
#' \href{https://github.com/rwehrens/alsace/blob/master/R/getPeakTable.R}{getPeakTable}
#' function in the alsace package by Ron Wehrens.
#' @md
#' @references
#' * Broeckling, C. D., Afsar F.A., Neumann S., Ben-Hur A., and Prenni J.E. 2014.
#' RAMClust: A Novel Feature Clustering Method Enables Spectral-Matching-Based
#' Annotation for Metabolomics Data. \emph{Anal. Chem.} \bold{86}:6812-6817. 
#' \doi{10.1021/ac501530d}.
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. \emph{
#' Metabolomics} \bold{11}:143-154. \doi{10.1007/s11306-014-0683-5}.
#' @examplesIf interactive()
#' data(Sa_pr)
#' pks <- get_peaks(Sa_pr, lambdas = c('210'))
#' get_peaktable(pks, response = "area")
#' @seealso \code{\link{peak_table-class}}, \code{\link{attach_ref_spectra}}, 
#' \code{\link{attach_metadata}}
#' @export get_peaktable

get_peaktable <- function(peak_list, chrom_list, response = c("area", "height"),
                          use.cor = NULL, hmax = 0.2,
                          plot_it = FALSE, ask = plot_it, 
                          clust = c("rt", "sp.rt"), 
                          sigma.t = NULL, sigma.r = 0.5,
                          deepSplit = FALSE, verbose = FALSE,
                          out = c('data.frame', 'matrix')){
  response <- match.arg(response, c("area", "height"))
  clust <- match.arg(clust, c("rt", "sp.rt"))
  out <- match.arg(out, c('data.frame', 'matrix'))
  if (is.null(use.cor)){
    use.cor <- "rt.cor" %in% colnames(peak_list[[1]][[1]])
  }
  rt <- ifelse(use.cor, "rt.cor", "rt")
  start <- ifelse(use.cor, "start.cor", "start")
  end <- ifelse(use.cor, "end.cor", "end")
  if (!inherits(peak_list, "peak_list"))
    stop("Peak list must be of the associated class.")
  if (clust == "sp.rt"){
    if (missing(chrom_list)){
      chrom_list <- get_chrom_list(peak_list)
    } else get_chrom_list(peak_list, chrom_list)
  }
  ncomp <- length(peak_list[[1]]) ## all elements should have the same length
  if (plot_it) {
    opar <- par(ask = ask, no.readonly = TRUE)
    on.exit(par(opar))
    myPalette <- colorRampPalette(c("green", "blue", "purple", "red", "orange"))
  }
  clusterPeaks <- function(comp, pkLst){
    pkLst <- lapply(pkLst, function(x) lapply(x, function(y){
      if (nrow(y) > 0){
      y[!is.na(y[, rt]), , drop = FALSE]
    } else {
      y
    }}))
    xx <- do.call(rbind, sapply(pkLst, function(samp) samp[comp]))
    file.idx <- xx$sample
    pkcenters <- xx[, rt]
    names(pkcenters) <- NULL
    if (length(pkcenters) < 2) 
      return(NULL)
    if (clust == 'rt'){
      pkcenters.hcl <- fastcluster::hclust(dist(pkcenters), method = "complete")
      pkcenters.cl <- cutree(pkcenters.hcl, h = hmax)
    } else if (clust == 'sp.rt'){
        if (is.null(sigma.t)){
          sigma.t <- 0.5 * mean(do.call(rbind, unlist(pkLst, recursive = FALSE))$end - 
                               do.call(rbind, unlist(pkLst, recursive = FALSE))$start)
        }
      ts <- as.numeric(rownames(chrom_list[[1]]))
      sp <- sapply(seq_along(pkcenters), function(i){
        rescale(t(chrom_list[[file.idx[i]]][
          which(elementwise.all.equal(ts, pkcenters[i])),]))
      }, simplify = TRUE)
      cor.matrix <- cor(sp, method = "pearson")
      mint <- abs(outer(unlist(pkcenters), unlist(pkcenters), FUN="-"))
      S <- (exp((-(1 - abs(cor.matrix))^2)/(2*sigma.r^2)))*exp(-(mint^2)/(2*sigma.t^2))
      D <- 1 - S
      linkage <- "average"
      pkcenters.hcl <- fastcluster::hclust(as.dist(D), method = linkage)
      pkcenters.cl <- dynamicTreeCut::cutreeDynamicTree(pkcenters.hcl, 
                                                        maxTreeHeight = hmax, 
                                                        deepSplit = deepSplit, 
                                                        minModuleSize = 2)
      sing <- which(pkcenters.cl == 0)
      pkcenters.cl[sing] <- max(pkcenters.cl) + seq_along(sing)
    }
    vars <- c(rt, start, end, "sd", "width", "tau", "FWHM", "r.squared", "purity")
    vars <- vars[vars %in% colnames(xx)]
    vars.idx <- match(vars, colnames(xx))
    cl.centers <- aggregate(xx[, vars.idx], by = list(pkcenters.cl),
                            FUN = "mean",
                            na.action = "na.pass")[, -1, drop = FALSE]
    ncl <- length(cl.centers[, rt])
    pkcenters.cl <- order(order(cl.centers[, rt]))[pkcenters.cl]
    cl.centers <- cl.centers[order(cl.centers[, rt]),]
    metaInfo <- cbind(lambda = rep(suppressWarnings(
      as.numeric(names(peak_list[[1]])[comp])), ncl),
                      peak = 1:ncl, 
                      round(cl.centers, 2)
                      )
    rownames(metaInfo) <- NULL
    if (plot_it){
      mycols <- myPalette(nrow(cl.centers))
      cl.df <- data.frame(peaks = pkcenters, 
                          files = factor(file.idx), 
                          cluster = pkcenters.cl)
      print(stripplot(files ~ peaks, data = cl.df, 
                      col = mycols[pkcenters.cl], 
                      pch = pkcenters.cl %% 14,
                      xlab = "Retention time", ylab = "",
                      main = paste("Component", comp),
                      panel = function(...) {
                        panel.stripplot(...)
                        panel.abline(v = cl.centers[,rt], col = mycols)
                      }))
    }
    if (verbose & max(clusCount <- table(file.idx, pkcenters.cl)) > 1){
      warning(paste("More than one peak of one injection in the same cluster", 
                paste("for component ", comp, ".", sep = ""), 
                "Keeping only the most intense one.", "", sep = "\n"))
    }
    allIs <- unlist(lapply(pkLst, function(samp) samp[[comp]][, response]))
    Iinfo <- matrix(0, ncl, length(pkLst), dimnames = list(NULL, names(pkLst)))
    for (i in seq(along = allIs)){
      Iinfo[pkcenters.cl[i],  file.idx[i]] <- 
      max(allIs[i], Iinfo[pkcenters.cl[i], file.idx[i]])
    }
    return(list(Iinfo, metaInfo))
  }
  as.structure <- switch(out, "data.frame" = as.data.frame,
             "matrix" = as.matrix)
  result <- lapply(seq_len(ncomp), clusterPeaks, peak_list)
  result <- list(tab = as.structure(t(do.call("rbind", lapply(result,    
                                                        function(x) x[[1]])))),
                 pk_meta = as.structure(t(do.call("rbind", lapply(result, 
                                                        function(x) x[[2]])))),
                 sample_meta = NA,
                 ref_spectra = NA,
                 args = list(peak_list = deparse(substitute(peak_list)),
                        chrom_list = attr(peak_list, "chrom_list"),
                        lambdas = list(names(peak_list[[1]])),
                        response = response,
                        use.cor = use.cor,
                        hmax = hmax,
                        clust = clust,
                        sigma.t = sigma.t,
                        sigma.r = sigma.r,
                        deepSplit = deepSplit,
                        reference_spectra = NA,
                        metadata_path = NA,
                        normalized = FALSE,
                        normalization_by = NA,
                        time_unit = get_metadata_attribute(peak_list, "time_unit"),
                        intensity_unit = get_metadata_attribute(peak_list, "intensity_unit")
                        ))
  class(result) <- c("peak_table", "list")
  attr(result, "pk_args") <- attr(peak_list,"meta")
  result
}

#' Reshape peaktable
#' 
#' Reshapes peak table from wide to long format
#' @name reshape_peaktable
#' @importFrom stats reshape
#' @param x A \code{peak_table} object.
#' @param peaks A character vector specifying the peaks to include. If the
#' character vector is named, the names of the vector elements will be used in
#' place of the original peak names.
#' @param metadata A character vector specifying the metadata fields to include.
#' @param fixed_levels Logical. Whether to fix factor levels of features in the
#' order provided. Defaults to \code{TRUE}.
#' @return A data.frame containing the information for the specified 
#' \code{peaks} in long format.
#' @author Ethan Bass
#' @family utility functions
#' @export
reshape_peaktable <- function(x, peaks, metadata, fixed_levels = TRUE){
  if (!missing(peaks)){
    if (is.numeric(peaks)){
      peaks <- colnames(x$tab)[peaks]
    }
    df <- x$tab[, match(peaks, colnames(x$tab)), drop = FALSE]
    if (!is.null(names(peaks))){
      colnames(df) <- names(peaks)
      peaks <- colnames(df)
    }
  } else {
    df <- x$tab
  }
  if (!missing(metadata)){
    meta_idx <- which(colnames(x$sample_meta) %in% metadata)
    x$sample_meta <- x$sample_meta[, meta_idx, drop = FALSE]
  }
  xx <- reshape(as.data.frame(chr = rownames(df), df), direction = "long",
                varying = list(seq_len(ncol(df))), v.names = x$args[["response"]],
                times = colnames(df), timevar = "peak",
                idvar = "sample", ids = rownames(df))
  rownames(xx) <- NULL
  xx <- merge(xx, data.frame(peak=colnames(x$pk_meta), 
                             t(x$pk_meta[c("lambda", "rt"),])),
              by = "peak", all.x = TRUE)
  xx <- xx[, c(1, 3, 4, 5, 2)]
  if (!is.null(dim(x$sample_meta))){
    xx <- merge(xx, data.frame(sample = row.names(df), x$sample_meta),
                by = "sample", all.x = TRUE)
  }
  if (fixed_levels){
    xx$peak <- factor(xx$peak, levels = peaks)
  }
  xx
}
