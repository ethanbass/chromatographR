#' Converts peak list into an ordered peak table.
#' 
#' The function performs a complete linkage clustering of retention times across
#' all samples, and cuts at a height given by the user (which can be understood
#' as the maximal inter-cluster retention time difference) in the simple case
#' based on retention times. Optionally, clustering can also incorporate 
#' information about spectral similarity.
#' 
#' By default, clustering is performed on retention times only (when 
#' `clust == "rt"`). Clustering can also incorporate information about
#' spectral similarity (when `clust == "sp.rt"`) using a distance function 
#' adapted from Broeckling et al., 2014:
#' \deqn{e^{-\frac{(1-c_{ij})^2}{2\sigma_r^2}} \cdot e^{-\frac{(1-(t_i-t_j)^2)}{2\sigma_t^2}}}
#' where \eqn{c_{ij}} is the spectral correlation coefficient between spectra
#' \eqn{i} and \eqn{j}, and \eqn{\sigma_t} and \eqn{\sigma_r} control the
#' relative contributions of retention time and spectral similarity,
#' respectively (see `sigma.t` and `sigma.r`).
#' 
#' If two peaks from the same sample are assigned to the same cluster, a warning
#' message is printed to the console. These warnings can usually be ignored, but
#' one could also consider reducing the `hmax` variable. However, this may 
#' lead to splitting of peaks across multiple clusters. Another option is to
#' filter the peaks by intensity to remove small features.
#' 
#' Once clusters are formed, peak metadata is summarized across all peaks within
#' each cluster according to `summarize_by`. By default, the peak metadata are 
#' summarized as the mean within each cluster weighted by the relative size of
#' each peak. This downweights the importance of small peaks that are more 
#' likely to represent noise, resulting in a more robust estimate of the cluster
#' center. Alternatively, if the `"max"` option is selected, the metadata
#' associated with the most intense peak in each cluster will be returned here.
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
#' @param peak_list A `peak_list` object created by [`get_peaks`], containing a
#' nested list of features, where the first level is the sample, and the second 
#' level is the spectral wavelength. Every sample x wavelength component is 
#' described by a  `data.frame` with a row for each peak and columns containing 
#' information on various peak parameters.
#' @param chrom_list A list of chromatographic matrices.
#' @param response Indicates whether peak `area` or peak `height` is to be used
#' as intensity measure. Defaults to `area` setting.
#' @param use.cor Logical. Indicates whether to use corrected retention times
#' (`rt.cor` column) or raw retention times `rt` column). Unless
#' otherwise specified, the `rt.cor` column will be used by default if it 
#' exists in the provided `peak_list`.
#' @param hmax Height at which the complete linkage dendrogram will be cut. Can
#' be interpreted as the maximal intercluster retention time difference.
#' @param summarize_by How to select the representative peak from each cluster.
#' Options are `"weighted.mean"` (default, weights peaks by their intensity as
#' given by `response`), `"mean"`, `"median"` (which aggregate metadata across 
#' all peaks in the cluster) or `"max"` (which selects the most intense peak in 
#' the cluster and uses its metadata directly).
#' @param plot_it Logical. If `TRUE`, for every component a strip plot will be
#' shown indicating the clustering.
#' @param ask Logical. Ask before showing new plot? Defaults to `TRUE`.
#' @param clust Specify whether to perform hierarchical clustering based on
#' spectral similarity and retention time (`sp.rt`) or retention time alone
#' (`rt`). Defaults to `rt`.
#' @param sigma.t Width of gaussian in retention time distance function.
#' Controls weight given to retention time if `sp.rt` is selected.
#' @param sigma.r Width of gaussian in spectral similarity function. Controls
#' weight given to spectral correlation if `sp.rt` is selected.
#' @param deepSplit Logical. Controls sensitivity to cluster splitting. If
#' `TRUE`, function will return more smaller clusters. See documentation for
#' [`cutreeDynamic`][dynamicTreeCut::cutreeDynamic] for additional information.
#' @param verbose Logical. Whether to print warning when combining peaks into 
#' single time window. Defaults to `FALSE`.
#' @param out Specify `data.frame` (default) or `matrix` as output.
#' @md
#' @return The function returns an S3 [`peak_table`] object, containing the
#' following elements:
#' * `tab`: The peak table itself -- a `data.frame` of intensities in a
#' sample x peak configuration.
#' * `pk_meta`: A `data.frame` containing peak meta-data (e.g., the spectral
#' component, peak number, and average retention time), summarized across all
#' peaks in each cluster according to `summarize_by`.
#' * `sample_meta`: A `data.frame` of sample meta-data. Must be added using
#' [`attach_metadata`].
#' * `ref_spectra`: A `data.frame` of reference spectra (in a wavelength × 
#' peak configuration). Must be added using [`attach_ref_spectra`].
#' * `args`: A vector of arguments given to [`get_peaktable`] to 
#' generate the peak table.
#' @author Ethan Bass
#' @note This function is adapted from the `getPeakTable` function in the alsace
#' package by Ron Wehrens: <https://github.com/rwehrens/alsace/blob/master/R/getPeakTable.R>.
#' @md
#' @references
#' * Broeckling, C. D., Afsar F.A., Neumann S., Ben-Hur A., and Prenni J.E. 2014.
#' RAMClust: A Novel Feature Clustering Method Enables Spectral-Matching-Based
#' Annotation for Metabolomics Data. *Anal. Chem.* 86:6812-6817. 
#' \doi{10.1021/ac501530d}.
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. 
#' *Metabolomics* 11:143-154. \doi{10.1007/s11306-014-0683-5}.
#' @examplesIf interactive()
#' data(Sa_pr)
#' pks <- get_peaks(Sa_pr, lambdas = c('210'))
#' get_peaktable(pks, response = "area")
#' @seealso [`peak_table-class`], [`attach_ref_spectra`], [`attach_metadata`]
#' @export get_peaktable

get_peaktable <- function(peak_list, chrom_list, response = c("area", "height"),
                          use.cor = NULL, hmax = 0.2, 
                          summarize_by = c("weighted.mean", "median", "mean", "max"),
                          plot_it = FALSE, ask = plot_it, 
                          clust = c("rt", "sp.rt"), 
                          sigma.t = NULL, sigma.r = 0.5,
                          deepSplit = FALSE, verbose = FALSE,
                          out = c('data.frame', 'matrix')){
  response <- match.arg(response, c("area", "height"))
  clust <- match.arg(clust, c("rt", "sp.rt"))
  out <- match.arg(out, c('data.frame', 'matrix'))
  summarize_by <- match.arg(summarize_by, choices = c("weighted.mean", 
                                                      "median", "mean", "max"))
  if (is.null(use.cor)){
    use.cor <- "rt.cor" %in% colnames(peak_list[[1]][[1]])
  }
  rt <- ifelse(use.cor, "rt.cor", "rt")
  start <- ifelse(use.cor, "start.cor", "start")
  end <- ifelse(use.cor, "end.cor", "end")
  if (!inherits(peak_list, "peak_list"))
    stop("Peak list must be of the associated class.")
  check_duplicated_names(names(peak_list))
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
    if (summarize_by == "weighted.mean"){
      weights <- xx[, response]
      cl.centers <- do.call(rbind,
                            tapply(seq_along(pkcenters),
                                   pkcenters.cl, function(i){
                                     vapply(vars.idx, function(v){
                                       stats::weighted.mean(xx[i, v], w = weights[i])
                                     }, numeric(1)) |>
                                       setNames(colnames(xx[, vars.idx]))
                                   }))
    } else if (summarize_by == "max"){
      cl.centers <- do.call(rbind, 
                            tapply(seq_along(pkcenters), 
                                   pkcenters.cl, function(i){
        xx[i[which.max(xx[i, response])], vars.idx]
      }))
    } else{
      cl.centers <- aggregate(xx[, vars.idx], by = list(pkcenters.cl),
                              FUN = summarize_by,
                              na.action = "na.pass")[, -1, drop = FALSE]
    }
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
    allIs[!is.finite(allIs) | allIs < 0] <- 0
    Iinfo <- matrix(0, ncl, length(pkLst), dimnames = list(NULL, names(pkLst)))
    for (i in seq_along(allIs)){
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
                        summarized_by = summarize_by,
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
#' @param x A `peak_table` object.
#' @param peaks A character vector specifying the peaks to include. If the
#' character vector is named, the names of the vector elements will be used in
#' place of the original peak names.
#' @param metadata A character vector specifying the metadata fields to include.
#' @param fixed_levels Logical. Whether to fix factor levels of features in the
#' order provided. Defaults to `TRUE`.
#' @return A data.frame containing the information for the specified peaks in 
#' long format.
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
  xx <- merge(xx, data.frame(peak = colnames(x$pk_meta), 
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
