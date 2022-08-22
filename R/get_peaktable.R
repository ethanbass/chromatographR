#' Convert peak list into an ordered peak table.
#' 
#' Returns a peak_table object. The first slot contains a matrix of
#' intensities, where rows correspond to samples and columns correspond to
#' aligned features. The rest of the slots contain various meta-data about peaks,
#' samples, and experimental settings. 
#' 
#' The function performs a complete linkage clustering of retention times across
#' all samples, and cuts at a height given by the user (which can be understood
#' as the maximal inter-cluster retention time difference) in the simple case
#' based on retention times. Clustering can also incorporate information about
#' spectral similarity using a distance function adapted from Broeckling et al.,
#' 2014:
#' \deqn{latexascii}{e^{-\frac{(1-c_{ij})^2}{2\sigma_r^2}}* e^{-\frac{(1-(t_i-t_j)^2}{2\sigma_t^2}}}
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
#' sample, and the second level is the spectral component. Every component is
#' described by a data.frame where every row is one peak, and the columns contain
#' information on various peak parameters.
#' @param chrom_list A list of chromatographic matrices.
#' @param response Indicates whether peak area or peak height is to be used
#' as intensity measure. Defaults to `area` setting.
#' @param use.cor Logical. Indicates whether to use corrected retention times
#' (by default) or raw retention times (not advised!).
#' @param hmax Height at which the complete linkage dendrogram will be cut. Can
#' be interpreted as the maximal inter-cluster retention time difference.
#' @param plot_it Logical. If TRUE, for every component a stripplot will be
#' shown indicating the clustering.
#' @param ask Logical. Ask before showing new plot?
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
#' \code{\link{cutreeDynamic}} for additional information.
#' @param verbose Logical. Whether to print warning when combining peaks into 
#' single time window. Defaults to \code{FALSE}.
#' @param out Specify `data.frame` or `matrix` as output. Defaults to
#' `data.frame`.
#' @md
#' @return The function returns a `peak_table` object, consisting of the following
#' elements:
#' * `tab`: the peak table itself -- a data-frame of intensities in a
#' sample x peak configuration.
#' * `pk_meta`: A data.frame containing peak meta-data (e.g. the spectral component,
#' peak number, and average retention time).
#' * `sample_meta`: A data.frame of sample meta-data. Must be added using
#' \code{\link{attach_metadata}}).
#' * `ref_spectra`: A data.frame of reference spectra (in a wavelength x peak
#' configuration). Must be added using \code{\link{attach_ref_spectra}}
#' * `args`: A vector of arguments given to \code{\link{get_peaktable}} to generate
#' the peak table.
#' @author Ethan Bass
#' @note Adapted from
#' \href{https://github.com/rwehrens/alsace/blob/master/R/getPeakTable.R}{getPeakTable}
#' function in the alsace package by Ron Wehrens.
#' @md
#' @references
#' * Broeckling, C. D., Afsar F.A., Neumann S., Ben-Hur A., and Prenni J.E. 2014.
#' RAMClust: A Novel Feature Clustering Method Enables Spectral-Matching-Based
#' Annotation for Metabolomics Data. \emph{Anal. Chem.}
#' \bold{86}:6812-6817. \doi{10.1021/ac501530d}
#' * Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R. \emph{
#' Metabolomics} \bold{11}:143-154. \doi{10.1007/s11306-014-0683-5}
#' @examplesIf interactive()
#' data(Sa_pr)
#' pks <- get_peaks(Sa_pr, lambdas = c('210'))
#' get_peaktable(pks, response = "area")
#' @seealso \code{\link{attach_ref_spectra}} \code{\link{attach_metadata}}
#' @export get_peaktable
#' 
get_peaktable <- function(peak_list, chrom_list, response = c("area", "height"),
                          use.cor = FALSE, hmax = 0.2, plot_it = FALSE,
                          ask = plot_it, clust = c("rt","sp.rt"),
                          sigma.t = NULL, sigma.r = 0.5,
                          deepSplit = FALSE, verbose = FALSE,
                          out = c('data.frame', 'matrix')){
  response <- match.arg(response, c("area", "height"))
  clust <- match.arg(clust, c("rt","sp.rt"))
  out <- match.arg(out, c('data.frame', 'matrix'))
  rt <- ifelse(use.cor, "rt.cor", "rt")
  if (!inherits(peak_list,"peak_list"))
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
    pkcenters <- xx$rt
    names(pkcenters) <- NULL
    if (length(pkcenters) < 2) 
      return(NULL)
    if (clust == 'rt'){
    pkcenters.hcl <- hclust(dist(pkcenters), method = "complete")
    pkcenters.cl <- cutree(pkcenters.hcl, h = hmax)
    }
    if (clust == 'sp.rt'){
      if (is.null(sigma.t)){
        sigma.t <- .5*mean(do.call(rbind,unlist(pkLst,recursive = FALSE))$end - 
                             do.call(rbind,unlist(pkLst,recursive = FALSE))$start)
      }
      ts<- as.numeric(rownames(chrom_list[[1]]))
      sp <- sapply(seq_along(pkcenters), function(i){
        rescale(t(chrom_list[[file.idx[i]]][
          which(elementwise.all.equal(ts, pkcenters[i])),]))
      }, simplify=TRUE)
      cor.matrix <- cor(sp, method = "pearson")
      mint <- abs(outer(unlist(pkcenters), unlist(pkcenters), FUN="-"))
      S <- (exp((-(1-abs(cor.matrix))^2)/(2*sigma.r^2)))*exp(-(mint^2)/(2*sigma.t^2))
      D <- 1-S
      linkage <- "average"
      pkcenters.hcl <- hclust(as.dist(D), method = linkage)
      pkcenters.cl <- cutreeDynamicTree(pkcenters.hcl, maxTreeHeight = hmax, 
                                      deepSplit = deepSplit, minModuleSize = 2)
      sing <- which(pkcenters.cl == 0)
      pkcenters.cl[sing] <- max(pkcenters.cl) + seq_along(sing)
    }
    cl.centers <- aggregate(cbind(rt,start,end,sd,tau,FWHM,r.squared) ~ 
                              pkcenters.cl, data=xx, "mean")[,-1]
    ncl <- length(cl.centers$rt)
    ## re-order clusters from small to large rt
    pkcenters.cl <- order(order(cl.centers$rt))[pkcenters.cl]
    cl.centers <- cl.centers[order(cl.centers$rt),]
    metaInfo <- cbind(component = rep(comp, ncl),
                      peak = 1:ncl, 
                      round(cl.centers,2)
                      )
    rownames(metaInfo) <- NULL
    if (plot_it){
      mycols <- myPalette(length(cl.centers))
      cl.df <- data.frame(peaks = pkcenters, files = factor(file.idx), 
                          cluster = pkcenters.cl)
      message(stripplot(files ~ peaks, data = cl.df, col = mycols[pkcenters.cl], 
                      pch = pkcenters.cl%%14, xlab = "Retention time", 
                      ylab = "", main = paste("Component", comp),
                      panel = function(...) {
                        panel.stripplot(...)
                        panel.abline(v = cl.centers, col = mycols)
                      }))
    }
    if (verbose & max(clusCount <- table(file.idx, pkcenters.cl)) > 1) 
      warning(paste("More than one peak of one injection in the same cluster", 
                paste("for component ", comp, ".", sep = ""), 
                "Keeping only the most intense one.", "", sep = "\n"))
    allIs <- unlist(lapply(pkLst, function(samp) samp[[comp]][, response]))
    Iinfo <- matrix(0, ncl, length(pkLst), dimnames = list(NULL, names(pkLst)))
    for (i in seq(along = allIs)){
      Iinfo[pkcenters.cl[i],  file.idx[i]] <- 
      max(allIs[i], Iinfo[pkcenters.cl[i], file.idx[i]])
    }
    return(list(Iinfo, metaInfo))
  }
  result <- lapply(seq_len(ncomp), clusterPeaks, peak_list)
  result <- list(tab = as.data.frame(t(do.call("rbind", lapply(result,    
                                                        function(x) x[[1]])))),
                 pk_meta = as.data.frame(t(do.call("rbind", lapply(result, 
                                                        function(x) x[[2]])))),
                 sample_meta = NA,
                 ref_spectra = NA,
                 args = c(peak_list = deparse(substitute(peak_list)),
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

#' @importFrom utils tail
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

#' @noRd
#' @export
row.names.peak_table <- function(x){
  row.names(x$tab)
}

#' Plot spectrum from peak table
#' 
#' Plots the trace and/or spectrum for a given peak in peak table. 
#' 
#' Can be used to confirm the identity of a peak or check that a particular
#' column in the peak table represents a single compound. Can also be used
#' to create simple box-plots to examine the distribution of a peak with respect
#' to variables defined in sample metadata.
#' 
#' @importFrom scales rescale
#' @importFrom graphics identify title text boxplot
#' @importFrom stats as.formula
#' @param x The peak table (output from \code{\link{get_peaktable}}
#' function).
#' @param ... Additional arguments.
#' @param loc The name of the peak or retention time that you wish to plot.
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param what What to look for. Either \code{peak} to extract spectral information
#' for a certain peak, \code{rt} to scan by retention time, or \code{click} to manually
#' select retention time by clicking on the chromatogram. Defaults to \code{peak}.
#' @param chr Numerical index of chromatogram you wish to plot; "max" to
#' plot the chromatogram with the largest signal; or "all" to plot spectra
#' for all chromatograms.
#' @param lambda The wavelength you wish to plot the trace at (if
#' \code{plot_chrom} is TRUE and/or the wavelength to be used for the determination
#' of signal abundance.
#' @param plot_spectrum Logical. If TRUE, plots the spectrum of the chosen
#' peak. Defaults to TRUE.
#' @param plot_trace Logical. If TRUE, plots the trace of the chosen peak at
#' lambda. Defaults to TRUE.
#' @param box_plot Logical. If TRUE, plots box plot using categories
#' defined by \code{vars}.
#' @param vars Independent variables for boxplot.
#' @param spectrum_labels Logical. If TRUE, plots labels on maxima in spectral
#' plot. Defaults to TRUE.
#' @param scale_spectrum Logical. If TRUE, scales spectrum to unit height.
#' Defaults to FALSE.
#' @param export_spectrum Logical. If TRUE, exports spectrum to console.
#' Defaults to FALSE.
#' @param verbose Logical. If TRUE, prints verbose output to console. Defaults
#' to TRUE.
#' @return If \code{export_spectrum} is TRUE, returns the spectrum as a \code{
#' data.frame} with wavelengths as rows and columns encoding the
#' absorbance (or normalized absorbance, if \code{scale_spectrum} is TRUE) for 
#' the specified sample(s). Otherwise, there is no return value.
#' @section Side effects:
#' If \code{plot_trace} is TRUE, plots the chromatographic trace of the specified
#' chromatogram (\code{chr}), at the specified wavelength (\code{lambda}) with a
#' dotted red line to indicate the retention time given by \code{loc}. The
#' trace is a single column from the chromatographic matrix.
#'
#' If \code{plot_spectrum} is TRUE, plots the spectrum for the specified chromatogram
#' at the specified retention time. The spectrum is a single row from the chromatographic
#' matrix.
#' 
#' If \code{box_plot} is TRUE, produces a \code{\link[graphics]{boxplot}} from the
#' specified peak with groups provided by \code{vars}.
#' @author Ethan Bass
#' @rdname plot.peak_table
#' @export
#'
plot.peak_table <- function(x, ..., loc, chrom_list, what="peak",
                            chr = 'max', lambda = 'max',
                            plot_spectrum = TRUE, plot_trace = TRUE,
                            box_plot = FALSE, vars = NULL,
                            spectrum_labels = TRUE, scale_spectrum = FALSE,
                            export_spectrum = FALSE, verbose = TRUE){
  if (what == "peak" & missing(loc)){
    loc <- readline(prompt="Which peak would you like to plot? \n")
    loc <- gsub('\\"', '', loc)
    loc <- gsub("\\'", "", loc)
    if (!(loc %in% colnames(x$tab)))
      stop("Peak not found.")
  }
  if (plot_spectrum == TRUE | plot_trace == TRUE){
    if (chr == "all"){
      out <- plot_all_spectra(loc, x, chrom_list,
                       plot_spectrum = plot_spectrum,
                       export_spectrum = export_spectrum,
                       verbose = verbose, what = what)
    } else{
      out <- plot_spectrum(loc, x, chrom_list, chr=chr,
                    lambda = lambda, plot_spectrum = plot_spectrum,
                    plot_trace = plot_trace, spectrum_labels = spectrum_labels,
                    scale_spectrum = scale_spectrum,
                    export_spectrum = export_spectrum,
                    verbose = verbose, what = what)
    }
  }
  if (box_plot == T){
    if (!is.data.frame(x$sample_meta)){
      stop("Attach metadata to peak_table to make mirror plot.")
    }
    if (is.null(vars)){
      stop("Must provide independent variable or variables (`var`) for boxplot.")
    }
    if (what != "peak"){
      stop("A peak name must be provided to `loc` to produce a boxplot.")
    }
    boxplot(as.formula(paste("x[['tab']][,loc]", vars, sep="~")),
            data = x$sample_meta,
            main = paste(loc, '\n', 'RT = ', round(x$pk_meta['rt', loc],2)),
            ylab="abs", xlab="", ...)
  }
  if (export_spectrum)
    out
}

#' Make mirror plot from peak table.
#' 
#' Plots chromatograms as a mirror plot.
#' 
#' Can be used to confirm the identity of a peak or check that a particular
#' column in the peak table represents a single compound. Can also be used
#' to create simple box-plots to examine the distribution of a peak with respect
#' to variables defined in sample metadata.
#' 
#' @importFrom graphics matplot legend plot.new plot.window par
#' @importFrom utils head tail
#' @param peak_table The peak table (output from \code{\link{get_peaktable}}
#' function).
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param lambdas The wavelength you wish to plot the traces at.
#' @param var Variable to index chromatograms.
#' @param subset Character vector specifying levels to use (if more than 2 levels
#' are present in \code{var}).
#' @param print_legend Logical. Whether to print legend. Defaults to \code{TRUE}.
#' @param legend_txt Character vector containing labels for legend.
#' @param legend_pos Legend position.
#' @param legend_size Legend size (\code{cex} argument). Default is 1.
#' @param mirror Logical. Whether to plot as mirror or stacked plots.
#' Defaults to \code{TRUE}.
#' @param xlim Numerical vector specifying limits for x axis.
#' @param ylim Numerical vector specifying limits for y axis.
#' @param ... Additional arguments to \code{\link{matplot}} function.
#' @return No return value, called for side effects.
#' @section Side effects:
#' If \code{mirror_plot} is TRUE, plots a mirror plot comparing two treatments
#' defined by \code{var} and \code{subset} (if more than two factors are present
#' in \code{var}).
#'
#' Otherwise, if \code{mirror_plot} is FALSE, the treatments are plotted in two
#' separate panes.
#' @author Ethan Bass
#' @examples
#' data(Sa_warp)
#' data(pk_tab)
#' path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
#' meta <- read.csv(path)
#' pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
#' mirror_plot(pk_tab,lambdas=c("210","260"), var="trt", mirror=TRUE, col=c("green","blue"))
#' @export

mirror_plot <- function(peak_table, chrom_list, lambdas, var, subset = NULL,
                        print_legend = TRUE, legend_txt = NULL, legend_pos = "topright",
                        legend_size = 1, mirror = TRUE, xlim=NULL, ylim=NULL, ...){
  meta <- peak_table$sample_meta
  if (!exists("var"))
    stop("Must provide independent variable or variables for mirror plot.")
  if (!is.data.frame(meta)){
    stop("Metadata must be attached to peak table to make mirror plot.")
  }
  if (!(var %in% colnames(meta))){
    stop(paste(var, "not found in sample meta-data."))
  }
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(peak_table)
  } else get_chrom_list(peak_table, chrom_list)
  new.ts <- round(as.numeric(rownames(chrom_list[[1]])),2)
  fac <- factor(meta[,var])
  if (is.null(subset) & length(levels(fac)) > 2)
    stop("The mirror plot requires two levels. Use the subset argument to 
         specify which levels to use.")
  if (!is.null(subset)){
    if (mean(subset %in% levels(fac)) < 1){
      stop(paste("The specified levels are not present in", var))
    }
  }
  if (is.null(subset)){
    trt1 <- levels(fac)[1]
    trt2 <- levels(fac)[2]
  } else {
    trt1 <- subset[1]
    trt2 <- subset[2]
  }
  set1 <- which(meta[,var] == trt1)
  set2 <- which(meta[,var] == trt2)
  if (print_legend & is.null(legend_txt))
    legend_txt <- c(trt1,trt2)
  if (is.null(xlim))
    xlim <- c(head(new.ts,1),tail(new.ts,1))
  y_max <- max(sapply(chrom_list, function(x){
    max(x[,as.character(min(as.numeric(lambdas)))], na.rm=TRUE)
    }))
  if (mirror){
    if (is.null(ylim))
      ylim <- c(-y_max, y_max)
    plot.new()
    plot.window(xlim = xlim,ylim = ylim)
    for (i in set1){
      matplot(new.ts, chrom_list[[i]][,lambdas], type='l', add=TRUE, ...)
    }
    if (print_legend)
      legend("topright", legend=legend_txt[[1]], cex=legend_size, bty = "n")
    for (i in set2){
      matplot(new.ts, -chrom_list[[i]][,lambdas], type='l', add=TRUE, ...)
    }
    if (print_legend)
      legend("bottomright", legend=legend_txt[[2]], cex=legend_size, bty = "n")
  } else{
    oldpar <- par(no.readonly = TRUE)
    par(mfrow=c(2,1))
    if (is.null(ylim))
      ylim <- c(0,y_max)
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    for (i in set1){
      matplot(new.ts, chrom_list[[i]][,lambdas], type='l', add=TRUE, ...)
    }
    if (print_legend)
      legend(legend_pos, legend=legend_txt[[1]], cex=legend_size, bty = "n")
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    for (i in set2){
      matplot(new.ts, chrom_list[[i]][,lambdas], type='l', add=TRUE, ...)
    }
    if (print_legend)
      legend(legend_pos, legend=legend_txt[[2]], cex=legend_size, bty = "n")
    par(oldpar) #reset par
  }
}
