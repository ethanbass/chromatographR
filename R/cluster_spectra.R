setClass("cluster", representation(peaks = "character", pval = "numeric"))

#' Function to cluster peaks by spectral similarity.
#' 
#' Function to cluster peaks by spectral similarity. A representative spectrum
#' is selected for each peak in the provided peak table and used to construct a
#' distance matrix based on spectral similarity (pearson correlation) between
#' peaks. Currently, representative spectrum is just selected from the
#' chromatogram with the highest absorbance at lambda max. Hierarchical
#' clustering with bootstrap resampling is performed on the resulting
#' correlation matrix, as implemented in the
#' \code{\link[pvclust:pvclust]{pvclust}} package. Bootstrap values can be used
#' to select clusters that exceed a certain confidence threshold as defined by
#' alpha.
#' 
#' @importFrom pvclust pvclust pvrect pvpick
#' @importFrom stats cor
#' @importFrom methods new
#' @importFrom graphics matplot
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param peak_no Minimum and maximum thresholds for the number of peaks a
#' cluster may have.
#' @param alpha Confidence threshold for inclusion of cluster.
#' @param nboot Number of bootstrap replicates for
#' \code{\link[pvclust:pvclust]{pvclust}}.
#' @param plot_dend Logical. If TRUE, plots dendrogram with bootstrap values.
#' @param plot_spectra Logical. If TRUE, plots overlapping spectra for each
#' cluster.
#' @param verbose Logical. If TRUE, prints progress report to console.
#' @param save Logical. If TRUE, saves pvclust object to current directory.
#' @param parallel Logical. If TRUE, use parallel processing for
#' \code{\link[pvclust:pvclust]{pvclust}}.
#' @param max.only Logical. If TRUE, returns only highest level for nested
#' dendrograms.
#' @param ... Additional arguments to \code{\link[pvclust:pvclust]{pvclust}}.
#' @return Returns S4 "cluster" object with the following components:
#' \item{peaks}{a character vector containing the names of all peaks contained
#' in the given cluster.} \item{pval}{a numeric vector of length 1 containing
#' the bootstrap p-value (au) for the given cluster.}
#' @author Ethan Bass
#' @references R. Suzuki, H. Shimodaira: Pvclust: an R package for assessing
#' the uncertainty in hierarchical clustering. Bioinformatics, 22-12:1540-1542
#' (2006). \doi{10.1093/bioinformatics/btl117}.
#' @examples \donttest{
#' data(pk_tab)
#' data(Sa_warp)
#' cl <- cluster_spectra(pk_tab, nboot=100, max.only = FALSE, save = FALSE, alpha = .97)
#' }
#' @export cluster_spectra

cluster_spectra <- function(peak_table, chrom_list, peak_no = c(5,100),
                            alpha=0.95, nboot=1000, plot_dend=TRUE,
                            plot_spectra=TRUE, verbose=TRUE, save=TRUE,
                            parallel=TRUE, max.only=FALSE,
                            ...){
  if (missing(chrom_list)){
    chrom_list <- try(get(peak_table$args["chrom_list"]))
    if (inherits(chrom_list, "try-error")) stop("Chromatograms not found!")
  }
  if (verbose) print('...collecting representative spectra')
  rep <- sapply(colnames(peak_table[[1]]), function(j){
    sp <- plot_spectrum(loc=j, peak_table=peak_table, chrom_list,
                        scale_spectrum=TRUE, plot_trace=FALSE,
                        export_spectrum=TRUE, plot_spectrum=FALSE, verbose=FALSE)
  })
  rep <- data.frame(do.call(cbind,rep))
  names(rep) <- paste0('V',seq_len(ncol(rep)))
  d<-1-abs(cor(rep,method="pearson"))
  
  if (verbose) print('...clustering spectra')
  result <- pvclust(rep, method.dist="cor",
                             nboot=nboot, parallel=parallel, ...)
  
  if (plot_dend){
  plot(result,labels = FALSE, cex.pv=0.5, print.pv='au',print.num = FALSE)
  pvrect(result, alpha=alpha, max.only = max.only)
  }
  if (save) saveRDS(result, 'pvclust.RDS')
  p <- pvpick(result, alpha=alpha, max.only=max.only)
  l <- sapply(p$clusters, length)
  sub <- p$clusters[which(l > peak_no[1] & l < peak_no[2])]
  pval <- 1-result$edges[p$edges[which(l > peak_no[1] & l < peak_no[2])],'au']
  sub <- lapply(seq_along(sub), function(i){
    new("cluster", peaks=sub[[i]], pval=pval[i])})
  pval <- format(round(
    result$edges[p$edges[which(l > peak_no[1] & l < peak_no[2])],'au'],2),
    nsmall=2)
  names(sub) <- paste0('c',seq_along(sub))
  
  if (plot_spectra){
    if (verbose) print('...plotting clustered spectra')
    new.lambdas <- colnames(chrom_list[[1]])
    sapply(seq_along(sub), function(i){ 
      matplot(new.lambdas,rep[,as.numeric(gsub('V','',sub[[i]]@peaks))],
              type='l', ylab='', yaxt='n', xlab=expression(lambda),
              main=paste0('cluster ', i, '; p = ',
                          format(round(sub[[i]]@pval,2),nsmall=2))
                          )})
  }
  return(sub)
}
