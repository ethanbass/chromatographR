setClass("cluster", representation(peaks = "character", pval = "numeric"))

#' Cluster peaks by spectral similarity.
#' 
#' Function to cluster peaks by spectral similarity. A representative spectrum
#' is selected for each peak in the provided peak table and used to construct a
#' distance matrix based on spectral similarity (pearson correlation) between
#' peaks. Hierarchical clustering with bootstrap resampling is performed on the resulting
#' correlation matrix to classify peaks into by their spectral similarity.
#'
#' A representative spectrum is selected for each peak in the provided peak table
#' and used to construct a distance matrix based on spectral similarity
#' (pearson correlation) between peaks. It is suggested to attach representative
#' spectra to the \code{peak_table} using \code{\link{attach_ref_spectra}}.
#' Otherwise, representative spectra are obtained from the chromatogram with the
#' highest absorbance at lambda max.
#'
#' Hierarchical clustering with bootstrap
#' resampling is performed on the resulting correlation matrix, as implemented in
#' \code{\link[pvclust:pvclust]{pvclust}}. Finally, bootstrap values can be used
#' to select clusters that exceed a certain confidence threshold as defined by
#' \code{alpha}. Clusters can also be filtered by the minimum and maximum
#' size of the cluster using the argument \code{peak_no}. If \code{max_only}
#' is TRUE, only the largest cluster in a nested dendrogram of clusters meeting
#' the confidence threshold will be returned.
#'
#' @name cluster_spectra
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
#' @param output What to return. Either \code{clusters} to return list of clusters,
#' \code{pvclust} to return pvclust object, or \code{both} to return both items.
#' @param ... Additional arguments to \code{\link[pvclust:pvclust]{pvclust}}.
#' @return Returns clusters and/or \code{pvclust} object according to the value
#' of the \code{output} argument.
#' * If \code{output = clusters},  returns a list of S4 \code{cluster} objects.
#' * If \code{output = pvclust}, returns a \code{\link[pvclust:pvclust]{pvclust}}
#' object.
#' * If \code{output = both}, returns a nested list containing \code{[[1]]} the
#' \code{\link[pvclust:pvclust]{pvclust}} object, and \code{[[2]]} the list of
#' S4 \code{cluster} objects.
#'
#' The \code{cluster} objects consist of the following components:
#' * \code{peaks}: a character vector containing the names
#' of all peaks contained in the given cluster.
#' * \code{pval}: a numeric vector of length 1 containing
#' the bootstrap p-value (au) for the given cluster.
#' @author Ethan Bass
#' @references R. Suzuki & H. Shimodaira. 2006. Pvclust: an R package for assessing
#' the uncertainty in hierarchical clustering. \emph{Bioinformatics},
#' \bold{22(12)}:1540-1542. \doi{10.1093/bioinformatics/btl117}.
#' @note
#' * Users should be aware that the clustering algorithm will often return nested
#' clusters. Thus, an individual peak could appear in more than one cluster.
#' * It is highly suggested to use more than 100 bootstraps if you run the 
#' clustering algorithm on real data even though we use \code{nboot = 100} in
#' the example to reduce runtime. The authors of \code{pvclust} suggest \code{
#' nboot = 10000}.
#' @examples \donttest{
#' data(pk_tab)
#' data(Sa_warp)
#' cl <- cluster_spectra(pk_tab, nboot=100, max.only = FALSE, save = FALSE, alpha = .97)
#' }
#' @export cluster_spectra
#' @md

cluster_spectra <- function(peak_table, chrom_list, peak_no = c(5,100),
                            alpha=0.95, nboot=1000, plot_dend=TRUE,
                            plot_spectra=TRUE, verbose=TRUE, save=TRUE,
                            parallel=TRUE, max.only=FALSE,
                            output=c("clusters", "pvclust", "both"),
                            ...){
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(peak_table)
  } else get_chrom_list(peak_table, chrom_list)
  output <- match.arg(output, c("clusters","pvclust","both"))
  if (is.data.frame(peak_table$ref_spectra) | is.matrix(peak_table$ref_spectra)){
    rep <- peak_table$ref_spectra
  } else{
    if (verbose)
      print('...collecting representative spectra')
    rep <- sapply(colnames(peak_table[[1]]), function(j){
      #print(j)
      sp <- plot_spectrum(loc=j, peak_table=peak_table, chrom_list = chrom_list,
                          scale_spectrum=TRUE, plot_trace=FALSE,
                          export_spectrum=TRUE, plot_spectrum=FALSE, verbose=FALSE)
    })
    rep <- data.frame(do.call(cbind,rep))
    names(rep) <- paste0('V',seq_len(ncol(rep)))
  }
  rm <- which(apply(rep,2,sd)==0)
  if (length(rm)>0)
    rep <- rep[,-rm]
  d <- 1 - abs(cor(rep, method="pearson"))
  if (verbose)
    print('...clustering spectra')
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
      matplot(new.lambdas,rep[,sub[[i]]@peaks],
              type='l', ylab='', yaxt='n', xlab=expression(lambda),
              main=paste0('cluster ', i, '; p = ',
                          format(round(sub[[i]]@pval,2),nsmall=2))
                          )})
  }
  switch(output, "clusters" = sub, "pvclust" = result, "both" = list(result,sub))
}
