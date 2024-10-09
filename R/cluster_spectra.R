setClass("cluster", representation(peaks = "character", pval = "numeric"))

#' Cluster peaks by spectral similarity.
#' 
#' Function to cluster peaks by spectral similarity. A representative spectrum
#' is selected for each peak in the provided peak table and used to construct a
#' distance matrix based on spectral similarity (pearson correlation) between
#' peaks. Hierarchical clustering with bootstrap resampling is performed on the 
#' resulting correlation matrix to classify peaks by spectral similarity.
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
#' @importFrom stats cor
#' @importFrom methods new
#' @importFrom graphics matplot
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
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
#' pk_tab <- attach_ref_spectra(pk_tab, Sa_warp, ref = "max.int")
#' cl <- cluster_spectra(pk_tab, nboot = 100, max.only = FALSE, 
#' save = FALSE, alpha = .97)
#' }
#' @export cluster_spectra
#' @md

cluster_spectra <- function(peak_table, peak_no = c(5, 100), alpha = 0.95, 
                            nboot = 1000, plot_dend = TRUE, plot_spectra = TRUE, 
                            verbose = getOption("verbose"), 
                            save = FALSE, parallel = TRUE, max.only = FALSE,
                            output = c("pvclust", "clusters"),
                            ...){
  check_for_pkg("pvclust")
  check_peaktable(peak_table)
  output <- match.arg(output, c("pvclust", "clusters"), several.ok = TRUE)
  if (is.data.frame(peak_table$ref_spectra) | is.matrix(peak_table$ref_spectra)){
    spectra <- peak_table$ref_spectra
  } else {
    stop("Please attach reference spectra (using the `attach_ref_spectra` function) before running `cluster_spectra`.")
  }
  rm <- which(apply(spectra, 2, sd) == 0)
  if (length(rm) > 0){
    if (verbose){
      message(paste0("Removing peaks due to bad spectra: ",
                     paste(sQuote(colnames(spectra)[rm]),collapse=", ")))
    }
    spectra <- spectra[, -rm]
  }
  if (verbose)
    message('...clustering spectra')
  result <- suppressWarnings(pvclust::pvclust(spectra, method.dist = "cor", 
                                              nboot = nboot, parallel = parallel,
                                              quiet = !verbose, ...)
  )
  if (plot_dend){
    plot(result, labels = FALSE, cex.pv = 0.5, print.pv = 'au',
         print.num = FALSE)
    pvclust::pvrect(result, alpha = alpha, max.only = max.only)
  }
  if (save){
    saveRDS(result, 'pvclust.RDS')
  }
  picks <- pvclust::pvpick(result, alpha = alpha, max.only = max.only)
  cl_size <- sapply(picks$clusters, length)
  
  ## filter clusters ##
  cl_idx <- which(cl_size > peak_no[1] & cl_size < peak_no[2])
  if (verbose && length(cl_idx) < length(cl_size)){
    message(paste0("Removing ", length(cl_size) - length(cl_idx),
                   " under- or oversized clusters from results."))
  }
  clusters <- picks$clusters[cl_idx]
  pval <- 1 - result$edges[picks$edges[cl_idx], 'au']
  clusters <- lapply(seq_along(clusters), function(i){
    new("cluster", peaks = clusters[[i]], pval = pval[i])
  })
  if (length(clusters) != 0){
    names(clusters) <- paste0('c', seq_along(clusters))
    
    if (plot_spectra){
      if (verbose) message('...plotting clustered spectra')
      lambdas <- as.numeric(rownames(spectra))
      sapply(seq_along(clusters), function(i){ 
        matplot(lambdas, spectra[, clusters[[i]]@peaks],
                type = 'l', ylab = '', yaxt = 'n', xlab = expression(lambda),
                main = paste0('cluster ', i, '; p = ',
                              format.pval(clusters[[i]]@pval, .015, eps=.001, 
                                          digits=2, nsmall=2))
        )
      })
    }
  }
  pvclust <- result
  mget(output)
}
