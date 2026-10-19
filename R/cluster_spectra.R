setClass("cluster", representation(peaks = "character", pval = "numeric"))

#' Cluster spectra
#' 
#' Cluster peaks by spectral similarity.
#' 
#' Before using this function, reference spectra must be attached to the 
#' `peak_table` using the [`attach_ref_spectra`] function. These spectra are 
#' then used to construct a distance matrix from pairwise Pearson correlations
#' among peaks, and [`pvclust`][pvclust::pvclust] clusters the peaks
#' hierarchically, with bootstrap resampling. Clusters whose bootstrap support
#' exceeds the confidence set by `alpha` are returned. Peaks whose reference
#' spectrum is constant are dropped before clustering.
#' 
#' Clusters can be filtered by size using the `min_size` and `max_size`
#' arguments. The clustering algorithm often returns nested clusters, so an
#' individual peak can appear in more than one. If `max.only` is `TRUE`, only
#' the largest cluster in a nested tree of clusters meeting the confidence
#' threshold is returned.
#' 
#' The example uses `nboot = 100` to reduce run time. For real data, the
#' `pvclust` authors recommend the default of 1000 replicates for a first look,
#' then 10,000 for smaller errors.
#' 
#' @importFrom stats cor
#' @importFrom methods new
#' @importFrom graphics matplot rect text
#' @importFrom grDevices palette.colors
#' @param peak_table A `peak_table` object from [`get_peaktable`].
#' @param alpha Significance threshold used to select clusters based on 
#' approximately unbiased (AU) bootstrap p-values from
#' [`pvclust`][pvclust::pvclust].
#' @param min_size Clusters must contain more than `min_size` peaks. Defaults
#' to `5`.
#' @param max_size Clusters must contain fewer than `max_size` peaks. Defaults
#' to the number of peaks in `peak_table`.
#' @param nboot Number of bootstrap replicates for [`pvclust`][pvclust::pvclust].
#' @param plot_dend Logical. If `TRUE`, plots dendrogram with bootstrap values.
#' @param plot_spectra Logical. If `TRUE`, plots overlapping spectra for each
#' cluster.
#' @param verbose Logical. If `TRUE`, prints progress report to console.
#' @param save The `save` argument is deprecated. Please use `outfile`
#' instead.
#' @param outfile Optional path to an `.RDS` file. If not `NULL`, the
#' `pvclust` object is saved to file using `saveRDS()`. Existing files are
#' overwritten.
#' @param parallel Logical. If `TRUE`, use parallel processing for 
#' [`pvclust`][pvclust::pvclust].
#' @param max.only Logical. If `TRUE`, returns only highest level cluster from
#' each nested set of clusters. Named to match the argument of the same name in 
#' [`pvrect`][pvclust::pvrect] and [`pvpick`][pvclust::pvpick]; passed to `pvpick`.
#' @param output What to return: `"pvclust"` for the `pvclust` object,
#' `"clusters"` for the list of clusters, or both (the default).
#' @param col Colors for each cluster's box and number on the dendrogram and
#' the title of its spectra panel, recycled if there are more clusters than
#' colors. Defaults to the Okabe-Ito palette, without black, yellow and grey.
#' @param ... Additional arguments to [`pvclust`][pvclust::pvclust].
#' @return Returns a named list with an element for each value of `output`:
#' `pvclust`, the `pvclust` object, and `clusters`, a list of S4 `cluster`
#' objects.
#'
#' The `cluster` objects consist of the following components:
#' * `peaks`: a character vector containing the names of all peaks contained in
#' the given cluster.
#' * `pval`: a numeric scalar containing one minus the approximately unbiased
#' (AU) value that [`pvclust`][pvclust::pvclust] reports for the given cluster,
#' so that small values indicate strong support.
#' @author Ethan Bass
#' @references R. Suzuki & H. Shimodaira. 2006. Pvclust: an R package for assessing
#' the uncertainty in hierarchical clustering. *Bioinformatics*,
#' 22(12):1540-1542. \doi{10.1093/bioinformatics/btl117}.
#' @examples \donttest{
#' if (requireNamespace("pvclust", quietly = TRUE)) {
#' data(pk_tab)
#' data(Sa_warp)
#' pk_tab <- attach_ref_spectra(pk_tab, Sa_warp, ref = "max.int")
#' cl <- cluster_spectra(pk_tab, nboot = 100, max.only = FALSE, alpha = 0.03)
#' }
#' }
#' @export cluster_spectra
#' @md

cluster_spectra <- function(peak_table, alpha = 0.05, 
                            min_size = 5, max_size = NULL,
                            nboot = 1000, plot_dend = TRUE, plot_spectra = TRUE, 
                            verbose = getOption("verbose"), outfile = NULL,
                            save = NULL, parallel = TRUE, max.only = FALSE,
                            output = c("pvclust", "clusters"),
                            col = palette.colors(palette = "Okabe-Ito")[-c(1, 5, 9)],
                            ...){
  check_for_pkg("pvclust")
  check_peaktable(peak_table)
  if (!is.null(save)){
    if (save){
      warning("The `save` argument is deprecated and will be removed in a future release. ",
      "Please just use `outfile` instead.", call. = FALSE, immediate. = TRUE)
      if (is.null(outfile)){
        outfile <- "pvclust.RDS"
      }
    }
  }
  if (!is.null(outfile)){
    outdir <- dirname(outfile)
    if (!dir.exists(outdir)) {
      stop("Directory does not exist: ", sQuote(outdir))
    }
    if (file.access(outdir, 2) != 0) {
      stop("Directory is not writable: ", sQuote(outdir))
    }
  }
  if (is.null(max_size)) max_size <- ncol(peak_table)
  output <- match.arg(output, c("pvclust", "clusters"), several.ok = TRUE)
  if (is.data.frame(peak_table$ref_spectra) | is.matrix(peak_table$ref_spectra)){
    spectra <- peak_table$ref_spectra
  } else {
    stop("Please attach reference spectra using the `attach_ref_spectra` ",
    "function before running `cluster_spectra`.")
  }
  rm <- which(apply(spectra, 2, sd) == 0)
  if (length(rm) > 0){
    if (verbose){
      message(paste0("Removing peaks due to bad spectra: ",
                     format_list(colnames(spectra)[rm])))
    }
    spectra <- spectra[, -rm]
  }
  if (verbose)
    message('...clustering spectra')
  result <- suppressWarnings(pvclust::pvclust(spectra, method.dist = "cor", 
                                              nboot = nboot, parallel = parallel,
                                              quiet = !verbose, ...)
  )
  if (!is.null(outfile)){
    tryCatch(saveRDS(result, outfile),
             error = function(e){
               warning("Could not save to '", outfile, "': ", e$message)
             })
  }
  picks <- pvclust::pvpick(result, alpha = (1-alpha), max.only = max.only)
  cl_size <- sapply(picks$clusters, length)
  
  ## filter clusters ##
  cl_idx <- which(cl_size > min_size & cl_size < max_size)
  if (verbose && length(cl_idx) < length(cl_size)){
    message(paste0("Removing ", length(cl_size) - length(cl_idx),
                   " under- or oversized clusters from results."))
  }
  clusters <- picks$clusters[cl_idx]
  pval <- 1 - result$edges[picks$edges[cl_idx], 'au']
  pal <- rep_len(col, length(cl_idx))
  if (plot_dend){
    plot(result, labels = FALSE, cex.pv = 0.5, 
         print.pv = 'au', print.num = FALSE)
    usr <- par("usr")
    mx <- (usr[2] - usr[1]) / length(result$hclust$order) / 3
    my <- (usr[4] - usr[3]) / 200
    for (i in seq_along(clusters)){
      pos <- match(match(clusters[[i]], result$hclust$labels),
                   result$hclust$order)
      yt <- result$hclust$height[picks$edges[cl_idx[i]]] + my
      rect(min(pos) - mx, usr[3] + my, max(pos) + mx, yt, border = pal[i])
      text(mean(range(pos)), yt, labels = i, col = pal[i], pos = 3,
           font = 2, xpd = TRUE)
    }
  }
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
                col.main = pal[i],
                main = paste0('cluster ', i, '; p = ',
                              format.pval(clusters[[i]]@pval, .015, eps = .001, 
                                          digits = 2, nsmall = 2))
        )
      })
    }
  }
  pvclust <- result
  mget(output)
}
