#' Attach experimental metadata
#' 
#' Attach experimental meta-data to peak table.
#' @aliases attach_metadata
#' @param peak_table Peak table
#' @param metadata Dataframe
#' @param column Name of the column with sample names.
#' @return Peak table with meta-data.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaktable}}
#' @examples \dontrun{
#' meta <- read.csv("meta.csv")
#' pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta,
#' column="vial")}
#' @export attach_metadata

attach_metadata <- function(peak_table, metadata, column){
  if (!(column %in% colnames(metadata)))
    stop(paste0("Column, ", column, ", is not found."))
  meta <- data.frame(rownames(peak_table$tab))
  names(meta) <- column
  metadata[,column] <- as.character(metadata[,column])
  meta <- merge(meta, metadata, all.x=TRUE, all.y=FALSE, sort=FALSE)
  peak_table$sample_meta <- meta
  return(peak_table)
}

#' Function to gather reference spectra.
#' 
#' @importFrom stats cor sd
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param ref What criterion to use to select reference spectra.
#' Current options are maximum correlation ("max.cor") or maximum signal
#' intensity ("max.int").
#' @return A matrix consisting of reference spectra for each peak in the
#' provided peak table.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @examples 
#' ref_m <- gather_reference_spectra(pk_tab, ref = "max.int")
#' ref_c <- gather_reference_spectra(pk_tab, ref="max.cor")
#' @noRd

gather_reference_spectra <- function(peak_table, chrom_list=NULL, ref = c("max.cor","max.int")){
  if (is.null(chrom_list)){
    try.out <- try(chrom_list <- get(peak_table$args["chrom_list"]))
    if (inherits(try.out, "try-error")) stop("Chromatograms not found!")
  }
  ref <- match.arg(ref, c("max.cor", "max.int"))
  X<-colnames(peak_table$tab)
  if (ref=="max.cor"){
    sp.l <- lapply(X,function(pk){
      plot_all_spectra(peak = pk, peak_table, chrom_list,
                       plot_spectrum = FALSE, export_spectrum = TRUE,
                       scale_spectrum = TRUE)
    })
    sp.ref <- sapply(seq_along(sp.l), function(i){
      sp.l[[i]][,which.max(colMeans(cor(sp.l[[i]][,which(apply((sp.l[[i]]),2,sd)!=0)])))]})
    # sp.ref <- data.frame(do.call(cbind,sp.ref))
  } else {
    sp.ref <- sapply(colnames(peak_table$tab), function(pk){
      plot_spectrum(loc = pk, peak_table, plot_trace=FALSE,
                    plot_spectrum = FALSE, export_spectrum = TRUE,
                    verbose = FALSE,
                    scale_spectrum = TRUE)})
    sp.ref <- do.call(cbind, sp.ref)
  }
  colnames(sp.ref) <- colnames(peak_table$tab)
  rownames(sp.ref) <- colnames(chrom_list[[1]])
  return(sp.ref)
}

#' Attach reference spectra
#' @aliases attach_ref_spectra
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param ref What criterion to use to select reference spectra.
#' Current options are maximum correlation ("max.cor") or maximum signal
#' intensity ("max.int").
#' @return peak_table object with reference spectra attached
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @examples \dontrun{
#' ref_m <- attach_ref_spectra(pk_tab, ref = "max.int")}
#' @seealso \code{\link{get_peaktable}}
#' @export attach_ref_spectra
attach_ref_spectra <- function(peak_table, chrom_list=NULL, ref = c("max.cor","max.int")){
  peak_table$ref_spectra <- gather_reference_spectra(peak_table, chrom_list, ref)
  return(peak_table)
}
