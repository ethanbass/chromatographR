#' Attach experimental metadata
#' 
#' Attaches experimental meta-data to `peak_table` object. One of the columns in
#' the supplied meta-data must match exactly the row names of the peak table.
#' 
#' @aliases attach_metadata
#' @param peak_table A `peak_table` object.
#' @param metadata A `data.frame` containing the sample meta-data.
#' @param column The name of the column containing the sample names.
#' @return A \code{peak_table} object with attached meta-data.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaktable}} \code{\link{normalize_data}}
#' @examples
#' data(pk_tab)
#' path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
#' meta <- read.csv(path)
#' pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
#' @export attach_metadata

attach_metadata <- function(peak_table, metadata, column){
  if (!(column %in% colnames(metadata)))
    stop(paste0("Column, ", column, ", is not found."))
  if (sum((duplicated(metadata[,column]))) > 0)
    stop(paste0("Sample names must be unique. Please check column '", column,
    "' for duplicates."))
  meta <- data.frame(rownames(peak_table$tab))
  names(meta) <- column
  metadata[,column] <- as.character(metadata[,column])
  missing_meta <- !(meta[,column] %in% metadata[,column])
  if (sum(missing_meta)>0)
    warning("The supplied metadata does not include all samples.")
  meta <- merge(meta, metadata, all.x=TRUE, all.y=FALSE, sort=FALSE)
  peak_table$sample_meta <- meta
  return(peak_table)
}

#' Gather reference spectra.
#' 
#' Defines reference spectra. Reference spectra are defined either as the spectrum
#' with the highest intensity (`ref == "max.int"`) or as the spectrum with the 
#' highest average correlation to the rest of the spectra in the peak_table
#' (`ref == "max.cor"`).
#' 
#' @importFrom stats cor sd
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param ref What criterion to use to select reference spectra.
#' Current options are maximum correlation (`max.cor`) or maximum signal
#' intensity (`max.int`).
#' @return A matrix consisting of reference spectra for each peak in the
#' provided peak table.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
#' @examples \dontrun{
#' ref_m <- gather_reference_spectra(pk_tab, ref = "max.int")
#' ref_c <- gather_reference_spectra(pk_tab, ref="max.cor")
#' }
#' @noRd

gather_reference_spectra <- function(peak_table, chrom_list, ref = c("max.cor","max.int")){
  if (missing(chrom_list)){
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
      plot_spectrum(loc = pk, peak_table, chrom_list, plot_trace=FALSE,
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
#' 
#' Gathers reference spectra and attaches them to peak_table object. Reference 
#' spectra are defined either as the spectrum with the highest intensity (
#' \code{max.int}) or as the spectrum with the  highest average correlation
#' to the other spectra in the peak_table (\code{max.cor}).
#' 
#' @aliases attach_ref_spectra
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param chrom_list A list of chromatograms in matrix form (timepoints x
#' wavelengths).
#' @param ref What criterion to use to select reference spectra.
#' Current options are maximum correlation (\code{max.cor}) or maximum signal
#' intensity (\code{max.int}).
#' @return A `peak_table` object with reference spectra attached
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}} \code{\link{get_peaktable}}
#' @examples
#' data(pk_tab)
#' pk_tab <- attach_ref_spectra(pk_tab, ref="max.int")
#' pk_tab <- attach_ref_spectra(pk_tab, ref = "max.cor")
#' @export attach_ref_spectra
#' 
attach_ref_spectra <- function(peak_table, chrom_list, ref = c("max.cor","max.int")){
  # if (missing(chrom_list)){
  #   chrom_list <- get_chrom_list(peak_table)
  # }
  peak_table$ref_spectra <- gather_reference_spectra(peak_table, chrom_list, ref)
  return(peak_table)
}

#' Normalize peak table or chromatograms
#' 
#' Normalizes peak table or list of chromatograms by specified column in sample
#' meta-data. Metadata must first be attached to \code{peak_table} using
#' \code{\link{attach_metadata}}.
#' 
#' @param peak_table A `peak_table` object
#' @param column The name of the column containing the weights.
#' @param chrom_list List of chromatograms for normalization. The samples must
#' be in same order as the peak_table.
#' @param what `peak_table` or list of chromatograms (`chrom_list`).
#' @return A normalized \code{peak_table} object.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaktable}} \code{\link{attach_metadata}}
#' @examples
#' data(pk_tab)
#' path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
#' meta <- read.csv(path)
#' pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
#' norm <- normalize_data(pk_tab, "mass", what = "peak_table")
#' @export normalize_data

normalize_data <- function(peak_table, column, chrom_list,
                           what=c('peak_table','chrom_list')){
  if (!is.data.frame(peak_table$sample_meta))
    stop("Meta-data must be attached to peak_table prior to normalization.")
  if (!(column %in% colnames(peak_table$sample_meta)))
    stop(paste0("Column, ", column, ", is not found."))
  what <- match.arg(what, c("peak_table", "chrom_list"))
  if (what == "peak_table"){
    peak_table$tab <- t(sapply(seq_len(nrow(peak_table$tab)), function(samp){
      peak_table$tab[samp,]/peak_table$sample_meta[samp,column]
    }))
    return(peak_table)
  } else if (what == "chrom_list"){
    if (missing(chrom_list)){
      chrom_list <- try(get(peak_table$args["chrom_list"]))
      if (inherits(chrom_list, "try-error")) stop("Chromatograms not found!")
    }
    if (mean(elementwise.all.equal(names(chrom_list),rownames(peak_table$tab))) < 1)
      stop("Names of chromatograms do not match the peak table.")
    chrom_list <- lapply(seq_len(nrow(peak_table$tab)), function(samp){
      chrom_list[[samp]]/peak_table$sample_meta[samp,column]
    })
  }
}