#' Attaches sample metadata to a `peak_table` object by matching sample names.
#' 
#' Metadata is provided as a `data.frame`, with one column containing sample
#' identifiers matching the row names of `peak_table$tab`.
#' 
#' @aliases attach_metadata
#' @param peak_table A `peak_table` object created by [`get_peaktable`].
#' @param metadata A `data.frame` of sample metadata.
#' @param column Name of the column in `metadata` containing sample identifiers.
#' Must match the row names of `peak_table$tab`.
#' @return A `peak_table` object with metadata stored in the `sample_meta` slot.
#' @author Ethan Bass
#' @seealso [`get_peaktable`], [`normalize_data`]
#' @examples
#' data(pk_tab)
#' path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
#' meta <- read.csv(path)
#' pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
#' @export attach_metadata

attach_metadata <- function(peak_table, metadata, column){
  check_peaktable(peak_table)
  if (any(grepl("tbl", class(metadata)))){
    metadata <- as.data.frame(metadata)
  }
  if (!inherits(metadata, "data.frame")){
    stop("Please provide metadata as a `data.frame`")
  }
  if (!(column %in% colnames(metadata)))
    stop(sprintf("Column %s could not be found.", sQuote(column)))
  if (sum((duplicated(metadata[,column], incomparables = NA))) > 0)
    stop(sprintf("Sample names must be unique. Please check column %s for duplicates.",
                 sQuote(column)))
  if (!inherits(peak_table,"peak_table"))
    stop(paste("Provided peak table object must be of the 'peak_table' class."))
  meta <- data.frame(rownames(peak_table$tab))
  names(meta) <- column
  metadata[, column] <- as.character(metadata[, column])
  missing_meta <- !(meta[, column] %in% metadata[, column])
  if (sum(missing_meta) > 0)
    warning("The supplied metadata does not include all samples.")
  meta <- keep_order(meta, merge, y = metadata, all.x = TRUE, all.y = FALSE,
                     sort = FALSE, by = column)
  peak_table$sample_meta <- meta
  return(peak_table)
}

#' Attach reference spectra to a `peak_table` object.
#' 
#' Reference spectra are selected either as the spectrum with maximum intensity
#' (`"max.int"`) or as the spectrum with the highest average correlation to all
#' other spectra in the peak table (`"max.cor"`).
#' 
#' @aliases attach_ref_spectra
#' @param peak_table A `peak_table` object created by [`get_peaktable`].
#' @param chrom_list Optional list of chromatograms (timepoints × wavelengths).
#' If `NULL`, the function attempts to retrieve the original `chrom_list`
#' used to generate `peak_table`.
#' @param ref Method for selecting reference spectra. Options are `"max.int"`
#' (maximum intensity) or `"max.cor"` (maximum average correlation).
#' @return A `peak_table` object with reference spectra stored in the
#' `$ref_spectra` slot.
#' @author Ethan Bass
#' @seealso [`get_peaks`], [`get_peaktable`]
#' @examples
#' data(pk_tab)
#' data(Sa_warp)
#' pk_tab <- attach_ref_spectra(pk_tab, ref = "max.int")
#' pk_tab <- attach_ref_spectra(pk_tab, ref = "max.cor")
#' @export attach_ref_spectra

attach_ref_spectra <- function(peak_table, chrom_list, 
                               ref = c("max.cor", "max.int")){
  check_peaktable(peak_table)
  ref <- match.arg(ref, c("max.cor", "max.int"))
  peak_table$ref_spectra <- get_reference_spectra(peak_table, chrom_list,
                                                  ref = ref)
  peak_table$args["reference_spectra"] <- ref
  return(peak_table)
}

#' note: convenience function from stackoverflow:
#' https://stackoverflow.com/questions/17878048/merge-two-data-frames-while-keeping-the-original-row-order
#' @noRd
keep_order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- seq_nrow(data)
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 

#' Get reference spectra.
#' 
#' Defines reference spectra. Reference spectra are defined either as the 
#' spectrum with the highest intensity (`"max.int"`) or as the spectrum 
#' with the highest average correlation to the rest of the spectra in the
#' `peak_table` (`"max.cor"`).
#' 
#' @importFrom stats cor sd
#' @param peak_table Peak table from [`get_peaktable`].
#' @param chrom_list A list of chromatograms in matrix format (timepoints x
#' wavelengths). If no argument is provided here, the function will try to find
#' the `chrom_list` object used to create the provided `peak_table`.
#' @param ref What criterion to use to select reference spectra.
#' Current options are maximum correlation (`"max.cor"`) or maximum signal
#' intensity (`"max.int"`).
#' @return A matrix consisting of reference spectra for each peak in the
#' provided peak table.
#' @author Ethan Bass
#' @seealso [`get_peaks`]
#' @noRd

get_reference_spectra <- function(peak_table, chrom_list,
                                     ref = c("max.cor", "max.int")){
  check_peaktable(peak_table)
  if (!inherits(peak_table, "peak_table"))
    stop("Provided peak_table object must be a `peak_table` object.")
  if (missing(chrom_list)){
    chrom_list <- get_chrom_list(peak_table)
  } else get_chrom_list(peak_table, chrom_list)
  ref <- match.arg(ref, c("max.cor", "max.int"))
  X <- colnames(peak_table$tab)
  if (ref == "max.cor"){
    sp.l <- lapply(X,function(pk){
      x <- plot_all_spectra(loc = pk, peak_table, chrom_list,
                       plot_spectrum = FALSE, export_spectrum = TRUE,
                       scale_spectrum = TRUE)
      apply(x, 2, as.numeric)
    })
    sp.ref <- sapply(seq_along(sp.l), function(i){
      sp.l[[i]][, which.max(
        colMeans(cor(sp.l[[i]][, which(apply((sp.l[[i]]), 2, sd) != 0), 
                                                          drop = FALSE])))]
    })
  } else {
    sp.ref <- sapply(colnames(peak_table$tab), function(pk){
      try(plot_spectrum(loc = pk, peak_table, chrom_list, plot_trace = FALSE,
                    plot_spectrum = FALSE, export_spectrum = TRUE,
                    verbose = FALSE,
                    scale_spectrum = TRUE, engine = "base")
      )
    })
    sp.ref <- do.call(cbind, sp.ref)
  }
  colnames(sp.ref) <- colnames(peak_table$tab)
  rownames(sp.ref) <- colnames(chrom_list[[1]])
  return(sp.ref)
}
