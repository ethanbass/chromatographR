#' Attach experimental metadata
#' 
#' Attaches sample metadata to `peak_table` object. Metadata should be
#' provided as a data.frame object. One of the columns in
#' the supplied metadata must match exactly the row names of the peak table.
#' 
#' @aliases attach_metadata
#' @param peak_table A `peak_table` object.
#' @param metadata A `data.frame` containing the sample metadata.
#' @param column The name of the column in your \code{metadata} object containing the
#' sample names. Sample names must match the row names of \code{peak_table$tab}.
#' @return A \code{peak_table} object with attached metadata in the \code{
#' $sample_meta} slot.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaktable}} \code{\link{normalize_data}}
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
    stop(paste0("Column, ", column, ", is not found."))
  if (sum((duplicated(metadata[,column], incomparables = NA))) > 0)
    stop(paste("Sample names must be unique. Please check column", sQuote(column),
    "for duplicates."))
  if (!inherits(peak_table,"peak_table"))
    stop(paste("Provided peak table object must be of class 'peak_table'."))
  meta <- data.frame(rownames(peak_table$tab))
  names(meta) <- column
  metadata[, column] <- as.character(metadata[, column])
  missing_meta <- !(meta[, column] %in% metadata[, column])
  if (sum(missing_meta) > 0)
    warning("The supplied metadata does not include all samples.")
  meta <- keep_order(meta, merge, y=metadata, all.x = TRUE, all.y = FALSE,
                     sort = FALSE, by = column)
  peak_table$sample_meta <- meta
  return(peak_table)
}

#' note: convenience function from stackoverflow:
#' https://stackoverflow.com/questions/17878048/merge-two-data-frames-while-keeping-the-original-row-order
#' @noRd
keep_order <- function(data, fn, ...) { 
  col <- ".sortColumn"
  data[,col] <- 1:nrow(data) 
  out <- fn(data, ...) 
  if (!col %in% colnames(out)) stop("Ordering column not preserved by function") 
  out <- out[order(out[,col]),] 
  out[,col] <- NULL 
  out 
} 


#' Get reference spectra.
#' 
#' Defines reference spectra. Reference spectra are defined either as the 
#' spectrum with the highest intensity (\code{max.int}) or as the spectrum 
#' with the highest average correlation to the rest of the spectra in the
#' \code{peak_table} (\code{max.cor}).
#' 
#' @importFrom stats cor sd
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param chrom_list A list of chromatograms in matrix format (timepoints x
#' wavelengths). If no argument is provided here, the function will try to find
#' the \code{chrom_list} object used to create the provided \code{peak_table}.
#' @param ref What criterion to use to select reference spectra.
#' Current options are maximum correlation (\code{max.cor}) or maximum signal
#' intensity (\code{max.int}).
#' @return A matrix consisting of reference spectra for each peak in the
#' provided peak table.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}}
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
      x <- plot_all_spectra(peak = pk, peak_table, chrom_list,
                       plot_spectrum = FALSE, export_spectrum = TRUE,
                       scale_spectrum = TRUE)
      apply(x, 2, as.numeric)
    })
    sp.ref <- sapply(seq_along(sp.l), function(i){
      sp.l[[i]][,which.max(colMeans(cor(sp.l[[i]][,which(apply((sp.l[[i]]), 2, sd)!=0), drop=FALSE])))]
    })
  } else {
    sp.ref <- sapply(colnames(peak_table$tab), function(pk){
      try(plot_spectrum(loc = pk, peak_table, chrom_list, plot_trace=FALSE,
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

#' Attach reference spectra
#' 
#' Gathers reference spectra and attaches them to peak_table object. Reference 
#' spectra are defined either as the spectrum with the highest intensity (
#' \code{max.int}) or as the spectrum with the highest average correlation
#' to the other spectra in the peak_table (\code{max.cor}).
#' 
#' @aliases attach_ref_spectra
#' @param peak_table Peak table from \code{\link{get_peaktable}}.
#' @param chrom_list A list of chromatograms in matrix format (timepoints x
#' wavelengths). If no argument is provided here, the function will try to find
#' the \code{chrom_list} object used to create the provided \code{peak_table}.
#' @param ref What criterion to use to select reference spectra.
#' Current options are maximum correlation (\code{max.cor}) or maximum signal
#' intensity (\code{max.int}).
#' @return A \code{peak_table} object with reference spectra attached in the
#' \code{$ref_spectra} slot.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaks}} \code{\link{get_peaktable}}
#' @examples
#' data(pk_tab)
#' pk_tab <- attach_ref_spectra(pk_tab, ref="max.int")
#' pk_tab <- attach_ref_spectra(pk_tab, ref = "max.cor")
#' @export attach_ref_spectra

attach_ref_spectra <- function(peak_table, chrom_list, ref = c("max.cor","max.int")){
  check_peaktable(peak_table)
  peak_table$ref_spectra <- get_reference_spectra(peak_table, chrom_list, ref = ref)
  peak_table$args["reference_spectra"] <- ref
  return(peak_table)
}

#' Normalize peak table or chromatograms
#' 
#' Normalizes peak table or list of chromatograms by specified column in sample
#' metadata. Metadata must first be attached to \code{peak_table} using
#' \code{\link{attach_metadata}}.
#' 
#' @param peak_table A `peak_table` object
#' @param column The name of the column containing the weights.
#' @param chrom_list List of chromatograms for normalization. The samples must
#' be in same order as the peak_table. If no argument is provided here, the
#' function will try to find the \code{chrom_list} object used to create the
#' provided \code{peak_table}.
#' @param what `peak_table` or list of chromatograms (`chrom_list`).
#' @param by Whether to normalize by a column in sample metadata (\code{meta}) or
#' by a column in the peak table itself (\code{peak}).
#' @return A \code{peak_table} object where the peaks are normalized by the mass
#' of each sample.
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
                           what = c('peak_table','chrom_list'),
                           by=c("meta", "peak")){
  check_peaktable(peak_table)
  if (!is.data.frame(peak_table$sample_meta))
    stop("Metadata must be attached to peak_table prior to normalization.")
  if (!(column %in% colnames(peak_table$sample_meta)))
    stop(paste0("The specified column (", sQuote(column), ") could not be found."))
  what <- match.arg(what, c("peak_table", "chrom_list"))
  by <- match.arg(by, c("meta", "peak"))
  df <- switch(by, meta = peak_table$sample_meta, peak = peak_table$tab)
  if (what == "peak_table"){
    pktab <- as.data.frame(t(sapply(seq_len(nrow(peak_table$tab)), function(samp){
      as.numeric(as.vector(peak_table$tab[samp,]))/df[samp,column]
    })))
    rownames(pktab) <- rownames(peak_table$tab)
    peak_table$tab <- pktab
    peak_table$args[c("normalized", "normalization_by")] <- c(TRUE, column)
    return(peak_table)
  } else if (what == "chrom_list"){
    if (missing(chrom_list)){
      chrom_list <- get_chrom_list(peak_table)
    } else get_chrom_list(peak_table, chrom_list)
    if (mean(elementwise.all.equal(names(chrom_list), rownames(peak_table$tab))) < 1)
      stop("Names of chromatograms do not match the peak table.")
    chrom_list <- lapply(seq_len(nrow(peak_table$tab)), function(samp){
      chrom_list[[samp]]/df[samp,column]
    })
    return(chrom_list)
  }
}
