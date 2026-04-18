#' Normalize peak table or chromatograms
#' 
#' Normalizes peak table or list of chromatograms by specified column in sample
#' metadata or in the peak table. For normalization by sample metadata, the 
#' metadata must first be attached to the \code{peak_table} using
#' \code{\link{attach_metadata}}.
#' 
#' @param peak_table A \code{peak_table} object.
#' @param column A string specifying the column containing the weights.
#' @param chrom_list List of chromatograms for normalization. The samples must
#' be in same order as the peak_table. If omitted, the function will attempt to
#' find it automatically using information stored in the \code{peak_table}.
#' @param what A \code{peak_table} or list of chromatograms (\code{chrom_list}).
#' @param by Whether to normalize by a column in sample metadata (\code{meta}) 
#' or by a column in the peak table (\code{peak}). Defaults to \code{NULL}. In
#' this case, this parameter is inferred based on the \code{column} name.
#' @param on_invalid How to handle invalid normalization values (i.e. zero,
#' negative, or \code{NA} values). One of \code{"warn"} (the default),
#' \code{"silent"}, or \code{"error"}. The former two options will replace
#' invalid values with \code{NA}. 
#' @importFrom stats setNames
#' @return A \code{peak_table} object where values are normalized as specified
#' by the specified \code{column}.
#' of each sample.
#' @author Ethan Bass
#' @seealso \code{\link{get_peaktable}} \code{\link{attach_metadata}}
#' @examples
#' data(pk_tab)
#' path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
#' meta <- read.csv(path)
#' 
#' # normalize by sample mass
#' pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
#' norm <- normalize_data(pk_tab, "mass", what = "peak_table")
#' 
#' # normalize by internal standard
#' norm <- normalize_data(pk_tab, column = "V16", by = "peak")
#' @export normalize_data

normalize_data <- function(peak_table, column, chrom_list=NULL,
                           what = c('peak_table', 'chrom_list'),
                           by = NULL, 
                           on_invalid = c("warn", "error", "silent")){
  on_invalid <- match.arg(on_invalid)
  check_peaktable(peak_table)
  if (is.null(by)){
    found_meta <- any(grepl(column, colnames(peak_table$meta)))
    found_tab <- any(grepl(column, colnames(peak_table$tab)))
    if (found_meta & found_tab){
      stop("Column could not be disambiguated. Please specify `by` argument.")
    }
    else if (found_meta){
      by <- "meta"
    } else if (found_tab){
      by <- "peak"
    }
  }
  by <- match.arg(by, c("meta", "peak"))
  if (by == "meta"){
    if (!is.data.frame(peak_table$sample_meta))
      stop("Metadata must be attached to peak table prior to normalization.")
  } 
  what <- match.arg(what, c("peak_table", "chrom_list"))
  df <- switch(by, meta = peak_table$sample_meta, peak = peak_table$tab)
  if (!(column %in% colnames(df))){
    stop(sprintf("The specified column (%s) could not be found.", 
                 sQuote(column)))
  }
  norm_vals <- check_norm_values(df[[column]], sample_names = rownames(df), 
                                 on_invalid = on_invalid)
  if (what == "peak_table"){
    peak_table$tab <- as.data.frame(
      sweep(as.matrix(peak_table$tab), 1, norm_vals, FUN = "/")
    )
    peak_table$args[c("normalized", "normalization_by", 
                      "normalization_column")] <- list(TRUE, by, column)
    return(peak_table)
  } else if (what == "chrom_list"){
    if (is.null(chrom_list)){
      chrom_list <- get_chrom_list(peak_table)
    } else {
      chrom_list <- get_chrom_list(peak_table, chrom_list)
    }
    if (!identical(names(chrom_list), rownames(peak_table$tab)))
      stop("Names of chromatograms do not match the peak table.")
    chrom_list <- setNames(lapply(seq_nrow(peak_table$tab), function(samp){
      chrom <- chrom_list[[samp]]/norm_vals[samp]
      attr(chrom, "normalized") <- TRUE
      attr(chrom, "normalized_by") <- by
      attr(chrom, "normalization_column") <- column
      chrom
    }), names(chrom_list))
    return(chrom_list)
  }
}

#' Check normalization values
#' @noRd
check_norm_values <- function(vals, sample_names,
                              on_invalid = c("warn", "error", "silent")) {
  on_invalid <- match.arg(on_invalid)
  
  invalid <- is.na(vals) | vals == 0 | vals < 0
  if (!any(invalid)) return(vals)
  
  reasons <- rep("negative", sum(invalid))
  reasons[is.na(vals[invalid])]  <- "NA"
  reasons[vals[invalid] == 0]    <- "0"
  
  msg <- sprintf("Invalid normalization values (%s)",
    paste(sQuote(sample_names[invalid]), reasons, sep = ": ", collapse = "; ")
  )
  switch(on_invalid,
         error  = stop(msg),
         warn   = {warning(msg); vals[invalid] <- NA},
         silent = {vals[invalid] <- NA}
  )
  vals
}