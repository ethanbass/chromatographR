#' Return first or last parts of a \code{peak_table}.
#'
#' Returns the first or last parts of the \code{\link{peak_table}}.
#'
#' @importFrom utils head
#' @param x A \code{\link{peak_table}} object.
#' @param ... Additional arguments to \code{\link{head}}.
#' @method head peak_table
#' @return The first or last \code{n} rows of \code{peak_table$tab}.
#' @rdname head.peak_table
#' @keywords internal
#' @export
head.peak_table <- function(x,...){
  head(x$tab, ...)
}

#' @importFrom utils tail
#' @param x A \code{\link{peak_table}} object.
#' @param ... Additional arguments to \code{\link{tail}}.
#' @rdname head.peak_table
#' @method tail peak_table
#' @keywords internal
#' @export
tail.peak_table <- function(x,...){
  tail(x$tab, ...)
}


#' Return dimensions of a \code{peak_table} object.
#'
#' Returns the dimensions of a \code{peak_table}, where the first dimension is
#' the number of samples and the second dimension is the number of peaks.
#' 
#' @param x A \code{\link{peak_table}} object.
#' @method dim peak_table
#' @keywords internal
#' @return Returns the number of rows and columns in \code{peak_table$tab}.
#' @export
dim.peak_table <- function(x){
  dim(x$tab)
}

#' Return row names from a \code{peak_table} object.
#'
#' These will be the names of the samples.
#'
#' @param x A \code{\link{peak_table}} object.
#' @param ... Additional arguments to \code{\link{row.names}}.
#' @method row.names peak_table
#' @keywords internal
#' @return Returns the row names of \code{peak_table$tab}.
#' @export
row.names.peak_table <- function(x, ...){
  row.names(x$tab, ...)
}

#' Subset peak table
#' 
#' Returns subset of \code{peak_table} object.
#' 
#' @param x A \code{peak_table} object.
#' @param subset Logical expression indicating rows (samples) to keep from
#' \code{peak_table}; missing values are taken as false.
#' @param select Logical expression indicating columns (peaks) to select from
#' \code{peak_table}.
#' @param drop Logical. Passed to indexing operator.
#' @param ... Additional arguments (placeholder).
#' @return A \code{peak_table} object with samples specified by \code{subset}
#' and peaks specified by \code{select}.
#' @author Ethan Bass
#' @method subset peak_table
#' @keywords internal
#' @export
subset.peak_table <- function(x, subset, select, drop = FALSE, ...){
  x$tab <- subset(x$tab, subset = subset, 
                  select = select, drop = drop)
  if (!is.null(dim(x$ref_spectra))){
    x$sample_meta <- subset(x$sample_meta, subset = subset, drop = drop)
  }
  if (!missing(select)){
    x$pk_meta <- subset(x$pk_meta, select = select, drop = drop)
    if (!is.null(dim(x$ref_spectra))){
      x$ref_spectra <- subset(x$ref_spectra, select = select, drop = drop)
    }
  }
  x
}

#' Summarize Peak Table
#' 
#' Prints basic statistics about \code{peak_table} to the console (e.g., the 
#' number of peaks, wavelengths, range of intensities, etc).
#' 
#' @param object A \code{\link{peak_table}} object.
#' @param ... Additional arguments (currently unused).
#' @author Ethan Bass
#' @method summary peak_table
#' @keywords internal
#' @export
#' @return A list containing basic information about the supplied 
#' \code{peak_table} object.
summary.peak_table <- function(object, ...) {
  
  # Calculate summary statistics
  intensities <- as.vector(as.matrix(object$tab))
  intensities <- intensities[!is.na(intensities)]
  
  summary_list <- list(
    dimensions = list(
      samples = nrow(object$tab),
      peaks = ncol(object$tab)
    ),
    retention_time = if (!is.null(object$pk_meta)) {
      range(object$pk_meta["rt",], na.rm = TRUE)
    } else NULL,
    wavelengths = object$args$lambdas[[1]],
    intensity = list(
      range = range(intensities),
      mean = mean(intensities),
      median = median(intensities),
      sd = sd(intensities)
    ),
    time_units = object$args$time_unit,
    response = object$args$response,
    has_sample_meta = !is.null(object$sample_meta),
    has_ref_spectra = !is.null(object$ref_spectra),
    normalized = object$args$normalized,
    normalization_method = if(object$args$normalized) object$args$normalization_by else NA
  )
  
  class(summary_list) <- c("summary.peak_table", "list")
  return(summary_list)
}

#' Peak table summary
#' @noRd
print.summary.peak_table <- function(x, ...) {
  cat("Peak Table Summary\n")
  cat("==================\n\n")
  
  cat("Dimensions:\n")
  cat("  Samples:", x$dimensions$samples, "\n")
  cat("  Peaks:  ", x$dimensions$peaks, "\n\n")
  
  if (!is.null(x$retention_time)) {
    cat("Retention Time:\n")
    cat("  Range:", paste0("[", round(x$retention_time[1], 2), ", ", 
                           round(x$retention_time[2], 2), "] min\n"))
  }
  cat("  Wavelengths:  ", paste(x$wavelengths, collapse = ", "),
      "\n")
  
  cat("\nIntensity:\n")
  cat("  Range: ", paste0("[", round(x$intensity$range[1], 2), ", ", 
                          round(x$intensity$range[2], 2), "]\n"))
  cat("  Mean:  ", round(x$intensity$mean, 2), "\n")
  cat("  Median:", round(x$intensity$median, 2), "\n")
  
  cat("\nArguments:\n")
  cat("  Response:", tools::toTitleCase(x$response), "\n")
  cat("  Sample Metadata:", ifelse(x$has_sample_meta, "Yes", "No"), "\n")
  cat("  Reference Spectra:", ifelse(x$has_ref_spectra, "Yes", "No"), "\n")
  cat("  Normalization:", ifelse(x$normalized, 
                                 paste("Normalized by", x$normalization_method), 
                                 "N/A"), "\n")
  
  invisible(x)
}
