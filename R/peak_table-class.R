#' Peak Table Object
#'
#' @name peak_table-class  
#' @aliases peak_table-class peak_table
#'
#' @description
#' S3 objects of class `peak_table` store chromatographic peak data along with
#' associated metadata and reference spectra.
#'
#' @section Components:
#' A `peak_table` object contains the following elements:
#' * `tab`: The peak table itself -- a `data.frame` of intensities in a
#' sample x peak configuration.
#' * `pk_meta`: A `data.frame` containing peak meta-data (e.g., the spectral
#' component, peak number, and average retention time).
#' * `sample_meta`: A `data.frame` of sample meta-data. Must be added using
#' [`attach_metadata`].
#' * `ref_spectra`: A `data.frame` of reference spectra (in a wavelength x peak
#' configuration). Must be added using [`attach_ref_spectra`].
#' * `args`: A vector of arguments given to [`get_peaktable`] to generate
#' the peak table.
#'
#' @section Methods:
#' The following methods are available:
#' - [print.peak_table()]
#' - [dim.peak_table()]
#' - [head.peak_table()]  
#' - [tail.peak_table()]
#' - [row.names.peak_table()]
#' - [subset.peak_table()]
#' - [summary.peak_table()]
#' - \code{\link{[.peak_table}}
#'
#' @seealso [get_peaktable()] for creating peak_table objects
#'
#' @examples
#' # See get_peaktable() for creation examples
#' # To access peak_table components:
#' data(pk_tab)
#' pk_tab$tab          # peak table
#' pk_tab$pk_meta      # peak metadata  
#' pk_tab$sample_meta  # sample metadata
#'
#' @md
NULL