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
#' * `tab`: The peak table itself — a `data.frame` of
#' intensities in a sample × peak configuration.
#' * `pk_meta`: A `data.frame` containing peak metadata (e.g., the spectral
#' component, peak number, and summarized retention time). See
#' [`get_peaktable`] for its rows.
#' * `sample_meta`: A `data.frame` of sample metadata: the fields added to the
#' chromatograms with `chromConverter::add_metadata()` (chromConverter 0.11.0
#' or later), and any attached with [`attach_metadata`]. `NA` if there are none.
#' * `instrument_meta`: A `data.frame` of instrumental metadata transferred
#' automatically from the chromatograms.
#' * `ref_spectra`: A matrix of reference spectra (in a wavelength × peak
#' configuration). `NA` until added with [`attach_ref_spectra`].
#' * `args`: A list of the arguments given to [`get_peaktable`] to generate
#' the peak table, with the peak list's fitting settings. See
#' [`get_peaktable`].
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
#' @seealso [get_peaktable()] for creating `peak_table` objects
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