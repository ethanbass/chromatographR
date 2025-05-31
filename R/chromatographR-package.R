#' chromatographR
#' 
#' Chromatographic Data Analysis Toolset
#' 
#' Tools for high-throughput analysis of HPLC-DAD/UV
#' chromatograms (or similar data). Includes functions for preprocessing, alignment,
#' peak-finding and fitting, peak-table construction, data-visualization, etc.
#' Preprocessing and peak-table construction follow the rough formula laid out
#' in alsace (Wehrens, R., Bloemberg, T.G., and Eilers P.H.C., 2015.
#' \doi{10.1093/bioinformatics/btv299}). Alignment of chromatograms is available
#' using parametric time warping (PTW) (Wehrens, R., Bloemberg, T.G., and Eilers
#' P.H.C. 2015. \doi{10.1093/bioinformatics/btv299}) or variable penalty dynamic
#' time warping (VPdtw) (Clifford, D., & Stone, G. 2012. \doi{10.18637/jss.v047.i08}).
#' Peak-finding relies on the algorithm suggested by Tom O'Haver in his
#' \href{https://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm}{
#' Pragmatic Introduction to Signal Processing}. Peaks are then fitted to a 
#' gaussian or exponential-gaussian hybrid peak shape using non-linear least 
#' squares (Lan, K. & Jorgenson, J. W. 2001. \doi{10.1016/S0021-9673(01)00594-5}). 
#' More details on package usage and a suggested workflow can be found in the 
#' vignette.
#' 
#' @section Analysis functions:
#' \itemize{
#'   * [read_chroms()]: Import chromatograms from a variety of vendor formats.
#'   * [preprocess()]: Preprocess chromatographic matrices.
#'   * [correct_rt()]: Align chromatograms.
#'   * [get_peaks()]: Find and integrate peaks.
#'   * [get_peaktable()]: Assemble peak table.
#'   * [attach_metadata()]: Attach metadata to peak table.
#'   * [attach_ref_spectra()]: Attach reference spectra to peak table.
#'   * [normalize_data()]: Normalize \code{peaktable} or \code{chrom_list}.
#' }
#' 
#' @section Visualization functions:
#' \itemize{
#'   * [boxplot.peak_table()]: Create boxplot from peaktable object.
#'   * [plot_chroms()]: Plot chromatograms as traces.
#'   * [plot_chroms_heatmap()]: Plot chromatograms as heatmap.
#'   * [plot_spectrum()]: Plot spectrum and/or trace of specified peak.
#'   * [plot_all_spectra()]: Plot all spectra for specified peak.
#'   * [mirror_plot()]: Plot chromatograms as mirror plot.
#'   * [scan_chrom()]: Plot spectrum at wavelength specified by clicking on a chromatogram.
#'   * [plot.peak_list()]: Plot fitted peaks over chromatographic trace.
#'   * [plot.peak_table()]: Plot the trace and/or spectrum of a specified peak from the peak table.
#'   * [plot.ptw_list()]: Plot PTW alignment object.
#' }
#' 
#' @section Utility functions:
#' \itemize{
#' * [combine_peaks()]: Combine duplicate peaks in peak table based on retention time and spectral similarity.
#' * [merge_peaks()]: Merge split peaks into a single column of a peak table.
#' * [get_times()]: Return retention times from a peak table or a list of chromatograms.
#' * [get_lambdas()]: Return wavelengths from a peak table or a list of chromatograms.
#' * [reshape_chroms()]: Reshape a list of chromatograms to long format.
#' * [reshape_peaktable()]: Reshape a \code{peak_table} object to long format.
#' * [write_peaktable()]: Export peak table in \code{csv} or \code{xlsx} format.
#' }
#' @section Example data:
#' \itemize{
#'   * [Sa]: A list of four goldenrod root chromatograms.
#'   * [Sa_pr]: Preprocessed goldenrod root chromatograms.
#'   * [Sa_warp]: Preprocessed and aligned goldenrod root chromatograms.
#'   * [pk_tab]: Peak table from aligned goldenrod root chromatograms.
#' }
#' @author Ethan Bass
#' @md
"_PACKAGE"

#' Raw goldenrod root chromatograms
#' 
#' A list of four HPLC-DAD data matrices of \emph{Solidago altissima} roots
#' extracted in 90\% methanol. Retention times are stored in rows and 
#' wavelengths are stored in columns. Data were collected on a Agilent 1100 HPLC.
#' 
#' @name Sa
#' @docType data
#' @keywords data
#' @usage data(Sa)
#' @format  A list of four matrices (1301 times x 60 wavelengths).
#' @family data objects
NULL

#' Preprocessed goldenrod root chromatograms
#' 
#' A list of four pre-processed HPLC-DAD chromatograms derived from the raw data
#' stored in \code{\link{Sa}}. Retention times are stored in rows and wavelengths
#' are stored in columns. The time axis is compressed to save space and 
#' processing time so the data are a little choppy. 
#'
#' @name Sa_pr
#' @keywords data
#' @docType data
#' @usage data(Sa_pr)
#' @format  A list of four pre-processed matrices (434 retention times x 60 wavelengths).
#' @family data objects
NULL

#' Warped goldenrod root chromatograms.
#' 
#' A list of four pre-processed and warped goldenrod root chromatograms derived
#' from the raw data stored in \code{\link{Sa}}.
#'
#' @name Sa_warp
#' @docType data
#' @keywords data
#' @usage data(Sa_warp)
#' @format  A list of four pre-processed and warped matrices (434 times x 60 wavelengths).
#' @family data objects
NULL

#' Goldenrod peak table
#' 
#' A peak table generated from the exemplary goldenrod root extracts stored in
#' \code{\link{Sa}}.
#'
#' @name pk_tab
#' @docType data
#' @keywords data
#' @usage data(pk_tab)
#' @format  A \code{peak_table} object.
#' @family data objects
NULL
