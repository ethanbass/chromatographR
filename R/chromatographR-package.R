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
#' using parametric time warping (ptw) (Wehrens, R., Bloemberg, T.G., and Eilers
#' P.H.C. 2015. \doi{10.1093/bioinformatics/btv299}) or variable penalty dynamic
#' time warping (VPdtw) (Clifford, D., & Stone, G. 2012. \doi{10.18637/jss.v047.i08}).
#' Peak-finding relies on the algorithm suggested by Tom O'Haver in his
#' \href{https://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm}{
#' Pragmatic Introduction to Signal Processing}. Peaks are then fitted to a 
#' gaussian or exponential-gaussian hybrid peak shape using non-linear least 
#' squares (Lan, K. & Jorgenson, J. W. 2001. \doi{10.1016/S0021-9673(01)00594-5}). 
#' More details on package usage and a suggested workflow can be found in the 
#' vignette.
#' @author Ethan Bass
"_PACKAGE"

#' Raw goldenrod root chromatograms
#' 
#' A list of four HPLC-DAD data matrices of \emph{Solidago altissima} roots
#' extracted in 90\% methanol. Retention times are stored in rows and 
#' wavelengths are stored in columns.
#' 
#' @name Sa
#' @docType data
#' @keywords data
#' @usage data(Sa)
#' @format  A list of four matrices (1301 times x 60 wavelengths).
NULL

#' Preprocessed goldenrod root chromatograms
#' 
#' A list of four pre-processed HPLC-DAD chromatograms derived from the raw data
#' stored in \code{Sa}. Retention times are stored in rows and wavelengths
#' are stored in columns. The time axis is compressed to save space and 
#' processing time so the data are a little choppy. 
#'
#' @name Sa_pr
#' @keywords data
#' @docType data
#' @usage data(Sa_pr)
#' @format  A list of four pre-processed matrices (434 retention times x 60 wavelengths).
NULL

#' Warped goldenrod root chromatograms.
#' 
#' A list of four pre-processed and warped goldenrod root chromatograms derived
#' from the raw data stored in \code{Sa}.
#'
#' @name Sa_warp
#' @docType data
#' @keywords data
#' @usage data(Sa_warp)
#' @format  A list of four pre-processed and warped matrices (434 times x 60 wavelengths).
NULL

#' Goldenrod peak table
#' 
#' Peak table generated from exemplary goldenrod root extracts.
#'
#' @name pk_tab
#' @docType data
#' @keywords data
#' @usage data(pk_tab)
#' @format  A \code{peak_table} object.
NULL
