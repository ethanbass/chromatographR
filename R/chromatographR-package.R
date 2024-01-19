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
#' Four HPLC-DAD data matrices of *Solidago altissima* roots extracted in 90%
#' methanol. The time axis is compressed to save space so the
#' data are a little choppy.
#' 
#' @name Sa
#' @docType data
#' @format  A list of four matrices (time x wavelength).

NULL

#' Preprocessed goldenrod root chromatograms
#' 
#' These are the pre-processed chromatograms derived from the raw data stored in
#' \code{Sa}.
#'
#' @name Sa_pr
#' @docType data
#' @format  Four pre-processed matrices (time x wavelength) to use in examples.
NULL

#' Warped goldenrod root chromatograms.
#' 
#' These are pre-processed and warped goldenrod root chromatograms derived from
#' the raw data stored in \code{Sa}.
#'
#' @name Sa_warp
#' @docType data
#' @format  Four pre-processed and warped matrices (time x wavelength) to use in
#' examples.
NULL

#' Goldenrod peak table
#' 
#' Peak table generated from exemplary goldenrod root extracts for use in examples.
#'
#' @name pk_tab
#' @docType data
#' @format  A \code{peak_table} object.
NULL


