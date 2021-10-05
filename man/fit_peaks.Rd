\name{fit_peaks}
\alias{find_peaks}
\alias{fit_peaks}
\title{Fit chromatographic peaks with a gaussian or exponential-gaussian hybrid profile}
\description{
  Find chromatographic peaks, and fit peak parameters using gaussian or exponential-gaussian hybrid peak fitting.}
\usage{
find_peaks(y, smooth_type='gaussian', smooth_window = 1, smooth_width = 0.1,
                               slope_thresh=.05, amp_thresh=0, bounds=F)
fit_peaks(y, pos, sd.max = 50, fit = c("egh", "gaussian", "emg"), 
                      max.iter = 100)
}
\arguments{
  \item{y}{response (numerical vector)}
  \item{smooth_type}{Type of smoothing. (Only "gaussian" is currently supported).}
  \item{smooth_window}{Window for smoothing.}
  \item{smooth_width}{Width for smoothing.}
  \item{slope_thresh}{Minimum threshold for peak slope.}
  \item{amp_thresh}{Minimum treshold for peak amplitude.}
  \item{pos}{locations of local maxima in vector y}
  \item{sd.max}{Maximum standard deviation of peak.}
  \item{fit}{Function for peak fitting. (Currently "gaussian" and 
  exponential gaussian hybrid ("egh") are supported.)}
  \item{bounds}{Logical. If TRUE, includes boundaries of peak in dataframe.}
  \item{max.iter}{Maximum number of iterations to use in nonlinear least squares peak-fitting.}

}
\details{
  Find peaks with function \code{find_peaks} by looking for "zero-crossings" in
  the smoothed first derivative of a signal that exceed a given slope threshold.

  Peak parameters are calculated using \code{fit_peaks}, which fits the data to a gaussian or exponential-gaussian hybrid curve using non-linear least squares estimation as implemented in \code{\link{stats::nls}}. 
}
\value{
  Function \code{find_peaks} simply returns the locations of the peak centers,
  expressed as indices.
  
  Function \code{fit_peaks} returns a matrix, whose columns contain the
  following information:
  \item{rt}{location of the maximum of the peak (x)}
  \item{start}{start of peak (only included in table if `bounds==T`)}
  \item{end}{end of peak (only included in table if `bounds==T`)}
  \item{sd}{width of the peak (x)}
  \item{tau}{tau parameter (only included in table if `fit=="egh"`)}
  \item{FWHM}{full width at half maximum (x)}
  \item{height}{height of the peak (y)}
  \item{area}{peak area}
  \item{r.squared}{r-squared value for fit of model to data.}

  Again, the first three elements (rt, sd and FWHM) are expressed as
  indices, so not in terms of the real retention times. The
  transformation to "real" time is done in function \code{get_peaks}.
}
\note{The \code{find_peaks} function is adapted from matlab code in
  Prof. Tom O'Haver's \href{http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm}{Pragmatic Introduction to Signal Processing}.}
\note{The \code{fit_peaks} function is adapted from Dr. Robert Morrison's \href{https://github.com/robertdouglasmorrison/DuffyTools}{DuffyTools} package as well as code published in Ron Wehrens' \code{\link{alsace}} package.}
\author{Ethan Bass}
\references{
Lan, K. & Jorgenson, J. W. A hybrid of exponential and gaussian functions as a simple model of asymmetric chromatographic peaks. Journal of Chromatography A 915, 1-13 (2001).

Naish, P. J. & Hartwell, S. Exponentially Modified Gaussian functions - A good model for chromatographic peaks in isocratic HPLC? Chromatographia 26, 285-296 (1988).
}
\seealso{\code{\link{get_peaks}}}
\keyword{manip}
