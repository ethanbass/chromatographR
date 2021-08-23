\name{fitpeaks}
\alias{findpeaks}
\alias{fitpeaks}
\title{Fit chromatographic peaks with a gaussian profile}
\description{
  Find chromatographic peaks, and fit peak parameters using gaussian
  peak fitting.}
\usage{
findpeaks(y, sd.max = 10)
fitpeaks(y, pos)
}
\arguments{
  \item{y}{response (numerical vector)}
  \item{smooth_type}{Type of smoothing. (Only gaussian is currently supported).}
  \item{smooth_window}{Window for smoothing.}
  \item{smooth_width}{Width for smoothing.}
  \item{slope_thresh}{Minimum threshold for peak slope.}
  \item{amp_thresh}{Minimum treshold for peak amplitude.}
  \item{pos}{locations of local maxima in vector y}
  \item{w}{Width to use for peak search.}
  \item{sd.max}{Maximum standard deviation of peak.}
  \item{fit}{Function for peak fitting. (Currently only "gaussian" is supported.)}
}
\details{
  Find peaks with function \code{findpeaks} by looking for "zero-crossings" in
  the smoothed first derivative of a signal that exceed a given slope threshold.

  Peak parameters are calculated using \code{fitpeaks}, assuming a
  normal distribution. 
}
\value{
  Function \code{findpeaks} simply returns the locations of the peak centers,
  expressed as indices.
  
  Function \code{fitpeaks} returns a matrix, whose columns contain the
  following information:
  \item{rt}{location of the maximum of the peak (x)}
  \item{sd}{width of the peak (x)}
  \item{FWHM}{full width at half maximum (x)}
  \item{height}{height of the peak (y)}
  \item{area}{peak area}

  Again, the first three elements (rt, sd and FWHM) are expressed as
  indices, so not in terms of the real retention times. The
  transformation to "real" time is done in function \code{getAllPeaks}.
}
\note{The \code{findpeaks} function is adapted from matlab code in
  Prof. Tom O'Haver's \href{http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm}{Pragmatic Introduction to Signal Processing}.}
\note{The \code{fitpeaks} function is adapted from Dr. Robert Morrison's \code{\href{https://github.com/robertdouglasmorrison/DuffyTools}{DuffyTools}} package as well as code published in Ron Wehrens' \code{\link{alsace}} package.}
\author{Ethan Bass}
\seealso{\code{\link{getAllPeaks}}}
\keyword{manip}
