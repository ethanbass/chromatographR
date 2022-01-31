\name{fit_peaks}
\alias{fit_peaks}
\title{Fit chromatographic peaks with a gaussian or exponential-gaussian hybrid profile}
\description{ Fit peak parameters using gaussian or exponential-gaussian hybrid peak fitting.}
\usage{
fit_peaks(y, pos, sd.max = 50, fit = c("egh", "gaussian"), max.iter = 1000)
}
\arguments{
  \item{y}{response (numerical vector)}
  \item{pos}{Locations of peaks in vector y.}
  \item{sd.max}{Maximum standard deviation of peak. Defaults to 50.}
  \item{fit}{Function for peak fitting. (Currently "gaussian" and 
  exponential gaussian hybrid ("egh") are supported. Defaults to egh.)}
  \item{max.iter}{Maximum number of iterations to use in nonlinear least squares peak-fitting. Defaults to 1000.}
}
\details{
  Find peaks with function \code{find_peaks} by looking for zero-crossings in
  the smoothed first derivative of a signal that exceed a given slope threshold.

  Peak parameters are calculated using \code{fit_peaks}, which fits the data to a gaussian or exponential-gaussian hybrid curve using non-linear least squares estimation as implemented in \code{\link[minpack.lm:nlsLM]{nlsLM}}. Area under the fitted curve is estimated using trapezoidal estimation.
}
\value{
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
  \item{r.squared}{r-squared value for linear fit of model to data.}

  Again, the first five elements (rt, start, end, sd and FWHM) are expressed as
  indices, so not in terms of the real retention times. The
  transformation to "real" time is done in function \code{get_peaks}.
}
\note{The \code{\link{fit_peaks}} function is adapted from Dr. Robert Morrison's \url{https://github.com/robertdouglasmorrison/DuffyTools}{DuffyTools} package as well as code published in Ron Wehrens' \url{https://github.com/rwehrens/alsace}{alsace} package.}
\author{Ethan Bass}
\references{
Lan, K. & Jorgenson, J. W. A hybrid of exponential and gaussian functions as a simple model of asymmetric chromatographic peaks. Journal of Chromatography A 915, 1-13 (2001).

Naish, P. J. & Hartwell, S. Exponentially Modified Gaussian functions - A good model for chromatographic peaks in isocratic HPLC? Chromatographia 26, 285-296 (1988).
}
\seealso{\code{\link{find_peaks}}, \code{\link{get_peaks}}}
\keyword{manip}
