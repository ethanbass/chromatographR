\name{plot_peaks}
\alias{plot_peaks}
\title{Function to visually assess accuracy of integration}
\description{
  Check integration by overlaying gaussian curves onto chromatogram.
}
\usage{
plot_peaks(chrome_list, peak_list, index=1, lambda='230', w=100, slope=.01, points=F, a=0.5, time=c('rt','raw'),fit=c("gaussian","egh"))
}
\arguments{
  \item{chrome_list}{List of chromatograms (retention time x wavelength matrices)}
  \item{peak_list}{Output from the \code{findpeaks} function.}
  \item{index}{Index of chromatogram to be plotted.}
  \item{lambda}{Wavelength for plotting.}
  \item{w}{Half-width of search window. Should match value used in \code{fitpeaks}.}
  \item{slope}{Slope threshold defining peaks.}
  \item{points}{Logical. If TRUE, plot peak maxima. Defaults to FALSE.}
  \item{a}{Select "alpha"" parameter controlling transparency of the gaussian shapes.}
  \item{time}{Choose "raw" indices or "rt" for conversion to "real" units (e.g. minutes).}
  \item{fit}{Currently only gaussian fit is supported.}
}
\author{Ethan Bass}
\seealso{\code{\link{getAllPeaks}}}
\keyword{manip}
