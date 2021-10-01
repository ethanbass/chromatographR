\name{plot_peaks}
\alias{plot_peaks}
\title{Function to visually assess accuracy of integration}
\description{
  Check integration by overlaying gaussian curves onto chromatogram.
}
\usage{
plot_peaks(chrome_list, peak_list, index=1, lambda=NULL,
  slope=.01, h=0.5, points=F, a=0.5)
}
\arguments{
  \item{chrom_list}{List of chromatograms (retention time x wavelength matrices)}
  \item{peak_list}{Output from the \code{findpeaks} function.}
  \item{index}{Index of chromatogram to be plotted.}
  \item{lambda}{Wavelength for plotting.}
  \item{w}{width of search window.}
  \item{slope}{Slope threshold defining peaks.}
  \item{h}{}
  \item{points}{Logical. If TRUE, plot peak maxima. Defaults to FALSE.}
  \item{a}{Select "alpha"" parameter controlling transparency of shapes.}
}
\author{Ethan Bass}
\seealso{\code{\link{get_peaks}}}
\keyword{manip}
