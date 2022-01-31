\name{plot_peaks}
\alias{plot_peaks}
\title{Function to visually assess accuracy of integration}
\description{
  Check integration by overlaying gaussian curves onto chromatogram.
}
\usage{
plot_peaks(chrom_list, peak_list, index=1, lambda=NULL, points=F,ticks=F, a=0.5)
}
\arguments{
  \item{chrom_list}{List of chromatograms (retention time x wavelength matrices)}
  \item{peak_list}{Output from the \code{findpeaks} function.}
  \item{index}{Index of chromatogram to be plotted.}
  \item{lambda}{Wavelength for plotting.}
  \item{points}{Logical. If TRUE, plot peak maxima. Defaults to FALSE.}
  \item{ticks}{Logical. If TRUE, mark beginning and end of each peak. Defaults to FALSE.}
  \item{a}{Select "alpha"" parameter controlling transparency of shapes.}
  \item{cex.points}{Size of points. Defaults to 0.5}
  \item{\dots}{Additional arguments to plot function.}
}
\author{Ethan Bass}
\seealso{\code{\link{get_peaks}}}
\keyword{manip}
