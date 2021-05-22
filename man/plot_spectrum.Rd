\name{plot_spectrum}
\alias{plot_spectrum}
\title{Plot spectrum from peak table}
\description{
  A function to plot the trace and/or spectrum for a given peak in peak table. Can be used to
  confirm the identity of a peak or check that a particular column in the peak table represents
  a single compound.
}
\usage{plot_spectrum(peak, peak_table, chrom_list)}
\arguments{
  \item{peak}{The name of the peak you wish to plot (must be in character format)}
  \item{peak_table}{The peak table (output from \code{\link{getPeakTable}} function)}
  \item{chrom_list}{A list of profile matrices, each of the same dimensions
    (timepoints times components).}
  \item{chr}{Index of chromatogram you wish to plot, or "max" if you want to automatically plot the chromatogram with the largest signal.}
  \item{lambda}{The wavelength you wish to plot the trace at if plot_chrom == T and/or the wavelength to be used for the determination of signal abundance.}
  \item{plot_spectrum}{Logical. If TRUE, plots the spectrum of the chosen peak.}
  \item{plot_trace}{Logical. If TRUE, plots the trace of the chosen peak at lambda.}
  \item{spectrum_labels}{Logical. If TRUE, plots labels on maxima in plot of spectrum}
  \item{verbose}{Logical. If TRUE, prints verbose output to console.}
  \item{scale_spectrum}{Logical. If TRUE, scales spectrum to unit height.}
}
\author{Ethan Bass}
\examples{
}
\keyword{manip}

install_github("https://github.com/ethanbass/alsace/tree/dev2")