\name{scan_spectrum}
\alias{scan_spectrum}
\title{Plot spectrum at given retention time}
\description{
  Plot spectrum at the specified retention time.
}
\usage{scan_spectrum(rt, chrom_list, chr = 'max', lambda = 'max',
                          plot_spectrum = T, plot_trace = T, export_spectrum = F,
                          spectrum_labels = T, verbose = F, scale_spectrum = F, ...)}
\arguments{
  \item{rt}{The retention time you wish to plot.)}
  \item{chrom_list}{A list of profile matrices, each of the same dimensions
    (timepoints times components).}
  \item{chr}{Index of chromatogram you wish to plot, or "max" if you want to automatically plot the chromatogram with the largest signal.}
  \item{lambda}{The wavelength you wish to plot the trace at if plot_chrom == T and/or the wavelength to be used for the determination of signal abundance.}
  \item{plot_spectrum}{Logical. If TRUE, plots the spectrum of the chosen peak.}
  \item{plot_trace}{Logical. If TRUE, plots the trace of the chosen peak at lambda.}
  \item{spectrum_labels}{Logical. If TRUE, plots labels on maxima in plot of spectrum}
  \item{verbose}{Logical. If TRUE, prints verbose output to console.}
  \item{scale_spectrum}{Logical. If TRUE, scales spectrum to unit height.}
  \item{export_spectrum}{Logical. If TRUE, exports spectrum to console. Defaults to FALSE.}
  \item{\dots}{Additional arguments.}
}
\author{Ethan Bass}
\seealso{\code{\link{plot_spectrum}}}
\keyword{manip}
