\name{plot_spectrum}
\alias{plot_spectrum}
\title{Plot spectrum from peak table}
\description{
  A function to plot the trace and/or spectrum for a given peak in peak table. Can be used to
  confirm the identity of a peak or check that a particular column in the peak table represents
  a single compound.
}
\usage{plot_spectrum(peak, peak_table, chrom_list, chr = 'max', lambda = 'max',
                          plot_spectrum = T, plot_trace = T, export_spectrum=FALSE,
                          spectrum_labels=T, verbose=T, scale_spectrum=F, ...)}
\arguments{
  \item{peak}{The name of the peak you wish to plot (must be in character format).}
  \item{peak_table}{The peak table (output from \code{\link{getPeakTable}} function).}
  \item{chrom_list}{A list of chromatograms in matrix form (timepoints x wavelengths).}
  \item{chr}{Numerical index of chromatogram you wish to plot, or "max"" if you want to automatically plot the chromatogram with the largest signal.}
  \item{lambda}{The wavelength you wish to plot the trace at if plot_chrom == T and/or the wavelength to be used for the determination of signal abundance.}
  \item{plot_spectrum}{Logical. If TRUE, plots the spectrum of the chosen peak. Defaults to TRUE.}
  \item{plot_trace}{Logical. If TRUE, plots the trace of the chosen peak at lambda. Defaults to TRUE.}
  \item{spectrum_labels}{Logical. If TRUE, plots labels on maxima in spectral plot. Defaults to TRUE.}
  \item{verbose}{Logical. If TRUE, prints verbose output to console. Defaults to TRUE.}
  \item{scale_spectrum}{Logical. If TRUE, scales spectrum to unit height. Defaults to FALSE.}
  \item{export_spectrum}{Logical. If TRUE, exports spectrum to console. Defaults to FALSE.}
  \item{...}{Additional arguments.}
}
\author{Ethan Bass}
\keyword{manip}