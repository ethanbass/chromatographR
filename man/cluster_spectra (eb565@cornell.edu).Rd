\name{cluster_spectra}
\alias{cluster_spectra}
\title{Function to cluster peaks by spectral similarity.}
\description{
Function that clusters peaks by spectral similarity.}
\usage{
cluster_spectra <- function(pkTab, chrom_list, 
                            deepSplit = FALSE, mn=5, mx=100,
                            alpha=0.95, nboot=1000, plot_dend=T, plot_spectra=T,
                            verbose=T, parallel=T)
                             }
\arguments{
  \item{pkTab}{Peak table from \code{\link{getPeakTable}}.}
  \item{chrom_list}{A list of chromatograms in matrix form (timepoints x wavelengths).}
  \item{mn}{Minimum threshold for the number of peaks a cluster may have.}
  \item{thresh_man}{Defines minimum spectral similarity threshold to trigger manual matching process.}
  \item{w}{Defines radius around focal peak to search for matches.}
  \item{plot_it}{Logical. If TRUE, plots .}
  \item{lambda}{Wavelength for plot.}
  \item{zeros}{Logical. If TRUE, only check zero rows in peak table.}
  \item{ref}{Defines chromatogram for reference spectrum. If "max" choose chromatogram with
  maximum absorbance.}
  \item{order_by}{If "distance", will order matches by distance from theoretical retention time.}
  \item{verbose}{Logical. If TRUE, prints verbose output to console.}
}
\value{A peak table similar to the input peak table, but with changes to entries selected column
as defined by the algorithm.}
\author{Ethan Bass}

\seealso{\code{\link{getAllPeaks}}}
\keyword{manip}
