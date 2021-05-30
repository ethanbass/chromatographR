\name{compare_spectra}
\alias{compare_spectra}
\title{Function to correct false zeros in peak table.}
\description{
Function that tries to correct false zeroes for a particular peak in the peaktable. In each chromatogram, the function compares all peaks within a certain radius around the focal peak on the basis of their spectral similarity to a reference spectrum.}
\usage{
compare_spectra(peak, peak_table, chrom_list, ts=new.ts, new.lambdas=new.lambdas,
                             thresh_auto=0.95, thresh_man=0.75, w=100, plot_it=F,
                             lambda='256', zeros=F, ref='max', order_by = "distance", verbose=T)}
\arguments{
  \item{peak}{Name of peak to be investigated.}
  \item{peak_table}{Peak table from \code{\link{getPeakTable}}.}
  \item{chrom_list}{Minimum spectral correlation to confirm a match.}
  \item{thresh_auto}{Defines minimum spectral similarity threshold for automatic match.}
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
