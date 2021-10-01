\name{filter_peaks}
\alias{filter_peaks}
\title{Filter peak lists}
\description{
Utility function to remove peaks from a peak list, e.g. because their
intensity is too low. Currently one can filter on peak height, peak
area, and width at half maximum.}
\usage{
filterPeaks(peakList, minHeight, minArea, minWHM, maxWHM)
}
\arguments{
  \item{peakList}{A nested list of peak tables: the first level is the
    sample, and the second level is the component. Every component is
    described by a matrix where every row is one peak, and the columns
    contain information on retention time, full width at half maximum
    (FWHM), peak width, height, and area.}
  \item{minHeight}{Minimum peak height.}
  \item{minArea}{Minimum peak area.}
  \item{minWHM}{Minimal width at half maximum.}
  \item{maxWHM}{Maximum width at half maximum.}
}
\value{A peak list similar to the input peakList, but with all rows
  removed from the peak tables that are not satisfying the criteria.}
\author{Ron Wehrens}

\seealso{\code{\link{getAllPeaks}}}
\examples{
data(teaMerged)
pks <- getAllPeaks(teaMerged$CList, span = 11)
## only retain peaks with a peak height of at least 2
pks.filtered <- filterPeaks(pks, minHeight = 2)
plot(teaMerged, mat.idx = 3, what = "profiles", comp.idx = 2,
     showWindows = FALSE, col = "blue")
abline(v = pks[[3]][[2]][,"rt"], col = "gray", lty = 2)
abline(v = pks.filtered[[3]][[2]][,"rt"])
}
\keyword{manip}
