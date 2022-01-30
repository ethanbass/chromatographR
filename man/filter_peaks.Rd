\name{filter_peaks}
\alias{filter_peaks}
\title{Filter peak lists}
\description{
Utility function to remove peaks from a peak list, e.g. because their
intensity is too low. Currently one can filter on peak height, peak
area, and width at half maximum.}
\usage{
filter_peaks(peak_list, min_height, min_area, min_sd, max_sd)
}
\arguments{
  \item{peak_list}{A nested list of peak tables: the first level is the
    sample, and the second level is the component. Every component is
    described by a matrix where every row is one peak, and the columns
    contain information on retention time, full width at half maximum
    (FWHM), peak width, height, and area.}
  \item{min_height}{Minimum peak height.}
  \item{min_area}{Minimum peak area.}
  \item{min_sd}{Minimal standard deviation.}
  \item{max_sd}{Maximum standard deviation.}
}
\value{A peak list similar to the input peakList, but with all rows
  removed from the peak tables that are not satisfying the criteria.}
\author{Ron Wehrens, EthanBass}

\seealso{\code{\link{get_peaks}}}
\keyword{manip}
