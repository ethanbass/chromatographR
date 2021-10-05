\name{correctPeaks}
\alias{correctPeaks}
\title{Correct peak positions according to a ptw warping model}
\description{Once an appropriate warping model has been established,
  corrected retention times can be predicted for each peak. These are
  stored in a separate column in the list of peak tables.
}
\usage{
correctPeaks(peakList, modList)
}
\arguments{
  \item{peakList}{A nested list of peak tables: the first level is the
    sample, and the second level is the component. Every component is
    described by a matrix where every row is one peak, and the columns
    contain information on retention time, full width at half maximum
    (FWHM), peak width, height, and area.}
  \item{modList}{A list of ptw models.}
}
\value{The input list of peak tables is returned with extra columns
  containing the corrected retention time.}
\author{Ron Wehrens}
\seealso{\code{\link{correctRT}}}
\keyword{manip}
