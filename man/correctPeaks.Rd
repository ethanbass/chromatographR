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
\examples{
data(teaMerged)
pks <- getAllPeaks(teaMerged$CList, span = 11)
warping.models <- correctRT(teaMerged$CList, reference = 2,
                            what = "models")
pks.corrected <- correctPeaks(pks, warping.models)

## original profiles and peaks, in black and gray
plot(teaMerged, mat.idx = 3, what = "profiles", comp.idx = 2,
     showWindows = FALSE, col = "gray")
abline(v = pks[[3]][[2]][,"rt"])
## shifted profiles and peaks, in red and pink
CList.corrected <- correctRT(teaMerged$CList, reference = 2)
lines(as.numeric(rownames(CList.corrected[[3]])),
      CList.corrected[[3]][,2], col = "pink")
abline(v = pks.corrected[[3]][[2]][,"rt.cor"], col = "red")
## note that the rightmost peak in the uncorrected data is no longer
## within the range of the data
}
\keyword{manip}
