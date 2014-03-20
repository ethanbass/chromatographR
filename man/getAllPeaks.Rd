\name{getAllPeaks}
\alias{getAllPeaks}
\title{Extract all peaks from the chromatographic profiles of an ALS object}
\description{
  Extractor function to find all peaks in the chromatographic profiles
  of an ALS object. Peaks are located as local maxima within the given
  span (function \code{\link{findpeaks}}) and at the given positions a
  gaussian curve is fit (function \code{\link{fitpeaks}}).
}
\usage{getAllPeaks(CList, span = NULL)}
\arguments{
  \item{CList}{A list of profile matrices, each of the same dimensions
    (timepoints times components).}
  \item{span}{The span used for identifying local maxima in the
    individual components. If not given, the default of
    \code{\link{findpeaks}} is used.}
}
\value{The result is a list, with each element corresponding to one data
  file, and containing data for the fitted peaks for each of the ALS
  components. Note that this function presents the "rt", "sd" and "FWHM"
  fields in real time units.}
\author{Ron Wehrens}
\examples{
data(teaMerged)
pks <- getAllPeaks(teaMerged$CList, span = 11)
## show component 2 from the second file
par(mfrow = c(2,1))
plot(teaMerged, what = "profiles", showWindows = FALSE,
     mat.idx = 2, comp.idx = 2)
## and show where the peaks are picked
abline(v = pks[[2]][[2]][,"rt"], col = "gray")

## same for component 6
plot(teaMerged, what = "profiles", showWindows = FALSE,
     mat.idx = 2, comp.idx = 6, col = "red")
abline(v = pks[[2]][[6]][,"rt"], col = "pink")
}
\keyword{manip}
