\name{getPeakTable}
\alias{getPeakTable}
\title{Convert MCR results into an ordered peak table}
\description{
Function returns a matrix of intensities, where rows correspond to
(aligned) features and columns to objects (samples, injections,
...). The function performs a complete linkage clustering of retention
times across all samples, and cuts at a height given by the user (which
can be interpreted as the maximal inter-cluster retention time
difference). If two peaks from the same sample are assigned to the same
cluster, and error message is given.
}
\usage{
getPeakTable(peakList, response = c("area", "height"),
             use.cor = TRUE, maxdiff = 0.2, plotIt = FALSE, ask = plotIt)
}
\arguments{
  \item{peakList}{A nested list of peak tables: the first level is the
    sample, and the second level is the component. Every component is
    described by a matrix where every row is one peak, and the columns
    contain information on retention time, full width at half maximum
    (FWHM), peak width, height, and area.}
  \item{response}{An indicator whether peak area or peak height is to be
    used as intensity measure. Default is peak area.}
  \item{use.cor}{Logical, indicating whether to use corrected retention
    times (by default) or raw retention times (not advised!).}
  \item{maxdiff}{Height at which the complete linkage dendrogram will be
    cut. Can be interpreted as the maximal inter-cluster retention time
    difference.}
  \item{plotIt}{Logical. If TRUE, for every component a stripplot will
    be shown indicating the clustering.}
  \item{ask}{Logical. Ask before showing new plot?}
}
\value{The function returns a data frame where the first couple of
  columns contain meta-information on the features (component, peak,
  retention time) and the other columns contain the intensities of the
  features in the individual injections.}
\details{If one sees warnings about peaks from the same sample sharing a
  cluster label, one option is to reduce the \code{maxdiff} variable -
  this, however, will increase the number of clusters. Another option is
  to filter the peaks on intensity: perhaps one of the two peaks in the
  cluster is only a very small feature.}
\author{Ron Wehrens}
\examples{
data(teaMerged)
pks <- getAllPeaks(teaMerged$CList, span = 11)
warping.models <- correctRT(teaMerged$CList, reference = 2,
                            what = "models")
pks.corrected <- correctPeaks(pks, warping.models)
pkTab <- getPeakTable(pks.corrected, response = "area")
}
\keyword{manip}
