\name{getAllPeaks}
\alias{getAllPeaks}
\title{Extract all peaks from the chromatographic profiles of an ALS object}
\description{
  Function to extract all peaks in a list of chromatograms. Peaks are located as local maxima within the given span (function \code{\link{findpeaks}}) and at the given positions a
  gaussian curve is fit (function \code{\link{fitpeaks}}).
}
\usage{getAllPeaks(chrome_list, lambdas=c('230','260'), span = NULL, eps = 1e-01)}
\arguments{
  \item{chrome_list}{A list of profile matrices, each of the same dimensions
    (timepoints times wavelengths).}
  \item{lambdas}{Character vector of wavelengths to find peaks at.}
  \item{span}{The span used for identifying local maxima in the
    individual components. If not given, the default of
    \code{\link{findpeaks}} is used.}
  \item{eps}{Minimal value for the peak width, basically used to
    eliminate peaks with zero width.}
}
\value{The result is a list, with each element corresponding to one data
  file, and containing data for the fitted peaks for each of the ALS
  components. Note that this function presents the "rt", "sd" and "FWHM"
  fields in real time units.}
\author{Ethan Bass}
\note{Function is adapted from the \code{\link[alsace::GetAllPeaks]{GetAllPeaks}} function authored by Ron Wehrens.}
\keyword{manip}
