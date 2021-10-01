\name{get_peaks}
\alias{getAllPeaks}
\title{Extract all peaks from the chromatographic profiles of a list of chromatograms.}
\description{
  Function to extract all peaks in a list of chromatograms. Peaks are located as local maxima within the given span (function \code{\link{find_peaks}}) and at the given positions a
  gaussian curve is fit (function \code{\link{fit_peaks}}).
}
\usage{get_peaks(chrome_list, lambdas=c('230','260'), span = NULL, eps = 1e-01)}
\arguments{
  \item{chrome_list}{A list of profile matrices, each of the same dimensions
    (timepoints times wavelengths).}
  \item{lambdas}{Character vector of wavelengths to find peaks at.}
  \item{span}{The span used for identifying local maxima in the
    individual components. If not given, the default of
    \code{\link{find_peaks}} is used.}
  \item{eps}{Minimal value for the peak width, basically used to
    eliminate peaks with zero width.}
}
\value{The result is a list, with each element corresponding to one data
  file, and containing data for the fitted peaks for each of the ALS
  components. Note that this function presents the "rt", "sd" and "FWHM"
  fields in real time units.}
\author{Ethan Bass & Ron Wehrens}
\note{Function is adapted from the \code{\link[alsace::getAllPeaks]{getAllPeaks}} function authored by Ron Wehrens.}
\keyword{manip}
