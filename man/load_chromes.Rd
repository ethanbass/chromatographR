\name{load_chromes}
\alias{load_chromes}
\title{Import folders of chromatograms in CSV format.}
\description{
Convenience function to import chromatograms from a list of folders. Chromatograms must be in CSV format}
\usage{
folders <- list.files(path, pattern = "EXPORT3D", full.names = T)
load_chromes(paths=folders)}
\arguments{
  \item{paths}{Paths of folders where chromatograms are stored.}
}
\value{A list of chromatograms in matrix format.}
\author{Ethan Bass}
\keyword{manip}
