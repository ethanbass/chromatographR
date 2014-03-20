\name{tea}
\alias{tea}
\alias{tea.raw}
\docType{data}
\title{HLPC-DAD data for grape extracts conserved with TEA}
\description{
Five (very much compressed) HPLC-DAD data matrices of grape extracts
after several storage times. All extracts come from the same pooled
sample. Since the raw data are given (no smoothing or
baseline subtraction has been done, only subsetting of the time and
wavelength axes), the object is called \code{tea.raw}.

}
\usage{data(tea)}
\format{The UV-Vis data (\code{tea.raw}) are given as a list of five
  matrices, each of dimension 97 times 209 (time x wavelength). The
  names of the list indicate the day of measurement - day 0 is
  represented by two measurements.
}
\source{
 Provided by Elisabete Carvalho.
}
\references{
  This is part of the data that have been used in:
  R. Wehrens, E. Carvalho, D. Masuero, A. de Juan and S. Martens:
  High-throughput carotenoid profiling using multivariate curve
  resolution. Anal. Bioanal. Chem. 15:5057-5086 (2013)
}
\examples{
data(tea)
tpoints <- as.numeric(rownames(tea.raw[[1]]))
lambdas <- as.numeric(colnames(tea.raw[[1]]))
contour(tpoints, lambdas, tea.raw[[1]], col = terrain.colors(15),
        xlab = "Retention time (min.)", ylab = "Wavelength (nm)")
}
\keyword{datasets}
