\name{Sa}
\alias{Sa}
\alias{Sa.raw}
\docType{data}
\title{HLPC-DAD data of goldenrod root extracts.}
\description{
Five HPLC-DAD data matrices of *Solidago altissima* roots extracted in 90 percent methanol. Since the raw data are given (no smoothing or baseline subtraction has been done, only subsetting of the time and wavelength axes), the object is called \code{Sa.raw}.
}
\usage{data(Sa.raw)}
\format{The UV-Vis data (\code{Sa.raw}) are given as a list of five
  matrices (time x wavelength).
}
\examples{
data(tea)
tpoints <- as.numeric(rownames(tea.raw[[1]]))
lambdas <- as.numeric(colnames(tea.raw[[1]]))
contour(tpoints, lambdas, tea.raw[[1]], col = terrain.colors(15),
        xlab = "Retention time (min.)", ylab = "Wavelength (nm)")
}
\keyword{datasets}
