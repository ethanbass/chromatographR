\name{opa}
\alias{opa}
\title{
Finding the most dissimilar variables in a data matrix: the Orthogonal
Projection Approach}
\description{
  This function finds the set of most dissimilar rows in a data
  matrix. If no initial selection is presented, the first object is
  selected by comparison with the vector of column means. As a distance
  function the determinant of the crossproduct matrix is used.
}
\usage{
opa(x, ncomp, initXref = NULL)
}
\arguments{
  \item{x}{Data matrix (numerical). May not contain missing values.}
  \item{ncomp}{Number of rows to be selected.}
  \item{initXref}{Optional matrix to be expanded - should be a subset of
    the rows to select.}
}
\value{
  The function returns a submatrix of X, where the columns contain the
  (unit-length scaled) spectra from the input data that are most dissimilar.
}
\references{
  F. Questa Sanchez et al.: Algorithm for the assessment of peak purity
  in liquid chromatography with photodiode-array detection. Analytica
  Chimica Acta 285:181-192 (1994)

  R. Wehrens: Chemometrics with R. Springer Verlag, Heidelberg (2011)
}
\author{Ron Wehrens}
\examples{
data(tea)

tea <- lapply(tea.raw, preprocess, maxI = 100)

ncomp <- 7
spectra <- opa(tea, ncomp)

myPalette <- colorRampPalette(c("black", "red", "blue", "green"))
mycols <- myPalette(ncomp)
matplot(as.numeric(rownames(spectra)), spectra, type = "l", lty = 1,
        xlab = expression(lambda), ylab = "", col = mycols)
legend("topright", legend = paste("Comp.", 1:ncomp), col = mycols,
       lty = 1, ncol = 2, bty = "n")
}
\keyword{manip}
