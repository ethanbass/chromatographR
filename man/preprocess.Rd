\name{preprocess}
\alias{preprocess}
\title{Preprocessing smooth time-wavelength data}
\description{
  Standard preprocessing of response matrices where the first axis is a
  time axis, and the second a spectral axis. An example is HPLC-DAD
  data. For smooth data, like UV-VIS data, there is the option to
  decrease the size of the matrix by interpolation. By default, the data
  are baseline-corrected in the time direction and smoothed in the
  spectral dimension.}
\usage{
preprocess(X, dim1, dim2, remove.time.baseline = TRUE,
           spec.smooth = TRUE, maxI, ...)
}
\arguments{
  \item{X}{A numerical data matrix, missing values are not allowed. If
    rownames or colnames attributes are used, they should be numerical
    and signify time points and wavelengths, respectively.}
  \item{dim1}{A new, usually shorter, set of time points (numerical). The
    range of these should not be outside the range of the original time
    points, otherwise the function stops with an error message.} 
  \item{dim2}{A new, usually shorter, set of wavelengths (numerical). The
    range of these should not be outside the range of the original
    wavelengths, otherwise the function stops with an error message.}
  \item{remove.time.baseline}{logical, indicating whether baseline
    correction should be done in the time direction. Default is TRUE.}
  \item{spec.smooth}{logical, indicating whether smoothing should be
    done in the spectral direction. Default is TRUE.}
  \item{maxI}{if given, the maximum intensity in the matrix is set to
    this value.}
  \item{\dots}{further optional arguments to the \code{baseline.corr}
    function.}
}
\value{
  The function returns the preprocessed data matrix, with rownames and
  colnames indicating the time points and wavelengths, respectively.
}
\author{Ron Wehrens}
\examples{
data(tea)
tpoints <- as.numeric(rownames(tea.raw[[1]]))
lambdas <- as.numeric(colnames(tea.raw[[1]]))

## limit retention time and wavelength ranges, and do smoothing and
## baseline correction
new.time <- seq(13, 14.1, by = .05)
new.wavelengths <- seq(400, 500, by = 2)
tea.raw1.processed <-
  preprocess(tea.raw[[1]], dim1 = new.time, dim2 = new.wavelengths)

plot(tpoints, tea.raw[[1]][,lambdas == 470],
     xlim = range(new.time), type = "l", col = "gray",
     main = "Chromatogram at 470 nm", xlab = "Time (min.)",
     ylab = "")
lines(new.time, tea.raw1.processed[,new.wavelengths == 470], col = "red")
legend("topleft", lty = 1, col = c("gray", "red"), bty = "n",
       legend = c("Original data", "Preprocessed data"))

plot(lambdas, tea.raw[[1]][tpoints == 13.7,],
     xlim = range(new.wavelengths),
     ylim = c(0, max(tea.raw[[1]][tpoints == 13.7,])),
     type = "l", col = "gray",
     main = "Spectrum at 13.7 min.", xlab = expression(lambda),
     ylab = "")
lines(new.wavelengths, tea.raw1.processed[new.time == 13.7,], col = "red")
legend("topleft", lty = 1, col = c("gray", "red"), bty = "n",
       legend = c("Original data", "Preprocessed data"))
}
\keyword{manip}
