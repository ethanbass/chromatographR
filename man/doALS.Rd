\name{doALS}
\alias{doALS}
\alias{print.ALS}
\alias{plot.ALS}
\alias{summary.ALS}
\alias{getTime}
\alias{getWavelength}
\title{Wrapper function for als, plus some support functions}
\description{Wrapper function for the \code{als} function in the ALS
  package, providing a simple interface with sensible defaults for
  hyphenated data.}
\usage{
doALS(Xl, PureS, maxiter = 100)
\method{print}{ALS}(x, \dots)
\method{summary}{ALS}(object, \dots)
\method{plot}{ALS}(x, what = c("spectra", "profiles"), showWindows = TRUE,
                   mat.idx, comp.idx, xlab, ylab, main, \dots)
getTime(x)
getWavelength(x)
}
\arguments{
  \item{Xl}{a list of (numerical) data matrices. Missing values are not
    allowed.}
  \item{x, object}{an object of class ALS.}
  \item{PureS}{Initial estimates of pure spectral components.}
  \item{maxiter}{maximum number of iterations in ALS.}
  \item{what}{Show spectra or elution profiles}
  \item{showWindows}{If showing elution profiles, the window borders and
    the overlap areas between the windows can be shown (by
    default). Simply set this parameter to FALSE if this is undesired.}
  \item{mat.idx}{If showing elution profiles, one can provide the index
    of the sample(s) that should be shown. For every sample one plot
    will be made. Default is to show all.}
  \item{comp.idx}{Indices of components to be shown. Default is to show
    all components.}
  \item{xlab, ylab, main, \dots}{self-explanatory optional arguments}
}
\value{Function \code{doALS} returns an object of class "ALS", a list
  with the following fields:
  \item{CList}{a list of matrices with the elution profiles in the
    columns. Every matrix in this list corresponds with a matrix in the
    input.}
  \item{S}{a matrix with the spectral components in the columns. These
    are common for all data matrices.}
  \item{rss}{residual sum of squares.}
  \item{resid}{a list of residual matrices.}
  \item{iter}{number of iterations.}
  \item{summ.stats}{summary statistics, providing more information about
    the fit quality.}
  
  See the \code{\link{als}} function for more details; only the
  \code{summ.stats} field is not part of the original \code{als} output.}
\details{The \code{plot} method can be used to plot the spectral
  components (one plot for the model) or the elution profiles (one plot
  for each data matrix, so usually several plots). The \code{summary}
  method also returns fit statistics like LOF, R2 and RMS. Extractor
  functions \code{getTime} and \code{getWavelength} provide the vectors
  of time points and wavelengths from the \code{ALS} object.}
\author{Ron Wehrens}
\seealso{\code{\link{als}},\code{\link{showALSresult}}}
\examples{
data(tea)
new.lambdas <- seq(260, 500, by = 2)
tea <- lapply(tea.raw, preprocess, dim2 = new.lambdas)
tea.split <- splitTimeWindow(tea, c(12, 14), overlap = 10)

Xl <- tea.split[[2]]
Xl.opa <- opa(Xl, 4)

Xl.als <- doALS(Xl, Xl.opa)
Xl.als
summary(Xl.als)
plot(Xl.als, "spectra")
par(mfrow = c(1, 3))
plot(Xl.als, "profiles", ylim = c(0, 600), mat.idx = 1:3)
}
\keyword{manip}
