\name{fitpeaks}
\alias{findpeaks}
\alias{fitpeaks}
\title{Fit chromatographic peaks with a gaussian profile}
\description{
  Find chromatographic peaks, and fit peak parameters using a gaussian
  profile. The algorithm is extremely simple and could be replaced by a
  more sophisticated algorithm. In particular one can expect bad fits if
  peaks are overlapping significantly.}
\usage{
findpeaks(y, span = NULL)
fitpeaks(y, pos)
}
\arguments{
  \item{y}{response (numerical vector)}
  \item{span}{number of points used in the definition of what
    constitutes a "local" maximum. If not given, a default value of 20
    percent of the number of time points is used.}
  \item{pos}{locations of local maxima in vector y}
}
\details{
  Finding peaks with function \code{findpeaks} is based on the position
  of local maxima within a window of width \code{span}.

  Peak parameters are calculated using \code{fitpeaks}, assuming a
  normal distribution. Peak width is given as a standard deviation,
  calculated from the full width at half maximum (FWHM); the peak area
  is given by the ratio of the peak height and the density.
}
\value{
  Function \code{findpeaks} simply returns the locations of the local
  maxima, expressed as indices.
  
  Function \code{fitpeaks} returns a matrix, whose columns contain the
  following information:
  \item{rt}{location of the maximum of the peak (x)}
  \item{sd}{width of the peak (x)}
  \item{FWHM}{full width at half maximum (x)}
  \item{height}{height of the peak (y)}
  \item{area}{peak area}

  Again, the first three elements (rt, sd and FWHM) are expressed as
  indices, so not in terms of the real retention times. The
  transformation to "real" time is done in function \code{getAllPeaks}.
}
\note{Function \code{findpeaks} was modelled after code suggested by
  Brian Ripley on the R help list.}
\author{Ron Wehrens}
\seealso{\code{\link{getAllPeaks}}}
\examples{
data(tea)
new.lambdas <- seq(260, 500, by = 2)
tea <- lapply(tea.raw, preprocess, dim2 = new.lambdas)
tea.split <- splitTimeWindow(tea, c(12, 14), overlap = 10)

Xl <- tea.split[[2]]
Xl.opa <- opa(Xl, 4)

Xl.als <- doALS(Xl, Xl.opa)

tpoints <- getTime(Xl.als)
plot(tpoints, Xl.als$CList[[2]][,2], type = "l", col = "gray")
pk.pos <- findpeaks(Xl.als$CList[[2]][,2], span = 11)
abline(v = tpoints[pk.pos], col = 4)

pks <- fitpeaks(Xl.als$CList[[2]][,2], pk.pos)
apply(pks, 1,
      function(pkmodel) {
        lines(tpoints,
              dnorm(1:length(tpoints), pkmodel["rt"], pkmodel["sd"]) *
              pkmodel["area"],
              col = 2)
        invisible()
      })
## reasonably close fit, apart from the small peak in the middle...
}
\keyword{manip}
