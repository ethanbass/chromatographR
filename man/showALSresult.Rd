\name{showALSresult}
\alias{showALSresult}
\title{Plot ALS results in a more elaborate way}
\description{
  Simultaneous visualization of pure components (spectra and time
  profiles) and either raw data, fitted data or residuals.
}
\usage{
showALSresult(xals, xlst,
              tp = getTime(xals), wl = getWavelength(xals),
              mat.idx = 1:length(xlst),
              img.col = terrain.colors(10), zlim, xlab, ylab,
              compound.col = 1:ncol(xals$S), logsc = TRUE,
              plotPureC = c("both", "spec", "conc", "none"),
              titles, annotation = show.img,
              PureChght = 0.33, PureCwdth = min(nplot, 5)/5 - 0.1,
              show.img = TRUE)
}
\arguments{
  \item{xals}{The fitted ALS object.}
  \item{xlst}{The data: the list of matrices on which the ALS object is
    based.}
  \item{tp}{Optional vector of time points. If missing, will be
    determined using function \code{getTime}.}
  \item{wl}{Optional vector of wavelengths. If missing, will be
    determined using function \code{getWavelength}.}
  \item{mat.idx}{Indices of the samples to be visualized.}
  \item{img.col}{Color vector for \code{image}.}
  \item{zlim}{Range of the image colors.}
  \item{xlab, ylab}{Axis annotation strings.}
  \item{compound.col}{Colors to be used for components in the pure
    specta/profile plots.}
  \item{logsc}{Logical, indicating whether the images should be
    logscaled vefore visualization. Default: TRUE.}
  \item{plotPureC}{Determines which part (if any) of the pure components
    is shown.}
  \item{titles}{Titles for the plots for the individual samples. If not
    given, the names of the \code{xlst} elements will be used.}
  \item{annotation}{if a text string, this will be shown in the top
    right corner panel. If anything else but FALSE, a color bar will be
    drawn. The default is to show the color bar for images and not to
    show it for contour plots.}
  \item{PureChght}{Height, relative to the height of the data panels, of
  the top row of pure concentration profiles}
  \item{PureCwdth}{Width, relative to the width of the data panels, of
  the right column of pure spectra}
  \item{show.img}{logical, indicating whether \code{image} is used (the
    default) or \code{contour}}
}
\author{Ron Wehrens}
\seealso{\code{\link{plot.ALS}}}
\examples{
data(teaMerged)

ncomp <- ncol(teaMerged$S)
myPalette <- colorRampPalette(c("black", "red", "blue", "green"))
mycols <- myPalette(ncomp)

maxResid <- max(abs(range(teaMerged$resid)))
showALSresult(teaMerged, teaMerged$resid,
              compound.col = mycols, logsc = FALSE, img.col = cm.colors(9),
              mat.idx = 2:4, zlim =c(-maxResid, maxResid))
}
\keyword{manip}
