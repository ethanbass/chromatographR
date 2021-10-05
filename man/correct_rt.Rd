\name{correct_rt}
\alias{correct_rt}
\alias{correctRT}
\alias{correctRT5}
\title{Retention time correction}
\description{
  Correction of retention time differences using parametric time warping.
}
\usage{correct_rt(chrom_list, models=NULL, lambdas, reference='best', 
        what = c("models", "corrected.values"), init.coef = c(0, 1, 0), n.traces=NULL,
        n.zeros=0, scale=T, trwdth=200, ...)}
\arguments{
  \item{chrom_list}{List of matrices containing concentration profiles.}
  \item{models}{List of models to warp by.}
  \item{lambdas}{Select wavelengths to use by name.}
  \item{reference}{Index of the sample that is to be considered the
    reference sample.}
  \item{what}{What to return: either the time-corrected profiles (useful
    for visual inspection) or the warping models (for further
    programmatic use).}
  \item{init.coef}{Starting values for the optimisation.}
  \item{n.traces}{Number of traces to use.}
  \item{n.zeros}{Number of zeros to add.}
  \item{scale}{Logical. If true, scale chromatograms before warping.}
  \item{trwdth}{width of the triangle in the WCC criterion.}
  \item{\dots}{Optional arguments for the \code{ptw} function. The only
    argument that cannot be changed is \code{warp.type}: this is always
    equal to \code{"global"}.}
}
\value{A list of `ptw` objects or a list of warped absorbance profiles.}
\author{Ron Wehrens, Ethan Bass}
\seealso{\code{\link{ptw}}, \code{\link{correctPeaks}}}
\examples{
data(teaMerged)
CList.corrected <- correctRT(teaMerged$CList, reference = 2)

original.profiles <- sapply(teaMerged$CList, identity, simplify = "array")
corrected.profiles <- sapply(CList.corrected, identity, simplify = "array")

def.par <- par(no.readonly = TRUE)
par(mfrow = c(2,4))
for (i in 1:4)
    matplot(dimnames(original.profiles)[[1]],
            original.profiles[,i,], type = "l", lty = 1,
            xlab = "Time (min.)", ylab = "Response",
            main = paste("Component", i))
for (i in 1:4)
    matplot(dimnames(original.profiles)[[1]],
            corrected.profiles[,i,], type = "l", lty = 1,
            xlab = "Time (min.)", ylab = "Response",
            main = paste("Component", i, "- warped"))
par(def.par) ## reset defaults
}
\references{
Eilers, P.H.C. 2004. \href{https://pubs.acs.org/doi/10.1021/ac034800e}{Parametric Time Warping.} \emph{Anal. Chem.} \bold{76}:404-411.

Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. \href{https://academic.oup.com/bioinformatics/article/31/18/3063/240579}{Fast parametric time warping of peak lists.} \emph{Bioinformatics} \bold{31}:3063-3065.
}
\keyword{manip}
