\name{components}
\alias{components}
\alias{smallComps}
\alias{removeComps}
\alias{combineComps}
\alias{suggestCompCombis}
\title{Functions to assess and refine ALS components}
\description{One of the inherent drawbacks of the MCR-ALS method is that
  in the vast majority of cases there is no one unique set of components
  describing the data, a situation known as "rotational ambiguity". This
  implies that in some cases a spectrum of a chemical compound can be
  described by a linear combination of two ALS components. This can
  sometimes be recognised by looking at elution profiles. In addition,
  in cases where the number of components is too large, some components
  may only describe noise or very small and irrelevant features. The
  functions clarified here allow one to find which components only
  correspond with minor features, to remove components, and to merge
  components. 
}
\usage{
smallComps(obj, Ithresh)
removeComps(obj, toRemove, ...)
combineComps(obj, compList, weights, ...)
suggestCompCombis(obj, indices, Ithresh = 0, corthresh = 0.9,
                  clusterHeight = 0.6)
}
\arguments{
  \item{obj}{The R object containing the als model}
  \item{Ithresh}{Intensity cutoff: all components with a maximal intensity
    (in the elution profiles) below this value will be termed "small".}
  \item{toRemove}{The indices of the components to remove from the ALS
    model. A new call to \code{doALS} will be done with the smaller set
    of components.}
  \item{\dots}{Additional arguments to \code{doALS}, e.g. \code{maxiter
      = 1} if no full set of iterations is required.}
  \item{compList}{A list indicating which components need to be
    combined. Using \code{list(c(1, c(2, 3), 4))} will lead to a
    three-component model, where components 1 and 4 are unchanged and
    components 2 and 3 are combined.}
  \item{weights}{Weights for the components to be combined. If not
    provided, equal weights will be assumed.}
  \item{indices}{A list indicating in which (groups of) samples
    correlations will be calculated. See details.}
  \item{corthresh}{Correlation threshold: components with elution
    profiles showing a higher correlation than this threshold may be
    candidates for merging.}
  \item{clusterHeight}{Similarity threshold at which to cut the
    dendrogram (see details).}
}
\details{Function \code{suggestCompCombis} checks correlations in
  elution profiles that could point to a situation where one chemical
  compound is described by two or more ALS components. For every sample
  in which this correlation is higher than the threshold, a "hit" will
  be recorded for these two components. After checking all samples and
  all combinations, the hit matrix will be used as a similarity measure
  in a hierarchical clustering. The dendrogram will be cut at a specific
  height, leading to groups of components, sometimes containing more
  than one element. In such a case, these components could be considered
  for merging.

  If injections of pure standards are present, they probably should not
  be used in isolation to check for coelution; rather, suggestions for
  combined components can be validated looking at the elution profiles
  of the standards.
}
\value{Functions \code{removeComps} and \code{combineComps} return
  \code{ALS} objects with fewer components than the original
  object. Function \code{smallComps} returns a list of two elements:
  \item{smallComps}{the indices of the small components}
  \item{maxCvalues}{the maximal values found in the concentration
    profiles across all samples for each of the components.}
}
\author{Ron Wehrens}
\examples{
data(tea)
new.lambdas <- seq(260, 500, by = 2)
tea <- lapply(tea.raw, preprocess)
tea.split <- splitTimeWindow(tea, c(12, 14), overlap = 10)

Xl <- tea.split[[3]]
Xl.opa <- opa(Xl, 10)

Xl.als <- doALS(Xl, Xl.opa)
smallC <- smallComps(Xl.als, 5)
smallC

Xl.als2 <- removeComps(Xl.als, smallC$smallC)
summary(Xl.als)
summary(Xl.als2)
## smaller models, but with a higher fit error...

## another way to decrease the number of components, this example
## not particularly deep, just to show how it can be done:
Xl.als3 <- combineComps(Xl.als, list(1, 2, 3:4, 5, c(6, 10), 6, 7:9))
summary(Xl.als3)
}
\keyword{manip}
