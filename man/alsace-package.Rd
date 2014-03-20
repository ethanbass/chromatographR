\name{alsace-package}
\alias{alsace-package}
\alias{alsace}
\docType{package}
\title{Alternating Least Squares in an Analytic Chemistry Environment}
\description{Add-on to the ALS package that implements Alternating Least
  Squares (or Multivariate Curve Resolution, MCR). This implementation
  is specifically geared to data from systems like HPLC-DAD where
  measurements are always positive. In addition, it provides extra
  functionality to deal with large data sets, and additional
  postprocessing tools including visualization.  
}
\details{
\tabular{ll}{
Package: \tab alsace\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-07-08\cr
License: GPL (>= 2)
}

The main function of the package is also the one that contributes the
least new material: \code{\link{doALS}} is simply a wrapper for
\code{\link[ALS:als]{ALS}}, provided by the ALS package. More important
novel material are the visualization, preprocessing and postprocessing
functions: \code{\link{showALSresult}},
\code{\link{preprocess}}, \code{\link{windows}},
\code{\link{combineComps}},... See the
manual pages of these functions for more information.

It should be noted that the examples are only meant to illustrate the
use of the functions in the package, and do not constitute the final
analysis of the data provided. Indeed, really meaningful results can
only be obtained by either careful definition of additional constraints,
or the addition of more data.
}
\author{
Ron Wehrens

Maintainer: Ron Wehrens <ron.wehrens@fmach.it>
}
\references{
  R. Wehrens: Chemometrics with R. Springer Verlag, Heidelberg (2011)

  R. Wehrens, E. Carvalho et al.: High-throughput carotenoid profiling using
  multivariate curve resolution. Anal. Bioanal. Chem., 15:5075-5086 (2013)
}
\keyword{ package }
