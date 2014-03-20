\name{teaMerged}
\alias{teaMerged}
\docType{data}
\title{Results of an ALS analysis on individual windows}
\description{
  Object of class \code{ALS}: the result of the analysis of the
  \code{tea} data, using three time windows with an overlap parameter of
  10. The three ALS models have been merged into one ALS object, which
  can be inspected and used for further analysis.
}
\usage{data(teaMerged)}
\examples{
## generation of the data
data(tea)
new.lambdas <- seq(260, 500, by = 2)
tea <- lapply(tea.raw,
              preprocess,
              dim2 = new.lambdas)
tea.split <- splitTimeWindow(tea, c(12, 14), overlap = 10)
tea.alslist <- lapply(tea.split,
                      function(Xl) {
                        Xl.opa <- opa(Xl, 4)
                        doALS(Xl, Xl.opa)
                      })
teaMerged <- mergeTimeWindows(tea.alslist)
## This is the object saved in teaMerged.RData

ncomp <- ncol(teaMerged$S)
myPalette <- colorRampPalette(c("black", "red", "blue", "green"))
mycols <- myPalette(ncomp)

plot(teaMerged, what = "spectra", col = mycols)
legend("top", col = mycols, lty = 1, bty = "n", ncol = 2,
       legend = paste("C", 1:ncol(teaMerged$S)))
}
\keyword{datasets}
