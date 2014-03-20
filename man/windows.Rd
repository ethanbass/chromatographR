\name{windows}
\alias{windows}
\alias{splitTimeWindow}
\alias{mergeTimeWindows}
\title{Splitting and merging of data across the time axis.}
\description{
Often MCR data sets can be analysed much more quickly and efficiently
when split into several smaller time windows. For interpretation
purposes, the results after analysis can be merged again.
}
\usage{
splitTimeWindow(datalist, splitpoints, overlap = 0)
mergeTimeWindows(obj, simSThreshold = .9, simCThreshold = .9, verbose = FALSE)
}
\arguments{
  \item{datalist}{A list of (numerical) data matrices}
  \item{splitpoints}{A numerical vector of cut points. In case the time
    axis extends beyond the range of the cut points, additional cut
    points are added at the beginning or at the end of the time axis to
    ensure that all time points are taken into account.}
  \item{overlap}{Number of points in the overlap region between two
    consecutive windows. Default: 0 (non-overlapping windows).}
  \item{obj}{Either experimental data that have been split up in
    different time windows (a list of matrices), or a list of ALS
    objects. See details section.}
  \item{simSThreshold, simCThreshold}{similarity thresholds to determine
    whether two patterns are the same (correlation). The two thresholds
    are checking the spectral and chromatographic components,
    respectively. If no overlap is present between time windows,
    \code{simCThreshold} is not used.}
  \item{verbose}{logical: print additional information?}
}
\value{
  Function \code{splitTimeWindows} splits every matrix in a list of data
  matrices into submatrices corresponding to time windows. This is
  represented as a list of lists, where each top level element is one
  time window. Such a time window can then be presented to the ALS
  algorithm.
  
  Function \code{mergeTimeWindows} can be used to merge data matrices as
  well as ALS result objects. In the first case, for each series of data
  matrices corresponding to different time windows, one big concatenated
  matrix will be returned. In the second case, exactly the same will be
  done for the residual matrices and concentration profiles in the ALS
  object. Spectral components are assumed to be different in different
  time windows, unless they have a correlation higher than
  \code{simSThreshold}, in which case they are merged. If overlapping
  time windows are used, an additional requirement is that the
  similarity between the concentration profiles in the overlap area must
  be at least \code{simCThreshold}. This similarity again is measured as
  a correlation.
}
\details{
  When splitting data files, the non-overlapping areas should be at
  least as big as the overlap areas. If not, the function stops with an
  error message. Note that the example below is only meant to show the
  use of the function: the data do not have enough time resolution to
  allow for a big overlap.}
\author{Ron Wehrens}
\examples{
## splitting and merging of data files
data(tea)
tea.split <- splitTimeWindow(tea.raw, c(12, 14))
names(tea.split)
sapply(tea.split, length)
lapply(tea.split, function(x) sapply(x, dim))
rownames(tea.split[[1]][[1]])[1:10]
rownames(tea.split[[2]][[1]])[1:10]

tea.merge <- mergeTimeWindows(tea.split)
all.equal(tea.merge, tea.raw)                    ## should be TRUE

tea.split2 <- splitTimeWindow(tea.raw, c(12, 14), overlap = 10)
lapply(tea.split2, function(x) sapply(x, dim))
tea.merge2 <- mergeTimeWindows(tea.split2)
all.equal(tea.merge2, tea.raw)                   ## should be TRUE

## merging of ALS results
data(teaMerged) 
ncomp <- ncol(teaMerged$S)
myPalette <- colorRampPalette(c("black", "red", "blue", "green"))
mycols <- myPalette(ncomp)

## show spectra - plotting only a few of them is much more clear...
plot(teaMerged, what = "spectra", col = mycols, comp.idx = c(2, 6))
legend("top", col = mycols[c(2, 6)], lty = 1, bty = "n",
       legend = paste("C", c(2, 6)))

## show concentration profiles - all six files
plot(teaMerged, what = "profiles", col = mycols)
## only the second file
plot(teaMerged, what = "profiles", mat.idx = 2, col = mycols)
legend("topleft", col = mycols, lty = 1, bty = "n",
       legend = paste("C", 1:ncol(teaMerged$S)))
## Note that components 2 and 6 are continuous across the window borders
## - these are found in all three windows
}
\keyword{manip}
