\name{cluster_spectra}
\alias{cluster_spectra}
\title{Function to cluster peaks by spectral similarity.}
\description{
Function that clusters peaks by spectral similarity.}
\usage{
cluster_spectra <- function(pkTab, chrom_list, 
                            deepSplit = FALSE, peak_no=c(5,100),
                            alpha=0.95, nboot=1000, plot_dend=T, plot_spectra=T,
                            verbose=T, save=T, parallel=T)
                             }
\arguments{
  \item{pkTab}{Peak table from \code{\link{getPeakTable}}.}
  \item{chrom_list}{A list of chromatograms in matrix form (timepoints x wavelengths).}
  \item{peak_no}{Minimum and maximum thresholds for the number of peaks a cluster may have.}
  \item{alpha}{Confidence threshold for inclusion of cluster.}
  \item{nboot}{Number of bootstrap replicates for \code{\link{pvclust::pvclust}}.}
  \item{parallel}{Logical. If TRUE, use parallel processing for \code{\link{pvclust::pvclust}}.}
  \item{plot_dend}{Logical. If TRUE, plots dendrogram with bootstrap values.}
  \item{plot_spectra}{Logical. If TRUE, plots overlapping spectra for each cluster.}
  \item{verbose}{Logical. If TRUE, prints progress report to console.}
  \item{save}{Logical. If TRUE, saves pvclust object to current directory.}
}
\value{S4 "cluster" object.}
\author{Ethan Bass}

\seealso{\code{\link{pvclust::pvclust}}}
