\name{cluster_spectra}
\alias{cluster_spectra}
\title{Function to cluster peaks by spectral similarity.}
\description{
Function that clusters peaks by spectral similarity. A representative spectrum is selected for each peak in the provided peak table and used to construct a distance matrix based on spectral similarity (pearson correlation) between peaks. Hierarchical clustering with bootstrap resampling is performed on the resulting correlation matrix, as implemented in the \code{\link[pvclust:pvclust]{pvclust}} package. Bootstrap values can be used to select clusters that exceed a certain confidence threshold as defined by alpha.
}
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
  \item{nboot}{Number of bootstrap replicates for \code{\link[pvclust::pvclust]{pvclust}}.}
  \item{parallel}{Logical. If TRUE, use parallel processing for \code{\link[pvclust::pvclust]{pvclust}}.}
  \item{plot_dend}{Logical. If TRUE, plots dendrogram with bootstrap values.}
  \item{plot_spectra}{Logical. If TRUE, plots overlapping spectra for each cluster.}
  \item{verbose}{Logical. If TRUE, prints progress report to console.}
  \item{save}{Logical. If TRUE, saves pvclust object to current directory.}
}
\value{Returns S4 "cluster" object with the following components:
\item{peaks}{a character vector containing the names of all peaks contained in the given cluster.}
  \item{pval}{a numeric vector of length 1 containing the boostrap p-value (au) for the given cluster.}
}
\author{Ethan Bass}

\references{
  R. Suzuki, H. Shimodaira: Pvclust: an R package for assessing the uncertainty in hierarchical clustering. Bioinformatics, 22-12:1540-1542 (2006). 
}