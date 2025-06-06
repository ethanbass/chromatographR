% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reshape_chroms.R
\name{reshape_peaktable}
\alias{reshape_peaktable}
\title{Reshape peaktable}
\usage{
reshape_peaktable(x, peaks, metadata, fixed_levels = TRUE)
}
\arguments{
\item{x}{A \code{peak_table} object.}

\item{peaks}{A character vector specifying the peaks to include. If the
character vector is named, the names of the vector elements will be used in
place of the original peak names.}

\item{metadata}{A character vector specifying the metadata fields to include.}

\item{fixed_levels}{Logical. Whether to fix factor levels of features in the
order provided. Defaults to \code{TRUE}.}
}
\value{
A data.frame containing the information for the specified 
\code{peaks} in long format.
}
\description{
Reshapes peak table from wide to long format
}
\seealso{
Other utility functions: 
\code{\link{combine_peaks}()},
\code{\link{get_lambdas}()},
\code{\link{get_times}()},
\code{\link{merge_peaks}()},
\code{\link{reshape_chroms}()}
}
\author{
Ethan Bass
}
\concept{utility functions}
