% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_chroms.R
\name{plot_chroms_heatmap}
\alias{plot_chroms_heatmap}
\title{Plot chromatograms as heatmap}
\usage{
plot_chroms_heatmap(
  chrom_list,
  idx = NULL,
  lambdas,
  engine = c("base", "ggplot", "plotly"),
  show_legend = TRUE,
  xlim = NULL,
  legend_position = "topright",
  title = "",
  show_ylabs = FALSE
)
}
\arguments{
\item{chrom_list}{List of chromatograms to plot}

\item{idx}{A vector representing the names or numerical indices of the 
chromatograms to plot.}

\item{lambdas}{A character or numeric vector specifying the wavelengths to 
plot.}

\item{engine}{Plotting engine. Either \code{base} (\code{\link[graphics]{matplot}}), 
\code{\link[plotly]{plotly}}, or \link[ggplot2:ggplot2-package]{ggplot}.}

\item{show_legend}{Logical. Whether to display legend or not. Defaults to
\code{TRUE}.}

\item{xlim}{Range of x axis.}

\item{legend_position}{Position of legend.}

\item{title}{Title for plot.}

\item{show_ylabs}{Logical. Whether to show y labels. Defaults to \code{FALSE}.}
}
\value{
No return value, called for side effects.
}
\description{
Plots the specified traces from a list of chromatograms as a heatmap.
}
\details{
Adapted from \code{\link[VPdtw]{plot.VPdtw}}.
}
\section{Side effects}{

Plots the traces of the specified chromatograms \code{idx} at the specified
wavelengths \code{lambdas} as a heatmap. Plots can be produced using base 
graphics engine, \code{ggplot2}, or \code{plotly}, according to the value of 
\code{engine}.
}

\examples{
data(Sa_warp)
plot_chroms_heatmap(Sa_warp, lambdas = 210)
}
\seealso{
Other visualization functions: 
\code{\link{boxplot.peak_table}()},
\code{\link{mirror_plot}()},
\code{\link{plot.peak_list}()},
\code{\link{plot.peak_table}()},
\code{\link{plot_all_spectra}()},
\code{\link{plot_chroms}()},
\code{\link{plot_spectrum}()},
\code{\link{scan_chrom}()}
}
\author{
Ethan Bass
}
\concept{visualization functions}
