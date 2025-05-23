% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_peaks.R
\name{get_peaks}
\alias{get_peaks}
\title{Get peak list.}
\usage{
get_peaks(
  chrom_list,
  lambdas,
  fit = c("egh", "gaussian", "raw"),
  sd.max = 50,
  max.iter = 100,
  time.units = c("min", "s", "ms"),
  estimate_purity = FALSE,
  noise_threshold = 0.001,
  show_progress = NULL,
  cl = 2,
  collapse = FALSE,
  ...
)
}
\arguments{
\item{chrom_list}{A list of profile matrices, each of the same dimensions
(timepoints × wavelengths).}

\item{lambdas}{Character vector of wavelengths to find peaks at.}

\item{fit}{What type of fit to use. Current options are exponential-gaussian
hybrid (\code{egh}), gaussian or raw. The \code{raw} setting performs
trapezoidal integration directly on the raw data without fitting a peak shape.}

\item{sd.max}{Maximum width (standard deviation) for peaks. Defaults to 50.}

\item{max.iter}{Maximum number of iterations for non-linear least squares
in \code{\link{fit_peaks}}.}

\item{time.units}{Units of \code{sd}, \code{FWHM}, \code{area}, and \code{tau}
(if applicable). Options are minutes (\code{"min"}), seconds (\code{"s"}), or
milliseconds (\code{"ms"}).}

\item{estimate_purity}{Logical. Whether to estimate purity or not. Defaults
to \code{FALSE}. (If \code{TRUE}, this will slow down the function significantly).}

\item{noise_threshold}{Noise threshold. Argument to \code{get_purity}.}

\item{show_progress}{Logical. Whether to show progress bar. Defaults to
\code{TRUE} if \code{\link[pbapply]{pbapply}} is installed.}

\item{cl}{Argument to \code{\link[pbapply]{pblapply}} or \code{\link[parallel]{mclapply}}.
Either an integer specifying the number of clusters to use for parallel
processing or a cluster object created by \code{\link[parallel]{makeCluster}}.
Defaults to 2. On Windows integer values will be ignored.}

\item{collapse}{Logical. Whether to collapse multiple peak lists per sample
into a single list when multiple wavelengths (\code{lambdas}) are provided.}

\item{\dots}{Additional arguments to \code{\link{find_peaks}}. Arguments
provided to \code{\link{find_peaks}} can be used to fine-tune the peak-finding
algorithm. Most importantly, the \code{smooth_window} should be increased if
features are being split into multiple bins. Other arguments that can be used
here include \code{smooth_type}, \code{slope_thresh}, and \code{amp_thresh}.}
}
\value{
The result is an S3 object of class \code{peak_list}, containing a
nested list of data.frames containing information about the peaks fitted for
each chromatogram at each of wavelengths specified by the \code{lamdas}
argument. Each row in these data.frames is a peak and the columns contain
information about various peak parameters:
\itemize{
\item \code{rt}: The retention time of the peak maximum.
\item \code{start}: The retention time where the peak is estimated to begin.
\item \code{end}: The retention time where the peak is estimated to end.
\item \code{sd}: The standard deviation of the fitted peak shape.
\item \code{tau} The value of parameter \eqn{\tau}. This parameter determines
peak asymmetry for peaks fit with an exponential-gaussian hybrid function.
(This column will only appear if \code{fit = egh}.
\item \code{FWHM}: The full-width at half maximum.
\item \code{height}: The height of the peak.
\item \code{area}: The area of the peak as determined by trapezoidal approximation.
\item \code{r.squared} The coefficient of determination (\eqn{R^2}) of the fitted
model to the raw data. (\strong{Note}: this value is calculated by fitting a
linear model of the fitted peak values to the raw data. This approach is
statistically questionable, since the models are fit using non-linear least
squares. Nevertheless, it can still be useful as a rough metric for
"goodness-of-fit").
\item \code{purity} The peak purity.
}
}
\description{
Finds and fits peaks and extracts peak parameters from a list of chromatograms
at the specified wavelengths.
}
\details{
Peaks are located by finding zero-crossings in the smoothed first derivative
of the specified chromatographic traces (function \code{\link{find_peaks}}).
At the given positions, an exponential-gaussian hybrid (or regular gaussian)
function is fit to the signal using \code{\link{fit_peaks}} according to the
value of \code{fit}. Finally, the area is calculated using trapezoidal
approximation.

Additional arguments can be provided to \code{\link{find_peaks}} to fine-tune
the peak-finding algorithm. For example, the \code{smooth_window} can be
increased to prevent peaks from being split into multiple features. Overly
aggressive smoothing may cause small peaks to be overlooked.

The standard deviation (\code{sd}), full-width at half maximum (\code{FWHM}),
tau \code{tau}, and \code{area} are returned in units determined by
\code{time.units}. By default, the units are in minutes. To compare directly
with 'ChemStation' integration results, the time units should be in seconds.
}
\note{
The bones of this function are adapted from the
\href{https://github.com/rwehrens/alsace/blob/master/R/getAllPeaks.R}{getAllPeaks}
function authored by Ron Wehrens (though the underlying algorithms for peak
identification and peak-fitting are not the same).
}
\examples{
\dontshow{if (interactive()) withAutoprint(\{ # examplesIf}
data(Sa_pr)
pks <- get_peaks(Sa_pr, lambdas = c('210'), sd.max=50, fit="egh")
\dontshow{\}) # examplesIf}
}
\references{
\itemize{
\item Lan, K. & Jorgenson, J. W. 2001. A hybrid of exponential and gaussian
functions as a simple model of asymmetric chromatographic peaks. \emph{Journal of
Chromatography A} \strong{915}:1-13. \doi{10.1016/S0021-9673(01)00594-5}.
\item Naish, P. J. & Hartwell, S. 1988. Exponentially Modified Gaussian functions - A
good model for chromatographic peaks in isocratic HPLC? \emph{Chromatographia},
\strong{26}: 285-296. \doi{10.1007/BF02268168}.
\item O'Haver, Tom. Pragmatic Introduction to Signal Processing:
Applications in scientific measurement.
\url{https://terpconnect.umd.edu/~toh/spectrum/} (Accessed January, 2022).
\item Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
LC–DAD using multivariate curve resolution: the alsace package for R. \emph{
Metabolomics} \bold{11}:143-154. \doi{10.1007/s11306-014-0683-5}.
}
}
\seealso{
\code{\link{find_peaks}}, \code{\link{fit_peaks}}
}
\author{
Ethan Bass
}
