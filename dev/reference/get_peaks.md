# Get peak list.

Finds and fits peaks and extracts peak parameters from a list of
chromatograms at the specified wavelengths.

## Usage

``` r
get_peaks(
  chrom_list,
  lambdas,
  fit = c("egh", "gaussian", "raw"),
  sd_max = 50,
  max_iter = 100,
  time_unit = c("min", "s", "ms"),
  estimate_purity = FALSE,
  noise_threshold = 0.001,
  show_progress = NULL,
  cl = 2,
  collapse = FALSE,
  time.units = NULL,
  sd.max = NULL,
  max.iter = NULL,
  ...
)
```

## Arguments

- chrom_list:

  A list of profile matrices, each of the same dimensions (timepoints ×
  wavelengths).

- lambdas:

  A character or numeric vector specifying the wavelengths to find peaks
  at. For one-dimensional chromatograms, this argument can be ignored.

- fit:

  What type of fit to use. Current options are exponential-gaussian
  hybrid (`egh`), gaussian or raw. The `raw` setting performs
  trapezoidal integration directly on the raw data without fitting a
  peak shape.

- sd_max:

  Maximum width (standard deviation) for peaks. Defaults to 50.

- max_iter:

  Maximum number of iterations for non-linear least squares in
  [`fit_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/fit_peaks.md).

- time_unit:

  Units of `sd`, `FWHM`, `area`, and `tau` (if applicable). Options are
  minutes (`"min"`), seconds (`"s"`), or milliseconds (`"ms"`).

- estimate_purity:

  Logical. Whether to estimate purity or not. Defaults to `FALSE`. (If
  `TRUE`, this will slow down the function significantly).

- noise_threshold:

  Noise threshold. Argument to `get_purity`.

- show_progress:

  Logical. Whether to show progress bar. Defaults to `TRUE` if
  [`pbapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  is installed.

- cl:

  Argument to
  [`pblapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  or [`mclapply`](https://rdrr.io/r/parallel/mclapply.html). Either an
  integer specifying the number of clusters to use for parallel
  processing or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to 2. On Windows integer values will be ignored.

- collapse:

  Logical. Whether to collapse multiple peak lists per sample into a
  single list when multiple wavelengths (`lambdas`) are provided.

- time.units:

  The `time.units` argument is deprecated. Please use `time_unit`
  instead.

- sd.max:

  The `sd.max` argument is deprecated. Please use `sd_max` instead.

- max.iter:

  The `max.iter` argument is deprecated. Please use `max_iter` instead.

- ...:

  Additional arguments to
  [`find_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/find_peaks.md).
  Arguments provided to
  [`find_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/find_peaks.md)
  can be used to fine-tune the peak-finding algorithm. Most importantly,
  the `smooth_window` should be increased if features are being split
  into multiple bins. Other arguments that can be used here include
  `smooth_type`, `slope_thresh`, and `amp_thresh`.

## Value

The result is an S3 object of class `peak_list`, containing a nested
list of data.frames containing information about the peaks fitted for
each chromatogram at each of wavelengths specified by the `lamdas`
argument. Each row in these data.frames is a peak and the columns
contain information about various peak parameters:

- `rt`: The retention time of the peak maximum.

- `start`: The retention time where the peak is estimated to begin.

- `end`: The retention time where the peak is estimated to end.

- `sd`: The standard deviation of the fitted peak shape.

- `tau` The value of parameter \\\tau\\. This parameter determines peak
  asymmetry for peaks fit with an exponential-gaussian hybrid function.
  (This column will only appear if `fit = egh`.

- `FWHM`: The full-width at half maximum.

- `height`: The height of the peak.

- `area`: The area of the peak as determined by trapezoidal
  approximation.

- `r.squared` The coefficient of determination (\\R^2\\) of the fitted
  model to the raw data. (**Note**: this value is calculated by fitting
  a linear model of the fitted peak values to the raw data. This
  approach is statistically questionable, since the models are fit using
  non-linear least squares. Nevertheless, it can still be useful as a
  rough metric for "goodness-of-fit").

- `purity` The peak purity as estimated by
  [`get_purity`](https://ethanbass.github.io/chromatographR/dev/reference/get_purity.md).

## Details

Peaks are located by finding zero-crossings in the smoothed first
derivative of the specified chromatographic traces (function
[`find_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/find_peaks.md)).
At the given positions, an exponential-gaussian hybrid (or regular
gaussian) function is fit to the signal using
[`fit_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/fit_peaks.md)
according to the value of `fit`. Finally, the area is calculated using
trapezoidal approximation.

Additional arguments can be provided to
[`find_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/find_peaks.md)
to fine-tune the peak-finding algorithm. For example, the
`smooth_window` can be increased to prevent peaks from being split into
multiple features. Overly aggressive smoothing may cause small peaks to
be overlooked.

The standard deviation (`sd`), full-width at half maximum (`FWHM`), tau
`tau`, and `area` are returned in units determined by `time_unit`. By
default, the units are in minutes. To compare directly with
'ChemStation' integration results, the time units should be changed to
seconds.

## Note

The bones of this function are adapted from the
[getAllPeaks](https://github.com/rwehrens/alsace/blob/master/R/getAllPeaks.R)
function authored by Ron Wehrens (though the underlying algorithms for
peak identification and peak-fitting are not the same).

## References

- Lan, K. & Jorgenson, J. W. 2001. A hybrid of exponential and gaussian
  functions as a simple model of asymmetric chromatographic peaks.
  *Journal of Chromatography A* **915**:1-13.
  [doi:10.1016/S0021-9673(01)00594-5](https://doi.org/10.1016/S0021-9673%2801%2900594-5)
  .

- Naish, P. J. & Hartwell, S. 1988. Exponentially Modified Gaussian
  functions - A good model for chromatographic peaks in isocratic HPLC?
  *Chromatographia*, **26**: 285-296.
  [doi:10.1007/BF02268168](https://doi.org/10.1007/BF02268168) .

- O'Haver, Tom. Pragmatic Introduction to Signal Processing:
  Applications in scientific measurement.
  <https://terpconnect.umd.edu/~toh/spectrum/> (Accessed January, 2022).

- Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
  LC–DAD using multivariate curve resolution: the alsace package for R.
  *Metabolomics* **11**:143-154.
  [doi:10.1007/s11306-014-0683-5](https://doi.org/10.1007/s11306-014-0683-5)
  .

## See also

[`find_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/find_peaks.md),
[`fit_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/fit_peaks.md)

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa_pr)
pks <- get_peaks(Sa_pr, lambdas = c('210'), sd_max = 50, fit = "egh")
}
```
