# Get peak list.

Detects chromatographic peaks, fits peak models, and extracts peak
parameters at one or more wavelengths.

## Usage

``` r
get_peaks(
  chrom_list,
  lambdas,
  fit = c("egh", "bemg", "gaussian", "raw"),
  sd_max = 50,
  max_iter = 100,
  time_unit = c("min", "s", "ms"),
  baseline = c("none", "flat", "sloped"),
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

  A list of chromatograms in matrix format (timepoints x wavelengths).

- lambdas:

  A character or numeric vector specifying the wavelengths to find peaks
  at. For one-dimensional chromatograms, this argument can be ignored.

- fit:

  What type of fit to use. Current options are bidirectional
  exponentially modified gaussian (`"bemg"`), exponential-gaussian
  hybrid (`"egh"`), `"gaussian"`, or `"raw"`. The `raw` setting performs
  trapezoidal integration directly on the raw data without fitting a
  model to the data. Defaults to `egh`.

- sd_max:

  Maximum width (standard deviation) for peaks. Defaults to `50`.

- max_iter:

  Maximum number of iterations for non-linear least squares in
  [`fit_peaks`](https://ethanbass.github.io/chromatographR/reference/fit_peaks.md).

- time_unit:

  Specifies the units for `sd`, `FWHM`, `area`, and `tau` (if
  applicable). Options are minutes (`"min"`), seconds (`"s"`), or
  milliseconds (`"ms"`). Defaults to `"min"`.

- baseline:

  Whether to fit a baseline offset beneath each peak model. Options are
  `"none"` (no baseline, default), `"flat"` (constant offset per peak),
  or `"sloped"` (linearly varying offset per peak, useful for gradient
  elution). Not available when `fit = "raw"`.

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

  Either an integer specifying the number of cores to use for parallel
  processing or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to `2`. On Windows integer values will be ignored.

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
  [`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md).
  Arguments provided to `find_peaks` can be used to fine-tune the
  peak-finding algorithm. Most importantly, the `smooth_window` should
  be increased if features are being split into multiple bins. Other
  arguments that can be used here include `smooth_type`, `slope_thresh`,
  and `amp_thresh`.

## Value

The result is an S3 object of class `peak_list`, containing a nested
list of data.frames containing information about the peaks fitted for
each chromatogram at each of wavelengths specified by the `lambdas`
argument. Each row in these data.frames is a peak and the columns
contain information about various peak parameters:

- `rt`: The retention time of the peak maximum.

- `start`: The retention time where the peak is estimated to begin.

- `end`: The retention time where the peak is estimated to end.

- `sd`: The standard deviation (\\\sigma\\) of the fitted peak shape.

- `tau`: Empirical skewness parameter (in units of time) determining
  peak asymmetry when `fit == "egh"`.

- `tau_right`: Exponential rate constant controlling right-sided tailing
  when `fit = "bemg"`. Note that unlike `tau` in the EGH model, this
  parameter has units of 1/time.

- `tau_left`: Exponential rate constant controlling left-sided tailing
  when `fit = "bemg"`. Note that unlike `tau` in the EGH model, this
  parameter has units of 1/time.

- `h`: The `height` parameter of the fitted BEMG shape function (only
  when `fit = "bemg"`). This is the value of the pure shape at `center`,
  before any baseline is added, and is distinct from `height`.

- `center`: The center parameter of the fitted BEMG shape function (only
  when `fit = "bemg"`). This is the location around which the
  exponential tails are applied. Because BEMG peaks can be asymmetric
  (`tau_right ≠ tau_left`) and a sloped baseline shifts the fitted
  curve, the actual peak maximum (`rt`) may differ from `center`.

- `FWHM`: The full-width at half maximum calculated as \\2.355 \sigma\\.

- `height`: The height of the peak above the fitted baseline.

- `area`: The area of the peak above the fitted baseline, estimated by
  trapezoidal approximation.

- `r.squared`: The coefficient of determination (\\R^2\\) of the fitted
  model to the raw data. (**Note**: this value is calculated by fitting
  a linear model of the fitted peak values to the raw data. This
  approach is statistically dubious, since the models are fit using
  non-linear least squares. Nevertheless, it can still be useful as a
  rough metric for "goodness-of-fit").

- `purity`: The peak purity as estimated by
  [`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md).

- `floor`: Constant baseline level (only when `baseline = "flat"`).

- `floor_start`: Baseline level at the left edge of the peak window
  (only when `baseline = "sloped"`).

- `floor_end`: Baseline level at the right edge of the peak window (only
  when `baseline = "sloped"`). Both `floor_start` and `floor_end` are
  constrained to be non-negative.

## Details

Peaks are detected by finding zero-crossings in the smoothed first
derivative of the specified chromatographic traces (using
[`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md)).
Peak models are then fit to the detected features using
[`fit_peaks`](https://ethanbass.github.io/chromatographR/reference/fit_peaks.md)
according to the value of `fit`. Available models include the
bidirectional exponentially modified gaussian (`"bemg"`),
exponential-gaussian hybrid (`"egh"`) or `"gaussian"`. Peak areas are
then calculated using a trapezoidal approximation.

Additional arguments can be provided to `find_peaks` to fine-tune the
peak-finding algorithm. For example, the `smooth_window` can be
increased to reduce splitting of broad peaks into multiple features.
Overly aggressive smoothing may cause small peaks to be overlooked.

The parameters `sd`, `FWHM`, `tau`, `tau_right`, `tau_left`, and `area`
are returned in units determined by `time_unit`, which defaults to
`minutes`. To compare directly with 'ChemStation' integration results,
the time units should be changed to seconds.

## Note

The bones of this function are adapted from the `getAllPeaks` function
(<https://github.com/rwehrens/alsace/blob/master/R/getAllPeaks.R>) by
Ron Wehrens (though the underlying algorithms for finding and
integrating peaks are not the same.

## References

- Lan, K. & Jorgenson, J. W. 2001. A hybrid of exponential and gaussian
  functions as a simple model of asymmetric chromatographic peaks.
  *Journal of Chromatography A* 915:1-13.
  [doi:10.1016/S0021-9673(01)00594-5](https://doi.org/10.1016/S0021-9673%2801%2900594-5)
  .

- Naish, P. J. & Hartwell, S. 1988. Exponentially Modified Gaussian
  functions — A good model for chromatographic peaks in isocratic HPLC?
  *Chromatographia*, 26: 285-296.
  [doi:10.1007/BF02268168](https://doi.org/10.1007/BF02268168) .

- Arase, Shuntaro, Kanta Horie, Takashi Kato, et al. 2016. Intelligent
  Peak Deconvolution through In-Depth Study of the Data Matrix from
  Liquid Chromatography Coupled with a Photo-Diode Array Detector
  Applied to Pharmaceutical Analysis. *Journal of Chromatography. A*,
  1469: 35–47.
  [doi:doi.org/10.1016/j.chroma.2016.09.037](https://doi.org/doi.org/10.1016/j.chroma.2016.09.037)
  .

- O'Haver, Tom. Pragmatic Introduction to Signal Processing:
  Applications in scientific measurement.
  <https://terpconnect.umd.edu/~toh/spectrum/> (Accessed January, 2022).

- Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
  LC–DAD using multivariate curve resolution: the alsace package for R.
  *Metabolomics* 11:143-154.
  [doi:10.1007/s11306-014-0683-5](https://doi.org/10.1007/s11306-014-0683-5)
  .

## See also

[`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md),
[`fit_peaks`](https://ethanbass.github.io/chromatographR/reference/fit_peaks.md),
[`print.peak_list`](https://ethanbass.github.io/chromatographR/reference/print.peak_list.md),
[`[.peak_list`](https://ethanbass.github.io/chromatographR/reference/sub-.peak_list.md)

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa_pr)
pks <- get_peaks(Sa_pr, lambdas = c('210'), sd_max = 50, fit = "egh")
}
```
