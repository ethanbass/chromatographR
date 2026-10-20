# Get peak list

Detects chromatographic peaks, fits peak models, and extracts peak
parameters at one or more wavelengths.

## Usage

``` r
get_peaks(
  chrom_list,
  lambdas,
  fit = c("egh", "bemg", "gaussian", "raw"),
  windows = NULL,
  targets = NULL,
  tol = 0.1,
  select_by = c("height", "proximity", "none"),
  mod_list = NULL,
  sd_max = NULL,
  max_iter = 1000,
  time_unit = NULL,
  baseline = c("none", "flat", "sloped"),
  estimate_purity = FALSE,
  purity_args = list(),
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

  A list of chromatograms in matrix format (timepoints × wavelengths).

- lambdas:

  A character or numeric vector specifying the wavelengths to find peaks
  at. Omit for one-dimensional chromatograms, or when `targets` has a
  `lambda` column.

- fit:

  Peak model. Options are bidirectional exponentially modified Gaussian
  (`"bemg"`), exponential-Gaussian hybrid (`"egh"`), `"gaussian"`, or
  `"raw"`. `"raw"` integrates the signal directly by trapezoidal
  approximation without fitting a model. Defaults to `"egh"`, or to
  `"raw"` when `windows` is provided, which allows no other value. See
  *Peak models*.

- windows:

  Optional two-column matrix or `data.frame` of start and end retention
  times (or a single pair of times), giving manual integration windows
  applied to every chromatogram. If provided, peak detection is skipped
  and each window is integrated in every sample. Requires `fit = "raw"`.
  Cannot be combined with `targets`. Defaults to `NULL`. See *Manual
  integration windows*.

- targets:

  Target retention times, as a numeric vector or a table of targets (see
  [`targeted_analysis`](https://ethanbass.github.io/chromatographR/reference/targeted_analysis.md)).
  If provided, only peaks whose apex falls within `tol` of a target are
  fit, and
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)
  groups peaks by target instead of clustering them. A `lambda` column
  replaces `lambdas`, and a `tol` column replaces the `tol` argument,
  which must then be omitted. Defaults to `NULL`, which fits every
  detected peak.

- tol:

  Numeric. The maximum distance between a peak's apex and its target in
  `targets`, in minutes, regardless of the time unit of `chrom_list`.
  Give one value for all targets, or one per target. Defaults to `0.1`.
  Ignored if `targets` is `NULL`.

- select_by:

  Character. Which peak to keep when more than one falls in a target's
  window: `"height"` (the default) keeps the most intense, `"proximity"`
  the one closest to the target, and `"none"` all of them. Ignored if
  `targets` is `NULL`.

- mod_list:

  Optional warping models from
  [`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md)
  (with `what = "models"`). If provided, the times in `targets` are read
  on the aligned time axis, so they can be copied from the `rt` row of
  an aligned peak table. Peaks are still detected and fitted on the
  unaligned data; use
  [`correct_peaks`](https://ethanbass.github.io/chromatographR/reference/correct_peaks.md)
  afterwards to add corrected retention times. Requires `targets`.

- sd_max:

  Deprecated. Filter peaks by width with
  [`filter_peaks`](https://ethanbass.github.io/chromatographR/reference/filter_peaks.md)
  instead, using `min_fwhm` and `max_fwhm`. A standard deviation in
  timepoints, converted to a `FWHM` cutoff of `2.355 * sd_max * dt` in
  the units of `time_unit`, where `dt` is the chromatogram's time step;
  wider peaks are discarded. `NULL` disables width filtering. Ignored if
  `windows` is provided.

- max_iter:

  Maximum number of iterations for non-linear least squares in
  [`fit_peaks`](https://ethanbass.github.io/chromatographR/reference/fit_peaks.md).
  Defaults to `1000`.

- time_unit:

  Units for `sd`, `FWHM`, `tau`, `tau_right`, `tau_left`, and `area`.
  Options are minutes (`"min"`), seconds (`"s"`), or milliseconds
  (`"ms"`). Defaults to the time unit of `chrom_list` (minutes if none
  is recorded). `rt`, `start` and `end` are not converted: they stay in
  the units of the chromatogram's own time axis. ChemStation reports
  areas with time in seconds, so use `"s"` to compare with it.

- baseline:

  Whether to measure each peak above a local baseline: `"none"` (the
  default), `"flat"` (a constant offset per peak), or `"sloped"` (an
  offset that varies linearly across the peak). See *Local baselines*.

- estimate_purity:

  Logical. Whether to estimate purity. Defaults to `FALSE`. When `TRUE`,
  calls
  [`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md)
  for every peak in every sample, at a cost that grows with the number
  of wavelengths.

- purity_args:

  A named list of additional arguments passed to
  [`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md),
  e.g. `list(background_correction = "none", metric = "max")`. A
  `lambdas` element here selects the wavelengths used to compare
  spectra, not the wavelengths at which peaks are integrated.

- show_progress:

  Logical. Whether to show progress bar. Defaults to `TRUE` if
  [`pbapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  is installed.

- cl:

  Number of cores (integer) or a cluster object from
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to `2`. On Windows, only a cluster object runs in parallel; an integer
  runs serially.

- collapse:

  Logical. Whether to combine the peaks from every wavelength into a
  single `data.frame` per sample. Defaults to `FALSE`.

- time.units:

  The `time.units` argument is deprecated. Please use `time_unit`
  instead.

- sd.max:

  The `sd.max` argument is deprecated. Please use `sd_max` instead.

- max.iter:

  The `max.iter` argument is deprecated. Please use `max_iter` instead.

- ...:

  Additional arguments passed to
  [`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md),
  such as `smooth_window`, `smooth_type`, `min_prominence_snr`,
  `min_curvature` and `min_amplitude`. See *Choosing the smoothing
  window* in
  [`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md).

## Value

An S3 object of class `peak_list`: a nested list (sample, then
wavelength) of `data.frame`s, or a list of one `data.frame` per sample
when `collapse = TRUE`, with one row per peak and columns:

- `rt`: The retention time of the peak maximum. For the fitted models it
  generally falls between timepoints, so match it to the time axis with
  a tolerance rather than `==` or
  [`match()`](https://rdrr.io/r/base/match.html);
  [`get_spectra`](https://ethanbass.github.io/chromatographR/reference/get_spectra.md)
  with `what = "rt"` returns the spectrum at the nearest timepoint.

- `start`: The retention time where the peak is estimated to begin.

- `end`: The retention time where the peak is estimated to end.

- `sd`: The standard deviation (\\\sigma\\) of the Gaussian component of
  the fitted peak shape. This is a parameter of the model rather than
  the width of the peak: for `"egh"` and `"bemg"` the asymmetry of the
  peak is carried by `tau` or by `tau_left`/`tau_right`, so \\\sigma\\
  describes only the core of the peak and is not comparable across
  values of `fit`. Use `FWHM` to compare widths. `NA` when
  `fit = "raw"`, which fits no model; the width of the integration
  window is available as `end - start`.

- `tau`: Empirical skewness parameter (in units of time) determining
  peak asymmetry when `fit == "egh"`.

- `tau_right`: Time constant (in units of time) of the exponential tail
  on the right side of the peak when `fit = "bemg"`. A larger value
  gives a longer tail, and a value near 0 means no tail on that side.

- `tau_left`: Time constant (in units of time) of the exponential tail
  on the left side of the peak when `fit = "bemg"`. A larger value gives
  a longer tail, and a value near 0 means no tail on that side.

- `h`: The amplitude parameter of the fitted BEMG shape function (only
  when `fit = "bemg"`). It scales the shape but is not its value at
  `center`, so use the `height` column for the height of the peak.

- `center`: The center parameter of the fitted BEMG shape function (only
  when `fit = "bemg"`). This is the location around which the
  exponential tails are applied. Because BEMG peaks can be asymmetric
  (`tau_right != tau_left`) and a sloped baseline shifts the fitted
  curve, the actual peak maximum (`rt`) may differ from `center`.

- `FWHM`: The full width at half maximum of the peak. For `"gaussian"`
  this is \\2\sqrt{2\ln 2}\\\sigma \approx 2.355 \sigma\\; for `"egh"`
  it is the closed form \\\sqrt{(\tau\ln 2)^2 + 8\sigma^2\ln 2}\\
  implied by the model (Lan & Jorgenson, 2001); and for `"bemg"`, which
  has no closed form, it is measured numerically from the fitted shape.
  For `"raw"` it is measured directly from the signal above the local
  `baseline`, and is `NA` when the integration window does not fall
  below half the height of the peak on both sides, as for a peak fused
  with a neighbor. Unlike `sd`, `FWHM` is comparable across values of
  `fit`.

- `height`: The height of the peak above the fitted baseline.

- `area`: The area of the peak above the fitted baseline, estimated by
  trapezoidal approximation.

- `r.squared`: Coefficient of determination (\\R^2\\) of the fitted
  curve over the peak window, `1 - SS_res / SS_tot`. Negative when the
  fit is worse than a flat line at the mean. Not returned when
  `fit = "raw"`.

- `purity`: The peak purity as estimated by
  [`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md).
  `NA` unless `estimate_purity = TRUE`.

- `window_limited`: `1` if a fitted shape parameter hit a limit set by
  the integration window, or a `"gaussian"` or `"egh"` fit narrowed to
  one sampling interval at half height, meaning the window or the
  sampling rather than the data determined the reported shape, and `0`
  otherwise. Not returned when `fit = "raw"`.

- `converged`: `1` if the nonlinear least squares optimizer reported
  convergence for the final stage of the fit, and `0` otherwise. A `0`
  means the reported parameters are wherever the optimizer stopped
  rather than a converged fit, so treat the corresponding `height` and
  `area` with caution. It also covers the case where part or all of the
  fit was skipped because the integration window holds too few points to
  determine it; the peak then keeps the last step that ran, or its
  starting values (see *Fitting procedure*). `get_peaks` warns when more
  than 20% of peaks sit in such a window. Not returned when
  `fit = "raw"`.

- `floor`: Constant baseline level (only when `baseline = "flat"`).

- `floor_start`: Baseline level at the left edge of the peak window
  (only when `baseline = "sloped"`).

- `floor_end`: Baseline level at the right edge of the peak window (only
  when `baseline = "sloped"`). For the peak models, `floor`,
  `floor_start` and `floor_end` are constrained to be non-negative.

- `target`: The row of `targets` the peak was selected for (only when
  `targets` is supplied).

## Details

Peaks are detected by
[`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md)
as zero-crossings in the smoothed first derivative of each trace, then
fitted by
[`fit_peaks`](https://ethanbass.github.io/chromatographR/reference/fit_peaks.md)
according to `fit` (see *Peak models* and *Fitting procedure*). Peak
areas are the trapezoidal integral of the fitted curve (or, for
`fit = "raw"`, of the signal) above the local `baseline` (see *Local
baselines*). To integrate fixed retention-time ranges instead, see
*Manual integration windows*.

## Note

Adapted from `getAllPeaks` by Ron Wehrens
(<https://github.com/rwehrens/alsace/blob/master/R/getAllPeaks.R>);
however, the underlying algorithms for finding and integrating peaks
differ.

## Peak models

The three peak models are fitted to each peak by nonlinear least
squares; `"raw"` uses no model.

- `"gaussian"`: a symmetric Gaussian, with a center, a width (`sd`) and
  a height.

- `"egh"`: the exponential-Gaussian hybrid of Lan & Jorgenson (2001), a
  Gaussian whose width grows on one side of the apex. It adds one
  parameter, `tau`, which is positive for a tailing peak and negative
  for a fronting one.

- `"bemg"`: the average of two exponentially modified Gaussians, one
  with an exponential tail to the left and one to the right. Their time
  constants, `tau_left` and `tau_right`, set the two tails
  independently, so a peak can tail, front, or both; a value near 0
  means no tail on that side.

- `"raw"`: no model. The signal itself is integrated across the window,
  so `r.squared`, `window_limited` and `converged` are not reported.

`"gaussian"` cannot follow tailing or fronting; `"egh"` follows one or
the other, but not both in the same peak; `"bemg"` follows both. The
more flexible models have more parameters (3, 4 and 5 respectively for
the peak shape, plus one more for a `"flat"` baseline or two more for a
`"sloped"` baseline).

## Local baselines

A peak's height and area are measured over its integration window, which
runs from the valley before the peak to the valley after it, as found by
[`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md),
or spans the manual window if the `windows` argument is given.

- With `baseline = "none"`, `height` and `area` are measured from zero,
  so the chromatogram should already be baseline corrected (see
  [`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md)).

- With the fitted peak models (`"gaussian"`, `"egh"`, `"bemg"`), a
  `"flat"` or `"sloped"` baseline is fitted together with the peak shape
  by nonlinear least squares. It cannot be negative; a flat baseline is
  capped at the lowest intensity in the window, and each end of a sloped
  one is capped at the median intensity of the first or last 2% of the
  window. `height` and `area` are those of the fitted shape above the
  baseline, with the area integrated across the integration window only.

- With `fit = "raw"`, nothing is fitted: a `"flat"` baseline is a
  horizontal line at the lower of the two boundary intensities, and a
  `"sloped"` one a straight line between them. Where that line lies
  above the signal, it subtracts area, so the area can come out
  negative;
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)
  records a negative area as `0`.

With `"flat"` or `"sloped"`, each peak gets its own baseline. This
differs from drop-line integration in vendor software such as
ChemStation, where a vertical line from each valley meets a baseline
shared by the whole cluster of peaks.

## Fitting procedure

A fit is built up in successive steps: first the peak shape, then a
constant baseline, and finally its slope. To fit the full model, a
peak's integration window must have at least as many timepoints as there
are parameters in the model. A step is skipped if the window has too few
timepoints to fit all the required parameters, and discarded if the
optimizer cannot start from its initial values or, for `"bemg"`, if a
fitted tail rises towards the edge of the window instead of falling back
towards the baseline. Either way, the peak keeps the last step that
worked. A `"sloped"` baseline can therefore come out flat, or a baseline
can be left at its starting estimate rather than fitted, and if even the
shape fails, every parameter will keep its starting value. In each of
these cases, `converged` is set to `0`.

Peaks that can't be fitted at all are dropped, unless they come from
`targets`. Those are kept, with `height` and `area` set to `NA` and
`converged` to `0`, so that every target found has a row.

## Manual integration windows

Manual windows integrate every sample over the same retention-time
ranges, so each compound is measured the same way in every sample, and
every sample gets the same set of peak table columns. The raw signal in
each window is integrated whether or not a peak is detected there, so a
window with no peak gets an area near zero. With `baseline = "none"`,
two adjacent windows sharing a boundary time act as a vertical drop
line: their areas sum exactly to the area of their union, so a peak can
be split at a chosen time. The `"flat"` and `"sloped"` baselines are
drawn within each window separately, so their areas do not sum to that
of the union. To fit peak models at known retention times, use `targets`
instead.

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
  [doi:10.1016/j.chroma.2016.09.037](https://doi.org/10.1016/j.chroma.2016.09.037)
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
[`[.peak_list`](https://ethanbass.github.io/chromatographR/reference/sub-.peak_list.md),
[`targeted_analysis`](https://ethanbass.github.io/chromatographR/reference/targeted_analysis.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_pr)
pks <- get_peaks(Sa_pr, lambdas = c('210'), fit = "egh")
pks <- filter_peaks(pks, max_fwhm = 1)
```
