# Preprocess time/wavelength data

Standard preprocessing of time × wavelength response matrices (e.g.
HPLC-DAD/UV data), consisting of: (i) optional subtraction of a
reference band of wavelengths, (ii) optional interpolation in either
dimension for dimensionality reduction, (iii) optional smoothing along
the spectral axis, and (iv) optional baseline correction along the time
axis. For densely sampled data (e.g., UV-VIS spectra), the size of the
matrix can be reduced by interpolation. By default, the data are
smoothed along the spectral axis using cubic smoothing splines
([`smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html)), then
baseline-corrected along the time axis by asymmetric least squares
([`asysm`](https://rdrr.io/pkg/ptw/man/asysm.html)). Three other
baseline estimators are available through `baseline`, and
`baseline = "none"` skips the step; see
[`baseline_estimators`](https://ethanbass.github.io/chromatographR/reference/baseline_estimators.md).

## Usage

``` r
preprocess(
  X,
  dim1,
  dim2,
  remove.time.baseline = NULL,
  baseline = c("als", "arpls", "aspls", "snip", "none"),
  attach_baseline = FALSE,
  zero_floor = FALSE,
  smooth_spectra = TRUE,
  reference_wavelengths = NULL,
  interpolate_rows = TRUE,
  interpolate_cols = TRUE,
  cl = 2,
  show_progress = NULL,
  outlier_cutoff = 5/60,
  ...,
  maxI = NULL,
  spec.smooth = NULL
)
```

## Arguments

- X:

  A numeric matrix (time × wavelength) or a list of such matrices. Row
  names must correspond to time points and column names to wavelengths.

- dim1:

  Numeric vector specifying the target time grid, usually built with
  [`seq`](https://rdrr.io/r/base/seq.html). The data will be
  interpolated onto these time points, which must be strictly
  increasing. Chromatograms that do not cover the whole of `dim1` are
  dropped with a warning. Writing `c(from, to, by = x)` in place of
  `seq(from, to, by = x)` yields a three-element vector rather than a
  sequence, so it raises an error. If `dim1` is omitted, the grid runs
  in steps of `0.01` (in the time units of the data) across the range
  shared by all chromatograms.

- dim2:

  Numeric vector specifying the target wavelength grid, usually built
  with [`seq`](https://rdrr.io/r/base/seq.html). The data will be
  interpolated onto these wavelengths, which must be strictly increasing
  and lie within the range of the original wavelengths; otherwise an
  error is raised. Ignored when `interpolate_cols = FALSE`. If omitted,
  the existing wavelengths are kept.

- remove.time.baseline:

  The `remove.time.baseline` argument is deprecated. Please use
  `baseline = "none"` instead.

- baseline:

  Baseline estimator to use along the time axis. Either asymmetric least
  squares (`"als"`, the default, via
  [`asysm`](https://rdrr.io/pkg/ptw/man/asysm.html)), asymmetrically
  reweighted penalized least squares (`"arpls"`), adaptive smoothness
  parameter penalized least squares (`"aspls"`), statistics-sensitive
  non-linear iterative peak-clipping (`"snip"`), or `"none"` to skip
  baseline correction along the time axis altogether. See
  [`baseline_estimators`](https://ethanbass.github.io/chromatographR/reference/baseline_estimators.md)
  for a description of each.

- attach_baseline:

  Logical. Whether to retain the estimated baseline alongside the
  corrected data as an `"estimated_baseline"` attribute. The parameter
  that produced it (`stiffness`, or `m` for `"snip"`) is always
  recorded, as a `"baseline_param"` attribute, whether or not the
  baseline itself is kept. The baseline is subtracted either way; this
  argument only controls whether the estimate is kept for inspection.
  Defaults to `FALSE`: the baseline is estimated separately for each
  wavelength, so it is the same size as the data and retaining it
  roughly doubles memory use. Adding the attribute back to the returned
  data recovers the uncorrected matrix only when `zero_floor = FALSE`
  and `maxI` and `reference_wavelengths` are unset: clipping alters the
  data without altering the stored estimate, `maxI` rescales both by a
  factor derived from the clipped data, and the subtracted reference
  band is not stored.

- zero_floor:

  Logical. Whether to clip negative values in the corrected data to
  zero. Defaults to `FALSE`. Clipping discards the half of the noise
  that falls below zero, which hides the negative dips that show where a
  baseline was over-subtracted, so it is turned off by default.

- smooth_spectra:

  Logical, indicating whether smoothing should be done along the
  spectral axis, according to
  [`smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html). Skipped
  for data with fewer than four wavelengths. Defaults to `TRUE`.

- reference_wavelengths:

  Numeric vector of length 2 giving a band of wavelengths,
  `c(from, to)`, whose mean is subtracted from every wavelength at each
  time point, like the reference wavelength of a diode-array detector.
  This removes offsets shared by all wavelengths, such as lamp drift or
  refractive-index effects, before the baseline is estimated. Choose a
  band where nothing in the samples absorbs. Either bound may be `NA` to
  leave that side open. Defaults to `NULL` (no subtraction).

- interpolate_rows:

  Logical. Whether to interpolate along the time axis (`dim1`). Defaults
  to `TRUE`.

- interpolate_cols:

  Logical. Whether to interpolate along the spectral axis (`dim2`).
  Defaults to `TRUE`.

- cl:

  Either an integer specifying the number of cores to use for parallel
  processing or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to `2`. On Windows integer values will be ignored.

- show_progress:

  Logical. Whether to show progress bar. Defaults to `TRUE` if
  [`pbapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  is installed.

- outlier_cutoff:

  Threshold (in minutes) for excluding chromatograms that end earlier
  than expected. Samples ending more than this value before the median
  end time are removed with a warning. The value is always in minutes,
  whatever time unit the chromatograms are stored in. Defaults to
  `5/60`, i.e. 5 seconds. Only applies when `dim1` is not specified.

- ...:

  Additional arguments to the selected `baseline` estimator: e.g.
  `stiffness` and `p` for `"als"`
  ([`asysm`](https://rdrr.io/pkg/ptw/man/asysm.html)), `stiffness` for
  `"arpls"`, `stiffness` and `k` for `"aspls"`, or `m` (required) for
  `"snip"`. Because `maxI` follows `...` in the signature, it must be
  named in full: an abbreviation such as `max` is passed to the
  estimator, which rejects it.

- maxI:

  If supplied, each chromatogram is rescaled so that its maximum
  intensity equals `maxI`.

- spec.smooth:

  The `spec.smooth` argument is deprecated. Please use `smooth_spectra`
  instead.

## Value

The function returns the preprocessed data matrix (or list of matrices),
with row names and column names indicating the time points and
wavelengths, respectively.

## Details

The operations are applied in the order given above: reference
subtraction first, on the original wavelengths, then interpolation, then
spectral smoothing, then baseline correction. Because the baseline is
estimated after the time axis is interpolated, `stiffness` is tied to
whichever `dim1` grid was chosen: the penalty applies on the index axis,
so the same value gives a different baseline on a differently sampled
grid. Choose `dim1` first, then tune `stiffness` against that grid.

## Note

Adapted from the `preprocess` function in the alsace package by Ron
Wehrens:
<https://github.com/rwehrens/alsace/blob/master/R/preprocess.R>.

## References

- Baek, S.-J., Park, A., Ahn, Y.-J., and Choo, J. 2015. Baseline
  correction using asymmetrically reweighted penalized least squares
  smoothing. *Analyst* 140:250-257.
  [doi:10.1039/C4AN01061B](https://doi.org/10.1039/C4AN01061B) .

- Zhang, F., Tang, X., Tong, A., Wang, B., Wang, J., Lv, Y., Tang, C.,
  and Wang, J. 2020. Baseline correction for infrared spectra using
  adaptive smoothness parameter penalized least squares method.
  *Spectroscopy Letters* 53:222-233.
  [doi:10.1080/00387010.2020.1730908](https://doi.org/10.1080/00387010.2020.1730908)
  .

- Ryan, C.G., Clayton, E., Griffin, W.L., Sie, S.H., and Cousens,
  D.R. 1988. SNIP, a statistics-sensitive background treatment for the
  quantitative analysis of PIXE spectra in geoscience applications.
  *Nuclear Instruments and Methods in Physics Research B* 34:396-402.
  [doi:10.1016/0168-583X(88)90063-8](https://doi.org/10.1016/0168-583X%2888%2990063-8)
  .

- Erb, D. 2025. *pybaselines: A Python library of algorithms for the
  baseline correction of experimental data*.
  [doi:10.5281/zenodo.5608581](https://doi.org/10.5281/zenodo.5608581) .

## See also

[`baseline_estimators`](https://ethanbass.github.io/chromatographR/reference/baseline_estimators.md),
[`correct_baseline`](https://ethanbass.github.io/chromatographR/reference/correct_baseline.md),
[`crop_chroms`](https://ethanbass.github.io/chromatographR/reference/crop_chroms.md),
[`get_baseline`](https://ethanbass.github.io/chromatographR/reference/get_baseline.md)

Other preprocessing functions:
[`correct_baseline()`](https://ethanbass.github.io/chromatographR/reference/correct_baseline.md),
[`crop_chroms()`](https://ethanbass.github.io/chromatographR/reference/crop_chroms.md),
[`subtract_blanks()`](https://ethanbass.github.io/chromatographR/reference/subtract_blanks.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa)
new.ts <- seq(10,18.66,by=.01) # choose time-points
new.lambdas <- seq(200, 318, by = 2) # choose wavelengths
Sa_pr <- preprocess(Sa[[1]], dim1 = new.ts, dim2 = new.lambdas)
```
