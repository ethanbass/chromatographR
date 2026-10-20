# Correct baseline along the time axis

Estimates a baseline along the time axis of each chromatogram and
subtracts it, using any of the estimators described in
[`baseline_estimators`](https://ethanbass.github.io/chromatographR/reference/baseline_estimators.md).

## Usage

``` r
correct_baseline(
  x,
  baseline = c("als", "arpls", "aspls", "snip"),
  attach_baseline = FALSE,
  zero_floor = FALSE,
  cl = 2,
  show_progress = NULL,
  ...
)
```

## Arguments

- x:

  A numeric matrix (time × wavelength), a list of such matrices, or a
  plain numeric vector for a single trace. Row names (or, for a vector,
  names) must correspond to time points and column names to wavelengths.

- baseline:

  Baseline estimator to use along the time axis. Either asymmetric least
  squares (`"als"`, the default, via
  [`asysm`](https://rdrr.io/pkg/ptw/man/asysm.html)), asymmetrically
  reweighted penalized least squares (`"arpls"`), adaptive smoothness
  parameter penalized least squares (`"aspls"`), or statistics-sensitive
  non-linear iterative peak-clipping (`"snip"`). See
  [`baseline_estimators`](https://ethanbass.github.io/chromatographR/reference/baseline_estimators.md)
  for a description of each.

- attach_baseline:

  Logical. Whether to retain the estimated baseline alongside the
  corrected data, as an `"estimated_baseline"` attribute on each
  returned matrix, retrievable with
  [`get_baseline`](https://ethanbass.github.io/chromatographR/reference/get_baseline.md).
  The baseline is subtracted either way; this argument only controls
  whether the estimate is kept for inspection. Defaults to `FALSE`: the
  baseline is estimated separately for each wavelength, so it is the
  same size as the data and retaining it doubles memory use.

- zero_floor:

  Logical. Whether to clip negative values in the corrected data to
  zero. Defaults to `FALSE`, because clipping also removes the negative
  dips that show where a baseline was over-subtracted. It can be useful
  with `baseline = "als"`, which can overshoot the signal; `"snip"`
  cannot overshoot by construction, so clipping is never needed there.
  The `zero_floor` argument of
  [`subtract_blanks`](https://ethanbass.github.io/chromatographR/reference/subtract_blanks.md),
  by contrast, defaults to `TRUE`, since a negative value there means
  the blank exceeded the sample rather than that an estimator overshot.

- cl:

  Either an integer specifying the number of cores to use for parallel
  processing or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to `2`. On Windows, integer values are ignored.

- show_progress:

  Logical. Whether to show a progress bar. Defaults to `TRUE` if
  [`pbapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  is installed.

- ...:

  Additional arguments to the selected `baseline` estimator: e.g.
  `stiffness` and `p` for `"als"`
  ([`asysm`](https://rdrr.io/pkg/ptw/man/asysm.html)), `stiffness` for
  `"arpls"`, `stiffness` and `k` for `"aspls"`, or `m` for `"snip"`.

## Value

The baseline-corrected data, in the same shape as `x`, with row and
column names (or, for a vector, names) unchanged.

## Details

Unlike
[`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md),
this function only subtracts the baseline; it does not interpolate,
smooth, or rescale.

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
[`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md),
[`get_baseline`](https://ethanbass.github.io/chromatographR/reference/get_baseline.md),
[`subtract_blanks`](https://ethanbass.github.io/chromatographR/reference/subtract_blanks.md)

Other preprocessing functions:
[`crop_chroms()`](https://ethanbass.github.io/chromatographR/reference/crop_chroms.md),
[`preprocess()`](https://ethanbass.github.io/chromatographR/reference/preprocess.md),
[`subtract_blanks()`](https://ethanbass.github.io/chromatographR/reference/subtract_blanks.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_pr)
bc <- correct_baseline(Sa_pr[1], cl = 1, show_progress = FALSE)
```
