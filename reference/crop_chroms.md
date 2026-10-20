# Crop chromatograms

Restricts chromatograms to a range of retention times and/or wavelengths
by subsetting the existing grid.

## Usage

``` r
crop_chroms(x, rt_range = NULL, lambda_range = NULL)
```

## Arguments

- x:

  A numeric matrix (time × wavelength) or a list of such matrices. Row
  names must correspond to time points and column names to wavelengths.

- rt_range:

  Numeric vector of length 2 giving the inclusive range of retention
  times to keep, in the time unit of the data. Either element may be
  `NA` for an open bound. Defaults to `NULL`, which keeps every time
  point.

- lambda_range:

  Numeric vector of length 2 giving the inclusive range of wavelengths
  to keep. Either element may be `NA` for an open bound. Defaults to
  `NULL`, which keeps every wavelength.

## Value

The cropped data, in the same shape as `x`.

## Details

Cropping is pure subsetting: the returned data are a submatrix of the
input, so values and dimnames pass through verbatim. This is the
difference from cropping via
[`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md),
which reaches a narrower range by interpolating onto a new axis —
resampling the signal unless the requested grid happens to coincide with
the existing time points, and reformatting the retention times either
way.

Bounds are inclusive, and either end of a range may be `NA` to leave it
open. Ranges are expressed in the units of the data itself: retention
times in whatever unit the chromatogram carries, as returned by
[`get_times`](https://ethanbass.github.io/chromatographR/reference/get_times.md),
and wavelengths as labeled in the column names, as returned by
[`get_lambdas`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md).

An `"estimated_baseline"` attribute attached by
[`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md)
or
[`correct_baseline`](https://ethanbass.github.io/chromatographR/reference/correct_baseline.md)
is cropped alongside the data, so it stays aligned and can still be
retrieved with
[`get_baseline`](https://ethanbass.github.io/chromatographR/reference/get_baseline.md).

A range that does not overlap the data at all is an error. A range that
merely extends past one end is honored as far as the data allow, with a
warning reporting what was actually kept, since this is usually a typo
or a unit mismatch. Cropping chromatograms whose time axes are offset
from one another — as unaligned chromatograms generally are — can leave
them with different numbers of rows, which also warns: functions that
read the time axis from the first chromatogram would then mislabel the
rest, so align before cropping if the result needs a common time axis.

## See also

[`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md),
[`correct_baseline`](https://ethanbass.github.io/chromatographR/reference/correct_baseline.md),
[`get_times`](https://ethanbass.github.io/chromatographR/reference/get_times.md),
[`get_lambdas`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md)

Other preprocessing functions:
[`correct_baseline()`](https://ethanbass.github.io/chromatographR/reference/correct_baseline.md),
[`preprocess()`](https://ethanbass.github.io/chromatographR/reference/preprocess.md),
[`subtract_blanks()`](https://ethanbass.github.io/chromatographR/reference/subtract_blanks.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_pr)
cropped <- crop_chroms(Sa_pr, rt_range = c(12, 16))
trimmed <- crop_chroms(Sa_pr, lambda_range = c(220, NA))
```
