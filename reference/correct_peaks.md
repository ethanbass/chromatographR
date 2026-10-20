# Correct peak positions according to a warping model

Corrects retention time differences in `peak_list` using warping models
fitted by
[`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md)
with any of its algorithms.

## Usage

``` r
correct_peaks(
  peak_list,
  mod_list,
  chrom_list,
  match_names = TRUE,
  chrom_list_aligned = NULL
)
```

## Arguments

- peak_list:

  A `peak_list` object created by
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md),
  containing a nested list of peak tables where the first level is the
  sample, and the second level is the spectral wavelength. Every
  wavelength is described by a matrix where each row corresponds to a
  feature, and the columns contain information on that feature (e.g.,
  retention time, peak width (FWHM), height, area, etc.)

- mod_list:

  Warping models, as returned by
  [`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md)
  with `what = "models"`: a `ptw_list`, an `lptw_list`, or a `VPdtw`
  object. Corrected times are on the time axis of the chromatograms
  `correct_rt` returns, which for `"vpdtw"` extends beyond the reference
  at either end.

- chrom_list:

  List of chromatograms from which the models are derived. If omitted,
  it is retrieved using the name stored in `mod_list`.

- match_names:

  Logical. If `TRUE` (default), models are matched to the samples of
  `peak_list` by name when both are named; otherwise by position.

- chrom_list_aligned:

  Optional aligned chromatograms, as returned by
  [`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md)
  with `what = "corrected.values"`. Only its name is recorded, so the
  object must still exist under that name when it is used.

## Value

The input list of peak tables is returned with extra columns containing
the corrected retention times.

## Details

Once an appropriate warping model has been established, corrected
retention times can be predicted for each peak. These are stored in an
`rt.cor` column, alongside `start.cor` and `end.cor` when the peak
tables have `start` and `end` columns.

The corrected times refer to the aligned time axis, so functions that
later read spectra or traces at those times need the aligned
chromatograms. The returned peak list stores the warping models for this
purpose, and
[`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)
carries them into the peak table, where they are re-applied to the
unaligned chromatograms when needed. To use a particular aligned object
instead, pass it as `chrom_list_aligned`.

## Note

This function is adapted from the `correctPeaks` function in the alsace
package by Ron Wehrens:
<https://github.com/rwehrens/alsace/blob/master/R/correctPeaks.R>.

## See also

[`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md)

## Author

Ron Wehrens, Ethan Bass

## Examples

``` r
data(Sa_pr)
warp <- correct_rt(chrom_list = Sa_pr, lambdas = 210, what = "models")
pks <- get_peaks(Sa_pr, lambda = "210")
correct_peaks(pks, mod_list = warp, chrom_list = Sa_pr)
#> A peak_list with 4 samples and 1 wavelength(s) (210 nm)
#> Fit method: egh | Time unit: min
#> Total peaks: 115 (mean 29 per sample)
#> Source: Sa_pr
```
