# Correct peak positions according to a PTW warping model

Corrects retention time differences in `peak_list` using parametric time
warping as implemented in the
[`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html) package.

## Usage

``` r
correct_peaks(peak_list, mod_list, chrom_list, match_names = TRUE)
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

  A list of `ptw` models.

- chrom_list:

  List of chromatograms from which the `ptw` models are derived.

- match_names:

  Logical. Whether to actively match the names of the `peak_list` to the
  list of models (`mod_list`). Defaults to `TRUE`.

## Value

The input list of peak tables is returned with extra columns containing
the corrected retention times.

## Details

Once an appropriate warping model has been established, corrected
retention times can be predicted for each peak. These are stored in a
separate column in the list of peak tables.

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
if (FALSE) { # interactive()
data(Sa_pr)
warp <- correct_rt(chrom_list = Sa_pr, lambdas = 210, what = "models")
pks <- get_peaks(Sa_pr, lambda = "210")
correct_peaks(pks, mod_list = warp, chrom_list = Sa_pr)
}
```
