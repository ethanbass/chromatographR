# Filter peak lists

Filters peaks from a `peak_list` object by height, area, width,
retention time, or fit diagnostics (`converged` and `window_limited`).

## Usage

``` r
filter_peaks(
  peak_list,
  min_height,
  min_area,
  min_sd,
  max_sd,
  min_rt,
  max_rt,
  min_fwhm,
  max_fwhm,
  converged,
  window_limited
)
```

## Arguments

- peak_list:

  A `peak_list` object, consisting of a nested list of peak tables,
  where the first level is the sample, and the second level is the
  spectral wavelength. Every sample × wavelength component is described
  by a `data.frame` where every row represents a peak, and the columns
  contain information on retention time, peak shape, height, and area.

- min_height:

  Minimum peak height. Peaks whose `height` is `NA`, such as failed fits
  kept by
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md),
  are dropped if `min_height` is supplied and retained otherwise.

- min_area:

  Minimum peak area. Peaks whose `area` is `NA` are treated the same
  way.

- min_sd:

  The `min_sd` argument is deprecated. Please use `min_fwhm` instead.

- max_sd:

  The `max_sd` argument is deprecated. Please use `max_fwhm` instead.

- min_rt:

  Minimum retention time.

- max_rt:

  Maximum retention time.

- min_fwhm:

  Minimum peak width, as a full width at half maximum, in the time units
  of the peak list. Peaks whose `FWHM` could not be measured (`NA`) are
  retained.

- max_fwhm:

  Maximum peak width, as a full width at half maximum, in the time units
  of the peak list. Peaks whose `FWHM` could not be measured (`NA`) are
  retained.

- converged:

  Logical. If `TRUE`, retains only peaks whose fit converged; if
  `FALSE`, retains only those that did not. A peak has `converged == 0`
  either because the optimizer stopped short of convergence or because
  no optimizer ran at all, the latter occurring when the integration
  window held fewer points than the model has parameters. Requires a
  `converged` column, which peak lists made with `fit = "raw"` or before
  chromatographR 0.9.0 lack. Peaks whose `converged` is `NA` are
  retained.

- window_limited:

  Logical. If `FALSE`, drops peaks whose fitted shape hit a limit set by
  the integration window or the sampling interval, meaning the window or
  the sampling rather than the data determined the reported shape; if
  `TRUE`, retains only those. Requires a `window_limited` column, which
  peak lists made with `fit = "raw"` or before chromatographR 0.9.0
  lack. Peaks whose `window_limited` is `NA` are retained.

## Value

A filtered `peak_list` object containing only peaks that satisfy the
specified criteria. The call is appended to the `filters` attribute, as
a line of code with numbers to 15 significant digits; `print` shows it,
and
[`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)
copies it into `args$filters`.

## Note

The structure of `filter_peaks` is adapted from the `filterPeaks`
function in Ron Wehrens' [alsace](https://github.com/rwehrens/alsace)
package.

## See also

[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md),
[`filter_peaktable`](https://ethanbass.github.io/chromatographR/reference/filter_peaktable.md)

## Author

Ron Wehrens, Ethan Bass

## Examples

``` r
data(Sa_warp)
pks <- get_peaks(Sa_warp[1], lambda = "210")
filter_peaks(pks, min_height = 100)
#> A peak_list with 1 samples and 1 wavelength(s) (210 nm)
#> Fit method: egh | Time unit: min
#> Filtered with: filter_peaks(min_height = 100)
#> Total peaks: 8 (mean 8 per sample)
#> Source: Sa_warp[1]
filter_peaks(pks, converged = TRUE, window_limited = FALSE)
#> A peak_list with 1 samples and 1 wavelength(s) (210 nm)
#> Fit method: egh | Time unit: min
#> Filtered with: filter_peaks(converged = TRUE, window_limited = FALSE)
#> Total peaks: 27 (mean 27 per sample)
#> Source: Sa_warp[1]
```
