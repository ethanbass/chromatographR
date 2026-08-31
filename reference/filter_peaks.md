# Filter peak lists

Utility function for filtering peaks from a `peak_list` object. Peaks
can be filtered according to height, area, standard deviation, and/or
retention time.

## Usage

``` r
filter_peaks(peak_list, min_height, min_area, min_sd, max_sd, min_rt, max_rt)
```

## Arguments

- peak_list:

  A `peak_list` object, consisting of a nested list of peak tables,
  where the first level is the sample, and the second level is the
  spectral component. Every component is described by a `data.frame`
  where every row represents a peak, and the columns contain information
  on retention time, peak shape, height, and area.

- min_height:

  Minimum peak height.

- min_area:

  Minimum peak area.

- min_sd:

  Minimal standard deviation.

- max_sd:

  Maximum standard deviation.

- min_rt:

  Minimum retention time.

- max_rt:

  Maximum retention time.

## Value

A filtered `peak_list` object containing only peaks that satisfy the
specified criteria.

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
#> Fit method: egh | Time unit: min | sd_max: 50
#> Total peaks: 8 (mean 8 per sample)
#> Source: Sa_warp[1]
```
