# Get retention times

Get retention times from a list of chromatograms or a `peak_table`
object.

## Usage

``` r
get_times(x, idx = 1)
```

## Arguments

- x:

  A list of chromatograms or `peak_table` object.

- idx:

  Index of chromatogram from which to extract times.

## Value

Numeric vector of retention times from the chromatogram specified by
`idx`.

## See also

Other utility functions:
[`combine_peaks()`](https://ethanbass.github.io/chromatographR/reference/combine_peaks.md),
[`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md),
[`merge_peaks()`](https://ethanbass.github.io/chromatographR/reference/merge_peaks.md),
[`reshape_chroms()`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md),
[`reshape_peaktable()`](https://ethanbass.github.io/chromatographR/reference/reshape_peaktable.md)
