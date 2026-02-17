# Reshape chromatograms

Reshapes a list of chromatograms from wide to long format.

## Usage

``` r
reshape_chroms(
  x,
  idx,
  time_resolution = NULL,
  sample_var = "sample",
  lambdas = NULL,
  rts = NULL,
  transfer_metadata = FALSE
)
```

## Arguments

- x:

  A list of chromatographic matrices in wide format.

- idx:

  Indices of chromatograms to convert.

- time_resolution:

  Time resolution for plot. This argument can be used to thin the time
  axis while reshaping. By default the time resoution is not altered.

- sample_var:

  String with name of new column containing sample IDs.

- lambdas:

  Vector specifying wavelength(s) to include.

- rts:

  Vector specifying retention times to include.

- transfer_metadata:

  Logical. Whether to transfer metadata attributes or not. Defaults to
  `FALSE`.

## Value

A list of chromatographic matrices in long format.

## See also

Other utility functions:
[`combine_peaks()`](https://ethanbass.github.io/chromatographR/dev/reference/combine_peaks.md),
[`get_lambdas()`](https://ethanbass.github.io/chromatographR/dev/reference/get_lambdas.md),
[`get_times()`](https://ethanbass.github.io/chromatographR/dev/reference/get_times.md),
[`merge_peaks()`](https://ethanbass.github.io/chromatographR/dev/reference/merge_peaks.md),
[`reshape_peaktable()`](https://ethanbass.github.io/chromatographR/dev/reference/reshape_peaktable.md)

## Author

Ethan Bass
