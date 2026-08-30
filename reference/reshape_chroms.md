# Reshape chromatograms

Converts a list of chromatographic matrices into a single long-format
data frame with one row per sample × retention time × wavelength
combination.

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

  Time resolution used when reshaping. Can be used to subsample the time
  axis. Defaults to full resolution.

- sample_var:

  Name of the new column containing sample identifiers

- lambdas:

  Wavelength(s) to include.

- rts:

  Retention times to include.

- transfer_metadata:

  Logical. If `TRUE`, metadata attributes are transferred to the output.
  Defaults to `FALSE`.

## Value

A data frame in long format with columns `rt`, `lambda`, `absorbance`,
and a sample identifier column specified by `sample_var`.

## Details

Each row corresponds to a single measurement of signal intensity at a
given retention time and wavelength for a specific sample.

## See also

Other utility functions:
[`combine_peaks()`](https://ethanbass.github.io/chromatographR/reference/combine_peaks.md),
[`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md),
[`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md),
[`merge_peaks()`](https://ethanbass.github.io/chromatographR/reference/merge_peaks.md),
[`reshape_peaktable()`](https://ethanbass.github.io/chromatographR/reference/reshape_peaktable.md)

## Author

Ethan Bass
