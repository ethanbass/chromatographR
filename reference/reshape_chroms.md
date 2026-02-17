# Reshape chromatograms

Reshapes a list of chromatograms from wide to long format.

## Usage

``` r
reshape_chroms(x, idx, sample_var = "sample", lambdas = NULL, rts = NULL)
```

## Arguments

- x:

  A list of chromatographic matrices in wide format.

- idx:

  Indices of chromatograms to convert.

- sample_var:

  String with name of new column containing sample IDs.

- lambdas:

  Vector specifying wavelength(s) to include.

- rts:

  Vector specifying retention times to include.

## Value

A list of chromatographic matrices in long format.

## See also

Other utility functions:
[`combine_peaks()`](https://ethanbass.github.io/chromatographR/reference/combine_peaks.md),
[`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md),
[`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md),
[`merge_peaks()`](https://ethanbass.github.io/chromatographR/reference/merge_peaks.md),
[`reshape_peaktable()`](https://ethanbass.github.io/chromatographR/reference/reshape_peaktable.md)

## Author

Ethan Bass
