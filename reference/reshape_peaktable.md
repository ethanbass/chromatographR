# Reshape peaktable

Reshapes peak table from wide to long format

## Usage

``` r
reshape_peaktable(x, peaks, metadata, fixed_levels = TRUE)
```

## Arguments

- x:

  A `peak_table` object.

- peaks:

  A character vector specifying the peaks to include. If the character
  vector is named, the names of the vector elements will be used in
  place of the original peak names.

- metadata:

  A character vector specifying the metadata fields to include.

- fixed_levels:

  Logical. Whether to fix factor levels of features in the order
  provided. Defaults to `TRUE`.

## Value

A data.frame containing the information for the specified peaks in long
format.

## See also

Other utility functions:
[`combine_peaks()`](https://ethanbass.github.io/chromatographR/reference/combine_peaks.md),
[`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md),
[`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md),
[`merge_peaks()`](https://ethanbass.github.io/chromatographR/reference/merge_peaks.md),
[`reshape_chroms()`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md)

## Author

Ethan Bass
