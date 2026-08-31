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

## Examples

``` r
data(pk_tab)
reshape_peaktable(pk_tab, peaks = c("V10", "V15"))
#>   peak sample lambda    rt      area
#> 1  V10    119    210 12.47  2.637016
#> 2  V10    121    210 12.47  1.042586
#> 3  V10    122    210 12.47  1.417466
#> 4  V10    458    210 12.47  4.157266
#> 5  V15    119    210 13.76 36.585326
#> 6  V15    121    210 13.76 20.566841
#> 7  V15    122    210 13.76 31.606156
#> 8  V15    458    210 13.76 42.991983
```
