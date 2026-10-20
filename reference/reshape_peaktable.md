# Reshape peaktable

Reshapes a peak table to long format, or to a wide data.frame with the
selected peaks alongside sample metadata.

## Usage

``` r
reshape_peaktable(
  x,
  peaks,
  metadata,
  fixed_levels = TRUE,
  format = c("long", "wide")
)
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
  provided. Defaults to `TRUE`. Ignored when `format = "wide"`.

- format:

  Either `"long"` (the default), with one row per sample and peak, or
  `"wide"`, with one row per sample and one column per peak.

## Value

A data.frame containing the information for the specified peaks. In long
format it has columns `sample`, `peak`, `lambda`, `rt`, the response and
any metadata; in wide format, `sample`, any metadata, then one column
per peak.

## See also

Other utility functions:
[`combine_peaks()`](https://ethanbass.github.io/chromatographR/reference/combine_peaks.md),
[`merge_peaks()`](https://ethanbass.github.io/chromatographR/reference/merge_peaks.md),
[`reshape_chroms()`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md)

## Author

Ethan Bass

## Examples

``` r
data(pk_tab)
reshape_peaktable(pk_tab, peaks = c("V8", "V13"))
#>   peak sample lambda       rt      area
#> 1  V13    119    210 13.76774 37.290324
#> 2  V13    121    210 13.76774 20.284852
#> 3  V13    122    210 13.76774 32.383822
#> 4  V13    458    210 13.76774 43.451017
#> 5   V8    119    210 12.46511  2.909701
#> 6   V8    121    210 12.46511  1.080224
#> 7   V8    122    210 12.46511  1.164881
#> 8   V8    458    210 12.46511  4.284673
reshape_peaktable(pk_tab, peaks = c(a = "V8", b = "V13"), format = "wide")
#>   sample        a        b
#> 1    119 2.909701 37.29032
#> 2    121 1.080224 20.28485
#> 3    122 1.164881 32.38382
#> 4    458 4.284673 43.45102
```
