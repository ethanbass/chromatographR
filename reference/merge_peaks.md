# Merge split peaks

Utility function to combine split peaks into a single column of the peak
table.

## Usage

``` r
merge_peaks(peak_table, peaks, method = c("max", "sum"))
```

## Arguments

- peak_table:

  A `peak_table` object from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- peaks:

  A vector specifying the names or indices of peaks to be merged.

- method:

  Method to merge peaks. Either `max` to select the largest peak from
  each sample or `sum` to add the peaks together.

## Value

A peak table similar to the input peak table, but where the specified
columns are combined.

## Details

Merges the specified peaks in peak table, by selecting the largest value
from each column if `method` is `"max"`. If `method` is `"sum"`, peaks
will be merged by summing their values.

## See also

Other utility functions:
[`combine_peaks()`](https://ethanbass.github.io/chromatographR/reference/combine_peaks.md),
[`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md),
[`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md),
[`reshape_chroms()`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md),
[`reshape_peaktable()`](https://ethanbass.github.io/chromatographR/reference/reshape_peaktable.md)

## Author

Ethan Bass

## Examples

``` r
data(pk_tab)
pk_tab <- merge_peaks(peak_table = pk_tab, peaks = c("V10","V11"))
```
