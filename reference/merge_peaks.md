# Merge split peaks

Combines split peaks into a single column of the peak table.

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

  Method to merge peaks. Either `"max"` to select the largest peak from
  each sample or `"sum"` to add the peaks together.

## Value

A peak table similar to the input peak table, but where the specified
columns are combined.

## Details

In each sample, `method = "max"` keeps the largest of the specified
peaks and `method = "sum"` adds them. The result replaces the peak with
the highest mean across samples, whose metadata is retained; the other
peaks are removed.

## See also

Other utility functions:
[`combine_peaks()`](https://ethanbass.github.io/chromatographR/reference/combine_peaks.md),
[`reshape_chroms()`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md),
[`reshape_peaktable()`](https://ethanbass.github.io/chromatographR/reference/reshape_peaktable.md)

## Author

Ethan Bass

## Examples

``` r
data(pk_tab)
pk_tab <- merge_peaks(peak_table = pk_tab, peaks = c("V8","V9"))
```
