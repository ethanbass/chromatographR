# Filter peak table

Utility function to remove peaks from peak table, e.g., because their
intensity is too low. Currently one can filter on `mean`, `median`, or
maximum (`"max"`) peak intensity or retention time.

## Usage

``` r
filter_peaktable(
  peak_table,
  rts,
  min_rt,
  max_rt,
  min_value,
  lambda,
  what = c("median", "mean", "max"),
  tol = 0
)
```

## Arguments

- peak_table:

  A peak_table object from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- rts:

  Vector of retention times to include in the peak table.

- min_rt:

  Minimum retention time to include in the peak table.

- max_rt:

  Maximum retention time to include in the peak table.

- min_value:

  Minimal cutoff for summarized peak intensity.

- lambda:

  Component(s) to include in peak table (e.g. wavelengths if you are
  using HPLC-DAD/UV).

- what:

  Whether to summarize intensities using `mean`, `median`, or `max`.
  Defaults to `median`.

- tol:

  Tolerance for matching of retention times to `rts`.

## Value

A peak table similar to the input, with all columns removed from the
peak table that do not satisfy the specified criteria.

## See also

[`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md),
[`filter_peaks`](https://ethanbass.github.io/chromatographR/reference/filter_peaks.md)

## Author

Ethan Bass

## Examples

``` r
data(pk_tab)
pk_tab <- filter_peaktable(pk_tab, min_rt = 10, max_rt = 16)
```
