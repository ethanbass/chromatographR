# Filter peak table

Utility function to filter features (columns) in a `peak_table`.
Filtering is applied consistently across all components of the peak
table, including the intensity matrix (`tab`), feature metadata
(`pk_meta`), and reference spectra (if present). Filtering can be based
on retention time, wavelength, feature intensity summarized across
samples (mean, median, or max), and/or feature sparsity (i.e.,
proportion of zero values across samples).

## Usage

``` r
filter_peaktable(
  peak_table,
  rts,
  min_rt,
  max_rt,
  min_value,
  max_zeros,
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

  Vector of retention times used to select features. Values are mapped
  to the closest retention times in the peak table. If no match is found
  within `tol`, the value is ignored and a warning is issued.

- min_rt:

  Minimum retention time for features to be retained.

- max_rt:

  Maximum retention time for features to be retained.

- min_value:

  Minimum threshold for feature intensity summarized across samples
  (using the method specified by `what`).

- max_zeros:

  Maximum allowed feature sparsity (proportion of zero values across
  samples).

- lambda:

  Wavelength(s) used to select features in the peak table. Only features
  matching the specified wavelength(s) are retained.

- what:

  Method used to summarize feature intensity across samples for
  filtering. One of `"mean"`, `"median"` (default), or `"max"`.

- tol:

  Tolerance for matching of retention times to `rts`. A feature is
  considered a match if the absolute difference is ≤ `tol`.

## Value

A filtered `peak_table` containing only features (columns) that satisfy
the specified criteria. The same feature indices are applied to all
`peak_table` components.

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
