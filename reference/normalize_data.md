# Normalize peak table or chromatograms

Normalizes peak table or list of chromatograms by a column in the sample
metadata or in the peak table. For normalization by sample metadata, the
metadata must first be attached to the `peak_table` using
[`attach_metadata`](https://ethanbass.github.io/chromatographR/reference/attach_metadata.md).

## Usage

``` r
normalize_data(
  peak_table,
  column,
  chrom_list = NULL,
  what = c("peak_table", "chrom_list"),
  by = NULL,
  on_invalid = c("warn", "error", "silent")
)
```

## Arguments

- peak_table:

  A `peak_table` object.

- column:

  The name of the column to be used for normalization.

- chrom_list:

  List of chromatograms for normalization. The samples must be in same
  order as the `peak_table`. If omitted, the function will attempt to
  find it automatically using the pointer from the `peak_table`.

- what:

  Output type to return: either `"peak_table"` (default) or
  `"chrom_list"`.

- by:

  Whether to normalize by a column in sample metadata (`meta`) or by a
  column in the peak table (`peak`). By default, this parameter is
  inferred based on the `column` name.

- on_invalid:

  How to handle invalid normalization values (i.e. zero, negative, or
  `NA` values). One of `"warn"` (the default), `"silent"`, or `"error"`.
  Invalid values are replaced with `NA` unless `"error"` is chosen.

## Value

Either a normalized `peak_table` object or a normalized `chrom_list`,
depending on the value of `what`.

## See also

[`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)
[`attach_metadata`](https://ethanbass.github.io/chromatographR/reference/attach_metadata.md)

## Author

Ethan Bass

## Examples

``` r
data(pk_tab)
path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
meta <- read.csv(path)

# normalize by sample mass
pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
norm <- normalize_data(pk_tab, "mass", what = "peak_table")

# normalize by internal standard
norm <- normalize_data(pk_tab, column = "V16", by = "peak")
```
