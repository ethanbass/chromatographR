# Normalize peak table or chromatograms

Normalizes peak table or list of chromatograms by specified column in
sample metadata. Metadata must first be attached to the `peak_table`
using
[`attach_metadata`](https://ethanbass.github.io/chromatographR/dev/reference/attach_metadata.md).

## Usage

``` r
normalize_data(
  peak_table,
  column,
  chrom_list,
  what = c("peak_table", "chrom_list"),
  by = c("meta", "peak")
)
```

## Arguments

- peak_table:

  A `peak_table` object.

- column:

  The name of the column containing the weights.

- chrom_list:

  List of chromatograms for normalization. The samples must be in same
  order as the peak_table. If no argument is provided here, the function
  will try to find the `chrom_list` object used to create the provided
  `peak_table`.

- what:

  \`peak_table\` or list of chromatograms (\`chrom_list\`).

- by:

  Whether to normalize by a column in sample metadata (`meta`) or by a
  column in the peak table itself (`peak`).

## Value

A `peak_table` object where the peaks are normalized by the mass of each
sample.

## See also

[`get_peaktable`](https://ethanbass.github.io/chromatographR/dev/reference/get_peaktable.md)
[`attach_metadata`](https://ethanbass.github.io/chromatographR/dev/reference/attach_metadata.md)

## Author

Ethan Bass

## Examples

``` r
data(pk_tab)
path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
meta <- read.csv(path)
pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
norm <- normalize_data(pk_tab, "mass", what = "peak_table")
```
