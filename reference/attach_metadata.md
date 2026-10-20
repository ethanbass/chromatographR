# Attach sample metadata to a `peak_table`

Matches rows of `metadata` to samples by name and stores them in the
`sample_meta` slot, in the same order as the rows of `peak_table$tab`.
Samples missing from `metadata` get `NA` values, with a warning. Columns
already in `sample_meta`, such as fields added to the chromatograms with
[`chromConverter::add_metadata()`](https://ethanbass.github.io/chromConverter/reference/add_metadata.html)
(chromConverter 0.11.0 or later), are kept unless `metadata` replaces
them. To replace `sample_meta` entirely, set
`peak_table$sample_meta <- NA` first.

## Usage

``` r
attach_metadata(peak_table, metadata, column)
```

## Arguments

- peak_table:

  A `peak_table` object created by
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- metadata:

  A `data.frame` (or tibble) of sample metadata.

- column:

  Name of the column in `metadata` containing sample identifiers, which
  must be unique and match the row names of `peak_table$tab`.

## Value

A `peak_table` object with metadata stored in the `sample_meta` slot.

## See also

[`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md),
[`normalize_data`](https://ethanbass.github.io/chromatographR/reference/normalize_data.md)

## Author

Ethan Bass

## Examples

``` r
data(pk_tab)
path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
meta <- read.csv(path)
pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
```
