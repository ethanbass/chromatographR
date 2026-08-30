# Attaches sample metadata to a `peak_table` object by matching sample names.

Metadata is provided as a `data.frame`, with one column containing
sample identifiers matching the row names of `peak_table$tab`.

## Usage

``` r
attach_metadata(peak_table, metadata, column)
```

## Arguments

- peak_table:

  A `peak_table` object created by
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- metadata:

  A `data.frame` of sample metadata.

- column:

  Name of the column in `metadata` containing sample identifiers. Must
  match the row names of `peak_table$tab`.

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
