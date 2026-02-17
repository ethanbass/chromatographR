# Attach experimental metadata

Attaches sample metadata to a `peak_table` object. Metadata should be
provided as a data.frame object. One of the columns in the supplied
metadata must match exactly the row names of the peak table.

## Usage

``` r
attach_metadata(peak_table, metadata, column)
```

## Arguments

- peak_table:

  A `peak_table` object.

- metadata:

  A `data.frame` containing the sample metadata.

- column:

  The name of the column in your `metadata` object containing the sample
  names. Sample names must match the row names of `peak_table$tab`.

## Value

A `peak_table` object with attached metadata in the ` $sample_meta`
slot.

## See also

[`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)
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
