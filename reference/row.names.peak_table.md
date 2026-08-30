# Return row names from a `peak_table` object.

These will be the names of the samples.

## Usage

``` r
# S3 method for class 'peak_table'
row.names(x, ...)
```

## Arguments

- x:

  A
  [`peak_table`](https://ethanbass.github.io/chromatographR/reference/peak_table-class.md)
  object.

- ...:

  Additional arguments to
  [`row.names`](https://rdrr.io/r/base/row.names.html).

## Value

Returns the row names of `peak_table$tab`.
