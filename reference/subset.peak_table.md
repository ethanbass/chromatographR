# Subset peak table

Returns subset of `peak_table` object.

## Usage

``` r
# S3 method for class 'peak_table'
subset(x, subset, select, drop = FALSE, ...)
```

## Arguments

- x:

  A `peak_table` object.

- subset:

  Logical expression indicating rows (samples) to keep from
  `peak_table`; missing values are taken as false.

- select:

  Logical expression indicating columns (peaks) to select from
  `peak_table`.

- drop:

  Logical. Passed to indexing operator.

- ...:

  Additional arguments (placeholder).

## Value

A `peak_table` object with samples specified by `subset` and peaks
specified by `select`.

## Author

Ethan Bass
