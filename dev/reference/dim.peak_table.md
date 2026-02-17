# Return dimensions of a `peak_table` object.

Returns the dimensions of a `peak_table`, where the first dimension is
the number of samples and the second dimension is the number of peaks.

## Usage

``` r
# S3 method for class 'peak_table'
dim(x)
```

## Arguments

- x:

  A
  [`peak_table`](https://ethanbass.github.io/chromatographR/dev/reference/peak_table-class.md)
  object.

## Value

Returns the number of rows and columns in `peak_table$tab`.
