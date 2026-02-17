# Return first or last parts of a `peak_table`.

Returns the first or last parts of the
[`peak_table`](https://ethanbass.github.io/chromatographR/dev/reference/peak_table-class.md).

## Usage

``` r
# S3 method for class 'peak_table'
head(x, ...)

# S3 method for class 'peak_table'
tail(x, ...)
```

## Arguments

- x:

  A
  [`peak_table`](https://ethanbass.github.io/chromatographR/dev/reference/peak_table-class.md)
  object.

- ...:

  Additional arguments to [`tail`](https://rdrr.io/r/utils/head.html).

## Value

The first or last `n` rows of `peak_table$tab`.
