# Subset peak table

Subsets a `peak_table` while preserving and subsetting all associated
slots (`tab`, `pk_meta`, `sample_meta`, `instrument_meta`,
`ref_spectra`).

## Usage

``` r
# S3 method for class 'peak_table'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  A
  [`peak_table`](https://ethanbass.github.io/chromatographR/reference/peak_table-class.md)
  object.

- i:

  Row indices (samples) to subset.

- j:

  Column indices (peaks) to subset.

- ...:

  Additional arguments (ignored).

- drop:

  Passed to the `[ `indexing operator.

## Value

A `peak_table` object.

## Note

Subsetting a `peak_table` does not modify the associated `chrom_list`.
Methods that rely on the `chrom_list` (e.g.,
[`plot_spectrum`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md))
may not work as expected on a subsetted `peak_table` unless the
`chrom_list` is subsetted separately.

## Author

Ethan Bass
