# Subtract blank chromatograms

Subtracts a blank chromatogram from each sample in a `chrom_list`, using
either the mean of all blanks, the nearest blank in run order, or the
most recent preceding blank.

## Usage

``` r
subtract_blanks(
  chrom_list,
  pattern = NULL,
  blank_id,
  method = c("mean", "nearest", "preceding"),
  zero_floor = TRUE
)
```

## Arguments

- chrom_list:

  A `chrom_list` object.

- pattern:

  A regular expression matched against `names(chrom_list)` to identify
  blanks (e.g. `"^Blank"` for samples named `Blank1`, `Blank2`, etc.).
  Takes precedence over `blank_id` if both are supplied.

- blank_id:

  Names or indices of the elements in `chrom_list` that are blanks.
  Ignored if `pattern` is supplied.

- method:

  How to select the blank to subtract for each sample: `"mean"`
  subtracts the elementwise mean of all blanks from every sample;
  `"nearest"` subtracts the blank closest in run order (list position);
  `"preceding"` subtracts the closest blank occurring at or before the
  sample in run order.

- zero_floor:

  Logical. If `TRUE` (default), values below zero after subtraction are
  set to zero.

## Value

A `chrom_list` with blank-subtracted samples. Blanks themselves are
removed from the returned list.
