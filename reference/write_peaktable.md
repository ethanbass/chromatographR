# Export peak table

Exports peak table in `csv` or `xlsx` format.

## Usage

``` r
write_peaktable(
  peak_table,
  path,
  filename = "peak_table",
  format = c("csv", "xlsx"),
  what = c("tab", "pk_meta", "sample_meta", "instrument_meta", "ref_spectra", "args")
)
```

## Arguments

- peak_table:

  Peak table object from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- path:

  Path to write file.

- filename:

  File name. Defaults to `"peak_table"`.

- format:

  File format to export. Either `csv` or `xlsx`.

- what:

  Which elements of the `peak_table` to export. One or more of `"tab"`,
  `"pk_meta"`, `"sample_meta"`, `"instrument_meta"`, `"ref_spectra"`,
  and `"args"`.

## Value

No return value. The function is called for its side effects.

## Side effects

Exports peak_table object as `.csv` or `.xlsx` file according to the
value of the `format` argument.

## Examples

``` r
# \donttest{
data(pk_tab)
path_out = tempdir()
write_peaktable(pk_tab, path = path_out, what = c("tab"))
# }
```
