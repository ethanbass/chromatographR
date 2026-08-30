# Attach reference spectra to a `peak_table` object.

Reference spectra are selected either as the spectrum with maximum
intensity (`"max.int"`) or as the spectrum with the highest average
correlation to all other spectra in the peak table (`"max.cor"`).

## Usage

``` r
attach_ref_spectra(peak_table, chrom_list, ref = c("max.cor", "max.int"))
```

## Arguments

- peak_table:

  A `peak_table` object created by
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- chrom_list:

  Optional list of chromatograms (timepoints × wavelengths). If `NULL`,
  the function attempts to retrieve the original `chrom_list` used to
  generate `peak_table`.

- ref:

  Method for selecting reference spectra. Options are `"max.int"`
  (maximum intensity) or `"max.cor"` (maximum average correlation).

## Value

A `peak_table` object with reference spectra stored in the
`$ref_spectra` slot.

## See also

[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md),
[`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)

## Author

Ethan Bass

## Examples

``` r
data(pk_tab)
pk_tab <- attach_ref_spectra(pk_tab, ref = "max.int")
pk_tab <- attach_ref_spectra(pk_tab, ref = "max.cor")
```
