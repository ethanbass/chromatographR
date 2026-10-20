# Attach reference spectra to a `peak_table` object

For each peak, selects one sample's spectrum as the reference: either
the spectrum with the highest mean correlation to that peak's spectra in
the other samples (`"max.cor"`), or the spectrum from the sample where
the peak is most intense (`"max.int"`). Spectra are scaled before they
are stored.

## Usage

``` r
attach_ref_spectra(peak_table, chrom_list, ref = c("max.cor", "max.int"))
```

## Arguments

- peak_table:

  A `peak_table` object created by
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- chrom_list:

  Optional list of chromatograms (timepoints × wavelengths). If omitted,
  the function attempts to retrieve the original `chrom_list` used to
  generate `peak_table`.

- ref:

  Method for selecting reference spectra: `"max.cor"` (maximum mean
  correlation, the default) or `"max.int"` (maximum intensity).

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
data(Sa_warp)
pk_tab <- attach_ref_spectra(pk_tab, ref = "max.int")
pk_tab <- attach_ref_spectra(pk_tab, ref = "max.cor")
```
