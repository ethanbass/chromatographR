# Attach reference spectra

Gathers reference spectra and attaches them to a `peak_table` object.
Reference spectra are defined either as the spectrum with the highest
intensity (`"max.int"`) or as the spectrum with the highest average
correlation to the other spectra in the peak table (`"max.cor"`).

## Usage

``` r
attach_ref_spectra(peak_table, chrom_list, ref = c("max.cor", "max.int"))
```

## Arguments

- peak_table:

  Peak table from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/dev/reference/get_peaktable.md).

- chrom_list:

  A list of chromatograms in matrix format (timepoints x wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object used to create the provided `peak_table`.

- ref:

  What criterion to use to select reference spectra. Current options are
  maximum correlation (`"max.cor"`) or maximum signal intensity
  (`"max.int"`).

## Value

A `peak_table` object with reference spectra attached in the
`$ref_spectra` slot.

## See also

[`get_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/get_peaks.md)
[`get_peaktable`](https://ethanbass.github.io/chromatographR/dev/reference/get_peaktable.md)

## Author

Ethan Bass

## Examples

``` r
data(pk_tab)
pk_tab <- attach_ref_spectra(pk_tab, ref = "max.int")
pk_tab <- attach_ref_spectra(pk_tab, ref = "max.cor")
```
