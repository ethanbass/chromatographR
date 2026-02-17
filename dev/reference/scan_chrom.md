# Plot spectra by clicking on the chromatogram.

Plot spectra by clicking on the chromatogram.

## Usage

``` r
scan_chrom(
  chrom_list,
  idx,
  lambda,
  plot_spectrum = TRUE,
  peak_table = NULL,
  scale_spectrum = FALSE,
  spectrum_labels = TRUE,
  export_spectrum = FALSE,
  chr = NULL,
  ...
)
```

## Arguments

- chrom_list:

  A list of chromatograms in matrix format (timepoints x wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object used to create the provided `peak_table`.

- idx:

  Numerical index of chromatogram you wish to plot.

- lambda:

  The wavelength to plot the trace at.

- plot_spectrum:

  Logical. Whether to plot the spectrum or not.

- peak_table:

  The peak table (output from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/dev/reference/get_peaktable.md)
  function).

- scale_spectrum:

  Logical. If TRUE, scales spectrum to unit height. Defaults to FALSE.

- spectrum_labels:

  Logical. If TRUE, plots labels on maxima in spectral plot. Defaults to
  TRUE.

- export_spectrum:

  Logical. If TRUE, exports spectrum to console. Defaults to FALSE.

- chr:

  Deprecated. Please use `idx` instead.

- ...:

  Additional arguments.

## Value

If `export_spectrum` is TRUE, returns the spectrum as a ` data.frame`
with wavelengths as rows and a single column encoding the absorbance (or
normalized absorbance, if `scale_spectrum` is TRUE) at each wavelength.
Otherwise, there is no return value.

## Side effects

Plots a chromatographic trace from the specified chromatogram (`idx`),
at the specified wavelength (`lambda`) with a dotted red line to
indicate the user-selected retention time. The trace is a single column
from the chromatographic matrix.

If `plot_spectrum` is TRUE, plots the spectrum for the specified
chromatogram at the user-specified retention time. The spectrum is a
single row from the chromatographic matrix.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/dev/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/dev/reference/mirror_plot.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/dev/reference/plot.peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/dev/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_all_spectra.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_chroms.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_chroms_heatmap.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_spectrum.md)

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa_pr)
scan_chrom(Sa_pr, lambda = "210", idx = 2, export_spectrum = TRUE)
}
```
