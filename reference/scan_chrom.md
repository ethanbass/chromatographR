# Plot spectra by clicking on the chromatogram

Plots a chromatographic trace from the specified chromatogram (`idx`),
at the specified wavelength (`lambda`) with a dotted red line to
indicate the user-selected retention time. The trace is a single column
from the chromatographic matrix.

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
  ...
)
```

## Arguments

- chrom_list:

  A list of chromatograms in matrix format (timepoints x wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object using the pointer in the `peak_table`.

- idx:

  Numerical index of chromatogram you wish to plot, or "max" to
  automatically select the chromatogram with the highest signal
  intensity at the specified peak or retention time.

- lambda:

  The wavelength used for plotting chromatographic traces and
  determining signal intensity.

- plot_spectrum:

  Logical. If `TRUE`, plots the spectrum of the chosen peak. Defaults to
  `TRUE`.

- peak_table:

  A `peak_table` object created by
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- scale_spectrum:

  Logical. If `TRUE`, scales spectrum to unit height. Defaults to
  `FALSE`.

- spectrum_labels:

  Logical. If `TRUE`, plots labels on maxima in spectral plot. Defaults
  to `TRUE`.

- export_spectrum:

  Logical. If `TRUE`, invisibly returns the spectrum as a `data.frame`.
  Defaults to `FALSE`.

- ...:

  Additional arguments to `plot_spectrum`.

## Value

Invisibly returns the numeric index of the scan selected by the user.

## Details

If `plot_spectrum` is `TRUE`, plots the spectrum for the specified
chromatogram at the user-specified retention time. The spectrum is a
single row from the chromatographic matrix.

## Side effects

If `engine == "base"`, plots are rendered to the active graphics device.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md)

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa_pr)
scan_chrom(Sa_pr, lambda = "210", idx = 2, export_spectrum = TRUE)
}
```
