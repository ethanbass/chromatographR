# Plot all spectra for chosen peak.

Plot multiple for a given peak in peak table. Wrapper for
[`plot_spectrum`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md).

## Usage

``` r
plot_all_spectra(
  peak,
  peak_table,
  chrom_list,
  idx = "all",
  chrs = NULL,
  engine = c("base", "ggplot2", "plotly"),
  plot_spectrum = TRUE,
  export_spectrum = TRUE,
  scale_spectrum = TRUE,
  overlapping = TRUE,
  verbose = FALSE,
  what = c("peak", "rt", "idx"),
  ...
)
```

## Arguments

- peak:

  The name of a peak to plot (in character format).

- peak_table:

  The peak table (output from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)
  function).

- chrom_list:

  A list of chromatograms in matrix format (timepoints x components). If
  no argument is provided here, the function will try to find the
  `chrom_list` object used to create the provided `peak_table`.

- idx:

  Vector of chromatograms to plot.

- chrs:

  Deprecated. Please use `idx` instead.

- engine:

  Which plotting engine to use: `base`, `ggplot2`, or `plotly`.

- plot_spectrum:

  Logical. If TRUE, plots the spectrum of the chosen peak.

- export_spectrum:

  Logical. If TRUE, exports spectrum to console. Defaults to FALSE.

- scale_spectrum:

  Logical. If TRUE, scales spectrum to unit height.

- overlapping:

  Logical. If TRUE, plot spectra in single plot.

- verbose:

  Logical. If TRUE, prints verbose output to console.

- what:

  What to look for. Either `peak` to extract spectral information for a
  certain peak, `rt` to scan by retention time, or `idx` to scan by
  numeric index. Defaults to "peak" mode.

- ...:

  Additional arguments to plot_spectrum.

## Value

If `export_spectrum` is TRUE, returns the spectra as a ` data.frame`
with wavelengths as rows and one column for each sample in the
`chrom_list` encoding the absorbance (or normalized absorbance, if
`scale_spectrum` is TRUE) at each wavelength. Otherwise, there is no
return value.

## Side effects

If `plot_spectrum` is TRUE, plots the spectra for the specified
chromatogram (`idx`) of the given `peak`. The spectrum is a single row
from the chromatographic matrix.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md)

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa_warp)
pks <- get_peaks(Sa_warp, lambda = "220")
pk_tab <- get_peaktable(pks)
plot_all_spectra(peak = "V13", peak_table = pk_tab, overlapping = TRUE)
}
```
