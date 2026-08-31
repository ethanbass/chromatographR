# Plot all spectra for chosen peak

Plot multiple for a given peak in peak table. Wrapper for
[`plot_spectrum`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md).

## Usage

``` r
plot_all_spectra(
  loc,
  peak_table,
  chrom_list,
  idx = "all",
  engine = c("base", "ggplot2", "plotly"),
  plot_spectrum = TRUE,
  export_spectrum = TRUE,
  scale_spectrum = TRUE,
  overlapping = TRUE,
  verbose = FALSE,
  what = c("peak", "rt", "idx"),
  peak = NULL,
  ...
)
```

## Arguments

- loc:

  The peak, retention time, or scan index from which to extract spectral
  data.

- peak_table:

  A `peak_table` object created by
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- chrom_list:

  A list of chromatograms in matrix format (timepoints x wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object using the pointer in the `peak_table`.

- idx:

  Vector of chromatograms to plot.

- engine:

  Which plotting engine to use: `base`, `ggplot2`, or `plotly`.

- plot_spectrum:

  Logical. If `TRUE`, plots the spectrum of the chosen peak. Defaults to
  `TRUE`.

- export_spectrum:

  Logical. If `TRUE`, invisibly returns the spectrum as a `data.frame`.
  Defaults to `FALSE`.

- scale_spectrum:

  Logical. If `TRUE`, scales spectrum to unit height. Defaults to
  `FALSE`.

- overlapping:

  Logical. If `TRUE`, plot spectra in single plot.

- verbose:

  Logical. If `TRUE`, prints verbose output to console. Defaults to
  `FALSE`.

- what:

  What to look for. Either `"peak"` to extract spectral information for
  a certain peak, `"rt"` to scan by retention time, or `"idx"` to scan
  by numeric index. Defaults to "peak" mode.

- peak:

  The name of a peak to plot (in character format).

- ...:

  Additional arguments to `plot_spectrum.`

## Value

- If `export_spectrum = FALSE` (default), a `plotly` or `ggplot` object,
  or nothing if `engine == "base"`.

- If `export_spectrum = TRUE`, invisibly returns a `data.frame` with
  wavelengths as rows and one column per sample encoding absorbance at
  each wavelength (normalized if `scale_spectrum = TRUE`).

## Side effects

If `engine == "base"`, the spectra are rendered to the active graphics
device.

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
