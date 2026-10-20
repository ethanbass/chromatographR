# Plot all spectra for chosen peak

Plot the spectra of a peak from several chromatograms. Wrapper for
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

  A list of chromatograms in matrix format (timepoints × wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object using the pointer in the `peak_table`.

- idx:

  Vector of chromatograms to plot, or `"all"` (the default).

- engine:

  Plotting backend to use. One of `"base"`,
  [`"ggplot"`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html),
  or [`"plotly"`](https://rdrr.io/pkg/plotly/man/plotly.html).
  `"ggplot2"` is also accepted here, for backward compatibility, but
  `"ggplot"` is the spelling the other plotting functions take.

- plot_spectrum:

  Logical. If `TRUE`, plots the spectrum of the chosen peak. Defaults to
  `TRUE`.

- export_spectrum:

  Logical. If `TRUE` (the default), invisibly returns the spectra as a
  `data.frame`.

- scale_spectrum:

  Logical. If `TRUE` (the default), scales each spectrum to the range
  \[0, 1\].

- overlapping:

  Logical. If `TRUE` (the default), plots the spectra on a single plot.
  Used only when `engine = "base"`.

- verbose:

  Logical. If `TRUE`, prints the index and name of each chromatogram
  plotted, with its retention time, wavelength and absorbance. Defaults
  to `FALSE`.

- what:

  What to look for. Either `"peak"` to extract spectral information for
  a certain peak, `"rt"` to scan by retention time, or `"idx"` to scan
  by numeric index. Defaults to "peak" mode.

- peak:

  The `peak` argument is deprecated. Please use `loc` instead.

- ...:

  Additional arguments to `plot_spectrum`.

## Value

- If `export_spectrum = TRUE` (the default), invisibly returns a
  `data.frame` with wavelengths as rows and one column per sample
  encoding absorbance at each wavelength (normalized if
  `scale_spectrum = TRUE`).

- If `export_spectrum = FALSE`, returns a `plotly` or `ggplot` object,
  or nothing if `engine == "base"`.

## Side effects

If `engine == "base"`, the spectra are rendered to the active graphics
device.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.chrom_list()`](https://ethanbass.github.io/chromatographR/reference/plot.chrom_list.md),
[`plot.matched_peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.matched_peak_list.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md),
[`plot_chrom_contour()`](https://ethanbass.github.io/chromatographR/reference/plot_chrom_contour.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_warp)
pks <- get_peaks(Sa_warp, lambda = "220")
pk_tab <- get_peaktable(pks)
plot_all_spectra(loc = "V11", peak_table = pk_tab, overlapping = TRUE)
```
