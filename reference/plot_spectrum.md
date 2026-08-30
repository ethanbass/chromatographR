# Plot spectrum from peak table

Visualizes chromatographic traces and/or UV spectra for a selected peak,
retention time, or scan index from a `peak_table` or list of
chromatograms. Spectra can also be exported as a `data.frame` when
`export_spectrum = TRUE`.

## Usage

``` r
plot_spectrum(
  loc = NULL,
  peak_table,
  chrom_list = NULL,
  idx = "max",
  lambda = "max",
  plot_spectrum = TRUE,
  plot_trace = TRUE,
  spectrum_labels = TRUE,
  scale_spectrum = FALSE,
  export_spectrum = FALSE,
  verbose = TRUE,
  what = c("peak", "rt", "idx", "click"),
  engine = c("base", "plotly", "ggplot2"),
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

  Numerical index of chromatogram you wish to plot, or "max" to
  automatically select the chromatogram with the highest signal
  intensity at the specified peak or retention time.

- lambda:

  The wavelength used for plotting chromatographic traces and
  determining signal intensity.

- plot_spectrum:

  Logical. If `TRUE`, plots the spectrum of the chosen peak. Defaults to
  `TRUE`.

- plot_trace:

  Logical. If `TRUE`, plots the trace of the chosen peak at lambda.
  Defaults to `TRUE`.

- spectrum_labels:

  Logical. If `TRUE`, plots labels on maxima in spectral plot. Defaults
  to `TRUE`.

- scale_spectrum:

  Logical. If `TRUE`, scales spectrum to unit height. Defaults to
  `FALSE`.

- export_spectrum:

  Logical. If `TRUE`, invisibly returns the spectrum as a `data.frame`.
  Defaults to `FALSE`.

- verbose:

  Logical. If `TRUE` (default), prints verbose output to console.

- what:

  What to look for. Either `peak` to extract spectral information for a
  certain peak, `rt` to scan by retention time, `idx` to scan by numeric
  index, or `click` to manually select retention time by clicking on the
  chromatogram. Defaults to `"peak"` mode.

- engine:

  Which plotting engine to use: `base`, `ggplot2`, or `plotly`.

- ...:

  Additional arguments.

## Value

- If `export_spectrum = FALSE` (default), returns a `plotly` or `ggplot`
  object according to the specified `engine`. No value is returned when
  `engine = "base"`.

- If `export_spectrum` is `TRUE`, invisibly returns the spectrum as a
  `data.frame` with wavelengths as rows and a single column encoding the
  absorbance (normalized, if `scale_spectrum = TRUE`).

## Details

This function can be used to confirm the identity of a peak or assess
whether a peak table column likely represents a single compound.
Retention times may also be selected interactively by clicking on a
chromatographic trace when `what = "click"`.

When `plot_trace` is `TRUE`, the chromatographic trace for the specified
chromatogram (`idx`) is plotted at wavelength `lambda`, with a dotted
red line indicating the selected retention time (`loc`). The trace
corresponds to a single column of the chromatographic matrix.

When `plot_spectrum` is `TRUE`, the UV spectrum at the specified
retention time is plotted. The spectrum corresponds to a single row of
the chromatographic matrix.

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
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md)

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa)
pks <- get_peaks(Sa, lambda = "220.00000")
pk_tab <- get_peaktable(pks)
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 1))
plot_spectrum(loc = "V10", peak_table = pk_tab, what = "peak")
par(oldpar)
}
```
