# Overlay fitted peak shapes on chromatograms

Draws fitted peak shapes over chromatographic traces to assess
integration accuracy.

## Usage

``` r
# S3 method for class 'peak_list'
plot(
  x,
  ...,
  chrom_list,
  idx = 1,
  lambda = NULL,
  targets = NULL,
  tol = 0.1,
  pad = 0.5,
  points = FALSE,
  ticks = FALSE,
  alpha = 0.5,
  color = NULL,
  cex.points = 0.5,
  numbers = FALSE,
  cex.font = 0.5,
  digits = 2,
  y.offset = 0.03,
  plot_purity = c("none", "ratio", "curves", "curves_local"),
  purity_args = list(),
  res
)
```

## Arguments

- x:

  A `peak_list` object. Output from the `get_peaks` function.

- ...:

  Additional arguments to main plot function. The plot is titled with
  the sample, wavelength and peak model by default; pass `main` to
  change the title, or `main = ""` to omit it, and `xaxt = "n"` or
  `yaxt = "n"` to hide an axis.

- chrom_list:

  List of chromatograms (retention time × wavelength matrices). If
  missing, extracted from environment using the pointer in `x`.

- idx:

  Index or name of chromatogram to plot from `chrom_list`.

- lambda:

  Wavelength (column) to use for plotting. Defaults to the first
  wavelength in `x`.

- targets:

  Optional target retention times, as a numeric vector or a table of
  targets (see
  [`targeted_analysis`](https://ethanbass.github.io/chromatographR/reference/targeted_analysis.md)).
  If supplied, the plot is zoomed to the span of every peak whose apex
  falls within `tol` of a target, or, for a target with no such peak,
  whose boundaries (widened by `tol`) enclose it. The vertical axis is
  rescaled to the contents of that window, so that small peaks are not
  flattened against the axis. Defaults to `NULL`.

- tol:

  Numeric. Half-width of the window around each target, in minutes,
  regardless of the time unit of `chrom_list`. Either a single value
  applied to every target, or a vector giving a separate tolerance for
  each target. A `tol` column in `targets` takes precedence. Defaults to
  `0.1`. Ignored if `targets` is `NULL`.

- pad:

  Fraction of the selected peak span to add on either side when zooming.
  Defaults to `0.5`, which keeps neighboring peaks in view. An explicit
  `xlim` or `ylim` passed through `...` takes precedence over the zoom;
  if `xlim` is given, the vertical axis is rescaled to that window
  instead.

- points:

  Logical. If `TRUE`, display peak apex locations as points. Defaults to
  `FALSE`.

- ticks:

  Logical. If `TRUE`, mark peak boundaries with tick marks. Defaults to
  `FALSE`.

- alpha:

  Transparency of fitted peak shapes. Defaults to `0.5`.

- color:

  Color used to fill fitted peak shapes. If `NULL`, the color depends on
  the fitted model type.

- cex.points:

  Size of points. Defaults to `0.5`.

- numbers:

  Labels drawn above each peak. `TRUE` labels peaks with their index in
  the peak list, and `"rt"` labels them with their retention time, which
  can be passed back to `targets`. Defaults to `FALSE`.

- cex.font:

  Font size if peaks are numbered. Defaults to `0.5`.

- digits:

  Number of decimal places for retention time labels when
  `numbers = "rt"`. Defaults to `2`.

- y.offset:

  Vertical offset of peak numbers, as a fraction of the height of the
  plot region. Defaults to `0.03`.

- plot_purity:

  Whether and how to draw peak purity traces below the chromatogram.
  `"none"` (default) draws nothing. `"ratio"` plots the purity ratio as
  lollipops on a log scale from 0.01 to 100, clamping values beyond it:
  blue for pure (ratio \< 1), red for impure, with a dotted reference
  line at ratio = 1. `"curves"` draws spectral dissimilarity (1 - r^2)
  against the noise threshold on a fixed \[0, 1\] scale. The threshold
  is a gray dotted line, and the dissimilarity line is blue where it
  falls below the threshold (pure) and red where it exceeds it (impure).
  `"curves_local"` is the same, but rescales each peak to its own
  maximum. `TRUE` maps to `"curves"` and `FALSE` to `"none"`.

- purity_args:

  A named list of additional arguments controlling the purity
  calculation, e.g. `list(background_correction = "none")`. See
  [`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md).
  Arguments that do not apply to the per-peak purity trace drawn here
  (`metric`, `min_run` and `try`) are ignored, so the same list can be
  passed to
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)
  and to this function.

- res:

  Time resolution for peak fitting. If missing, inferred from
  `chrom_list`.

## Value

No return value, called for side effects.

## Details

The appearance of fitted peaks depends on the `"fit"` attribute of `x`,
which may be `"gaussian"`, `"egh"`, `"bemg"`, or `"raw"`.

Peak lists must contain columns `rt`, `height`, `start`, `end`, and
`sd`, with additional columns depending on the fit type.

Width parameters are rescaled using the time units of `x`.

A peak whose shape cannot be drawn is skipped with a warning.

## Side effects

Plots the trace of chromatogram `idx` at wavelength `lambda`, with the
fitted peak shapes from `x` filled in beneath it.

## See also

[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.chrom_list()`](https://ethanbass.github.io/chromatographR/reference/plot.chrom_list.md),
[`plot.matched_peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.matched_peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md),
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
pks <- get_peaks(chrom_list = Sa_warp[1], lambdas = 210)
plot(pks, points = TRUE, ticks = TRUE)

plot(pks, targets = 11.12, ticks = TRUE)

op <- par(mfrow = c(1, 2))
for (t in c(12.19, 14.41)) plot(pks, targets = t)

par(op)
```
