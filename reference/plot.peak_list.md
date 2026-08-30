# Overlay fitted peak shapes on chromatograms

Visually assess peak integration accuracy by overlaying fitted peak
shapes over chromatographic traces.

## Usage

``` r
# S3 method for class 'peak_list'
plot(
  x,
  ...,
  chrom_list,
  idx = 1,
  lambda = NULL,
  points = FALSE,
  ticks = FALSE,
  alpha = 0.5,
  color = NULL,
  cex.points = 0.5,
  numbers = FALSE,
  cex.font = 0.5,
  y.offset = 25,
  plot_purity = FALSE,
  res
)
```

## Arguments

- x:

  A `peak_list` object. Output from the `get_peaks` function.

- ...:

  Additional arguments to main plot function.

- chrom_list:

  List of chromatograms (retention time x wavelength matrices). If
  missing, extracted from environment using the pointer in `x`.

- idx:

  Index or name of chromatogram to plot from `chrom_list`.

- lambda:

  Wavelength (column) to use for plotting.

- points:

  Logical. If `TRUE`, display peak apex locations as points. Defaults to
  `FALSE`.

- ticks:

  Logical. If `TRUE`, mark peak boundaries with tick marks. Defaults to
  `FALSE`.

- alpha:

  Transparency of fitted peak shapes. Defaults to `0.5`.

- color:

  Color used to fill fitted peak shapes. If `NULL`, a default color is
  chosen based on the fitted model type.

- cex.points:

  Size of points. Defaults to `0.5`.

- numbers:

  If `TRUE`, label peaks with numeric identifiers. Defaults to `FALSE`.

- cex.font:

  Font size if peaks are numbered. Defaults to `0.5`.

- y.offset:

  Y offset for peak numbers. Defaults to `25`.

- plot_purity:

  Logical. If `TRUE`, overlays peak purity traces based on peak
  boundaries. Defaults to `FALSE`.

- res:

  Time resolution for peak fitting. If missing, inferred from
  `chrom_list`.

## Value

No return value, called for side effects.

## Details

The appearance of fitted peaks depends on the `"fit"` attribute of `x`,
which may be `"gaussian"`, `"egh"`, `"bemg"`, or `"raw"`.

Peak lists are expected to contain columns such as `rt`, `height`,
`start`, `end`, and `sd`, with additional parameters depending on the
fit type.

Time units in `x` are used to rescale width parameters for plotting.

Peak rendering errors are silently ignored.

## Side effects

Plots a chromatographic trace from the specified chromatogram (`chr`) at
the specified wavelength (`lambda`) with fitted peak shapes from the
provided `peak_list` drawn underneath the curve.

## See also

[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md),
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
```
