# Plot fitted peak shapes.

Visually assess integration accuracy by plotting fitted peaks over
trace.

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
  a = 0.5,
  color = NULL,
  cex.points = 0.5,
  numbers = FALSE,
  cex.font = 0.5,
  y.offset = 25,
  plot_purity = FALSE,
  res,
  index = NULL
)
```

## Arguments

- x:

  A `peak_list` object. Output from the `get_peaks` function.

- ...:

  Additional arguments to main plot function.

- chrom_list:

  List of chromatograms (retention time x wavelength matrices)

- idx:

  Index or name of chromatogram to be plotted.

- lambda:

  Wavelength for plotting.

- points:

  Logical. If TRUE, plot peak maxima. Defaults to FALSE.

- ticks:

  Logical. If TRUE, mark beginning and end of each peak. Defaults to
  FALSE.

- a:

  Alpha parameter controlling the transparency of fitted shapes.

- color:

  The color of the fitted shapes.

- cex.points:

  Size of points. Defaults to 0.5

- numbers:

  Whether to number peaks. Defaults to FALSE.

- cex.font:

  Font size if peaks are numbered. Defaults to 0.5.

- y.offset:

  Y offset for peak numbers. Defaults to 25.

- plot_purity:

  Whether to add visualization of peak purity.

- res:

  time resolution for peak fitting

- index:

  This argument is deprecated. Please use `idx` instead.

## Value

No return value, called for side effects.

## Side effects

Plots a chromatographic trace from the specified chromatogram (`chr`) at
the specified wavelength (`lambda`) with fitted peak shapes from the
provided `peak_list` drawn underneath the curve.

## See also

[`get_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/get_peaks.md)

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/dev/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/dev/reference/mirror_plot.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/dev/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_all_spectra.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_chroms.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_chroms_heatmap.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_spectrum.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/dev/reference/scan_chrom.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_warp)
pks <- get_peaks(chrom_list = Sa_warp[1], lambdas = 210)
plot(pks, points = TRUE, ticks = TRUE)
```
