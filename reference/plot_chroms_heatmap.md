# Plot chromatograms as heatmap

Plots the traces from a list of chromatograms as a heatmap, with one row
per chromatogram and color encoding absorbance at a single wavelength.
Adapted from
[`VPdtw::plot.VPdtw`](https://rdrr.io/pkg/VPdtw/man/plot.VPdtw.html).

## Usage

``` r
plot_chroms_heatmap(
  x,
  idx = NULL,
  lambdas = NULL,
  engine = c("base", "ggplot", "plotly"),
  show_legend = TRUE,
  xlim = NULL,
  zlim = NULL,
  legend_position = "topright",
  title = "",
  show_ylabs = FALSE
)
```

## Arguments

- x:

  A list of chromatograms in matrix format (timepoints × wavelengths) or
  a `peak_table` containing a pointer to a list of chromatograms
  accessible in the current environment.

- idx:

  A vector of names or numerical indices of the chromatograms to plot.
  Defaults to all chromatograms.

- lambdas:

  A single wavelength to plot. Required for multi-wavelength
  chromatograms and optional for one-dimensional ones, where a supplied
  wavelength must match the column name.

- engine:

  Plotting backend to use. One of
  [`"base"`](https://rdrr.io/r/graphics/matplot.html),
  [`"ggplot"`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html),
  or [`"plotly"`](https://rdrr.io/pkg/plotly/man/plotly.html).

- show_legend:

  Logical. Whether to display a color bar. Defaults to `TRUE`.

- xlim:

  Range of retention times to plot.

- zlim:

  Range of intensities spanned by the color scale. Values outside it are
  drawn in the end colors. Defaults to the range within `xlim`. Pass the
  same `zlim` to several plots to put them on one scale.

- legend_position:

  Position of legend. Defaults to `"topright"`. Cannot be repositioned
  when using the `base` plotting engine.

- title:

  Title for plot.

- show_ylabs:

  Logical. Whether to label each row with its chromatogram name.
  Defaults to `FALSE`.

## Value

A `plotly` or `ggplot` object when `engine = "plotly"` or
`engine = "ggplot"`, respectively. No value is returned when
`engine = "base"`.

## Side effects

If `engine == "base"`, the plot is rendered to the active graphics
device.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.chrom_list()`](https://ethanbass.github.io/chromatographR/reference/plot.chrom_list.md),
[`plot.matched_peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.matched_peak_list.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md),
[`plot_chrom_contour()`](https://ethanbass.github.io/chromatographR/reference/plot_chrom_contour.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_warp)
plot_chroms_heatmap(Sa_warp, lambdas = 210)


# zoomed in, on the color scale of the whole run
zlim <- range(sapply(Sa_warp, function(m) m[, "210"]), na.rm = TRUE)
plot_chroms_heatmap(Sa_warp, lambdas = 210, xlim = c(13, 15), zlim = zlim)
```
