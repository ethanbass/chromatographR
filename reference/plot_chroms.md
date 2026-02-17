# Plot traces from list of chromatograms.

Plots the specified traces from a list of chromatograms.

## Usage

``` r
plot_chroms(
  x,
  lambdas,
  idx,
  xlim = NULL,
  ylim = NULL,
  xlab = "",
  ylab = "Absorbance",
  engine = c("base", "ggplot", "plotly"),
  linewidth = 1,
  show_legend = TRUE,
  legend_position = "topright",
  title = "",
  ...
)
```

## Arguments

- x:

  A list of chromatograms in matrix format (timepoints x wavelengths).

- lambdas:

  A character or numeric vector specifying the wavelengths to plot.

- idx:

  A vector representing the names or numerical indices of the
  chromatograms to plot.

- xlim:

  Range of x axis.

- ylim:

  Range of y axis.

- xlab:

  X label.

- ylab:

  Y label. Defaults to "Absorbance".

- engine:

  Plotting engine. Either `base`
  ([`matplot`](https://rdrr.io/r/graphics/matplot.html)),
  [`plotly`](https://rdrr.io/pkg/plotly/man/plotly.html), or
  [ggplot](https://ggplot2.tidyverse.org/reference/ggplot2-package.html).

- linewidth:

  Line width.

- show_legend:

  Logical. Whether to display legend or not. Defaults to TRUE.

- legend_position:

  Position of legend.

- title:

  Title for plot.

- ...:

  Additional arguments to plotting function specified by `engine`.

## Value

No return value, called for side effects.

## Side effects

Plots the traces of the specified chromatograms `idx` at the specified
wavelengths `lambdas`. Plots can be produced using base graphics,
ggplot2, or plotly, according to the value of `engine`.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_warp)
plot_chroms(Sa_warp, lambdas = 210)
```
