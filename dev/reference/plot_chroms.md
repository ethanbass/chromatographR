# Plot traces from list of chromatograms.

Plots the specified traces from a list of chromatograms.

## Usage

``` r
plot_chroms(
  x,
  lambdas,
  idx,
  time_resolution = 0.01,
  time_unit = NULL,
  xlim = NULL,
  ylim = NULL,
  xlab = "",
  ylab = "Absorbance",
  engine = c("base", "ggplot", "plotly"),
  linewidth = 1,
  show_legend = FALSE,
  legend_position = "topright",
  title = "",
  ...
)
```

## Arguments

- x:

  A list of chromatograms in matrix format (timepoints x wavelengths).

- lambdas:

  A character or numeric vector specifying the wavelengths to plot. For
  one-dimensional chromatograms, this argument can be ignored.

- idx:

  A vector representing the names or numerical indices of the
  chromatograms to plot.

- time_resolution:

  Time resolution for plot in minutes. Defaults to `0.01`. Thinning the
  time axis dramatically improved speed when plotting many
  chromatograms.

- time_unit:

  Time units of the provided chromatograms. Units will be detected
  automatically if possible from chromatogram metadata. If `time_unit`
  attribute is not present, the time units will default to to `min`.

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

  Logical. Whether to display legend or not. Defaults to `FALSE`.

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
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/dev/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/dev/reference/mirror_plot.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/dev/reference/plot.peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/dev/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_all_spectra.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_chroms_heatmap.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_spectrum.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/dev/reference/scan_chrom.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_warp)
plot_chroms(Sa_warp, lambdas = 210)
```
