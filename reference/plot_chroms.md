# Plot traces from list of chromatograms

Visualizes absorbance traces from a list of one- or multi-wavelength
chromatograms using base R, ggplot2, or plotly graphics. For
multi-wavelength chromatograms, one or more wavelengths can be selected
using `lambdas`. Traces are thinned to `time_resolution` before
plotting. Each chromatogram specified by `idx` is plotted as a separate
trace in its own color.

## Usage

``` r
plot_chroms(
  x,
  lambdas = NULL,
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
  what = NULL,
  ...
)
```

## Arguments

- x:

  A list of chromatograms in matrix format (timepoints × wavelengths) or
  a `peak_table` containing a pointer to a list of chromatograms
  accessible in the current environment.

- lambdas:

  A character or numeric vector specifying the wavelengths to plot.
  Required for multi-wavelength chromatograms and optional for
  one-dimensional ones, where a supplied wavelength must match the
  column name.

- idx:

  A vector representing the names or numerical indices of the
  chromatograms to plot. Defaults to all chromatograms.

- time_resolution:

  Desired temporal resolution for plotting (in minutes). Chromatograms
  are thinned to approximately this interval prior to plotting, which
  speeds up rendering of large datasets. Defaults to `0.01`.

- time_unit:

  Time units of the provided chromatograms: either `min`, `s`, or `ms`.
  If possible, units will be detected automatically from chromatogram
  metadata. If `time_unit` attribute is not present and no argument is
  provided, the time units will default to minutes.

- xlim:

  Range of x axis values. Unless `ylim` is given, the y axis is then
  fitted to the traces within this range.

- ylim:

  Range of y axis values.

- xlab:

  X label.

- ylab:

  Y label. Defaults to `"Absorbance"`.

- engine:

  Plotting backend to use. One of
  [`"base"`](https://rdrr.io/r/graphics/matplot.html),
  [`"ggplot"`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html),
  or [`"plotly"`](https://rdrr.io/pkg/plotly/man/plotly.html).

- linewidth:

  Line width. Defaults to `1`.

- show_legend:

  Logical. Whether to display a legend. Defaults to `FALSE`.

- legend_position:

  Position of legend. Defaults to `"topright"`. For `engine = "base"`, a
  top or bottom keyword spreads the entries over columns and, unless
  `ylim` is given, widens the y-range so the legend clears the traces.
  For `engine = "ggplot"`, any value accepted by the `legend.position`
  argument of
  [`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html);
  the default leaves ggplot2's own placement.

- title:

  Title for plot.

- what:

  Which element to plot when each chromatogram in `x` is itself a list,
  as returned by
  [`chromConverter::read_chroms()`](https://ethanbass.github.io/chromConverter/reference/read_chroms.html)
  with several `what` values. A single name (e.g. `"dad"`), or a vector
  of names for a nested element (e.g.
  `c("instrument", "PMP1, Pressure")`). A nested element can also be
  written as one string split at the first colon, as in
  `"instrument:PMP1, Pressure"`, unless that string is itself an element
  name. Defaults to `"dad"` for such lists; ignored when `x` is a list
  of matrices.

- ...:

  Additional arguments to
  [`matplot`](https://rdrr.io/r/graphics/matplot.html)
  (`engine = "base"`),
  [`ggplot2::geom_line()`](https://ggplot2.tidyverse.org/reference/geom_path.html)
  (`engine = "ggplot"`), or the `line` list of
  [`plotly::plot_ly()`](https://rdrr.io/pkg/plotly/man/plot_ly.html)
  (`engine = "plotly"`).

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
