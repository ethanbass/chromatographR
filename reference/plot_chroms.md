# Plot traces from list of chromatograms

Visualizes absorbance traces from a list of one- or multi-wavelength
chromatograms using base R, ggplot2, or plotly graphics. For
multi-wavelength chromatograms, one or more wavelengths can be selected
using `lambdas`. When plotting many chromatograms, traces are
automatically thinned according to `time_resolution` to improve
rendering performance. Each chromatogram specified by `idx` is plotted
as a separate trace and assigned a distinct color.

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

  A list of chromatograms in matrix format (timepoints × wavelengths) or
  a `peak_table` containing a pointer to a list of chromatograms
  accessible in the current environment.

- lambdas:

  A character or numeric vector specifying the wavelengths to plot.
  Ignored for one-dimensional chromatograms.

- idx:

  A vector representing the names or numerical indices of the
  chromatograms to plot.

- time_resolution:

  Desired temporal resolution for plotting (in minutes). Chromatograms
  are thinned to approximately this interval prior to plotting, which
  can substantially improve performance for large datasets. Defaults to
  `0.01`.

- time_unit:

  Time units of the provided chromatograms: either `min`, `s`, or `ms`.
  If possible, units will be detected automatically from chromatogram
  metadata. If `time_unit` attribute is not present and no argument is
  provided, the time units will default to minutes.

- xlim:

  Range of x axis values.

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

  Logical. Whether to display legend or not. Defaults to `FALSE`.

- legend_position:

  Position of legend. Defaults to `"topright"`

- title:

  Title for plot.

- ...:

  Additional arguments to plotting function specified by `engine`.

## Value

A `plotly` or `ggplot` object when `engine = "plotly"` or
`engine = "ggplot"`, respectively. No value is returned when
`engine = "base"`.

## Details

For multi-wavelength chromatograms, wavelength indices are resolved
internally using `get_lambda_idx()`. Time values are converted according
to `time_unit`, either supplied explicitly or inferred from chromatogram
metadata.

When using the `"ggplot"` or `"plotly"` engines, chromatograms are first
reshaped into long format using
[`reshape_chroms()`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md).

## Side effects

If `engine == "base"`, the plot is rendered to the active graphics
device.

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
