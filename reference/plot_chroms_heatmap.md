# Plot chromatograms as heatmap

Plots the specified traces from a list of chromatograms as a heatmap.

## Usage

``` r
plot_chroms_heatmap(
  chrom_list,
  idx = NULL,
  lambdas,
  engine = c("base", "ggplot", "plotly"),
  show_legend = TRUE,
  xlim = NULL,
  legend_position = "topright",
  title = "",
  show_ylabs = FALSE
)
```

## Arguments

- chrom_list:

  List of chromatograms to plot

- idx:

  A vector representing the names or numerical indices of the
  chromatograms to plot.

- lambdas:

  A character or numeric vector specifying the wavelengths to plot.

- engine:

  Plotting engine. Either `base`
  ([`matplot`](https://rdrr.io/r/graphics/matplot.html)),
  [`plotly`](https://rdrr.io/pkg/plotly/man/plotly.html), or
  [ggplot](https://ggplot2.tidyverse.org/reference/ggplot2-package.html).

- show_legend:

  Logical. Whether to display legend or not. Defaults to `TRUE`.

- xlim:

  Range of x axis.

- legend_position:

  Position of legend.

- title:

  Title for plot.

- show_ylabs:

  Logical. Whether to show y labels. Defaults to `FALSE`.

## Value

No return value, called for side effects.

## Details

Adapted from
[`plot.VPdtw`](https://rdrr.io/pkg/VPdtw/man/plot.VPdtw.html).

## Side effects

Plots the traces of the specified chromatograms `idx` at the specified
wavelengths `lambdas` as a heatmap. Plots can be produced using base
graphics engine, `ggplot2`, or `plotly`, according to the value of
`engine`.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_warp)
plot_chroms_heatmap(Sa_warp, lambdas = 210)
```
