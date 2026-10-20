# Make mirror plot from peak table

Plots chromatograms as a mirror plot.

## Usage

``` r
mirror_plot(
  x,
  chrom_list,
  lambdas = NULL,
  var,
  subset = NULL,
  print_legend = TRUE,
  legend_txt = NULL,
  legend_pos = "topright",
  legend_size = 1,
  mirror = TRUE,
  xlim = NULL,
  ylim = NULL,
  ...
)
```

## Arguments

- x:

  The peak table (output from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)).

- chrom_list:

  A list of chromatograms in matrix format (timepoints × wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object used to create the `peak_table`.

- lambdas:

  The wavelength(s) at which to plot the traces. Required for
  multi-wavelength chromatograms; ignored for one-dimensional
  chromatograms.

- var:

  Name of the column in the sample metadata that defines the two groups.

- subset:

  Character vector of the two levels of `var` to compare. Required if
  `var` has more than two levels.

- print_legend:

  Logical. Whether to print legend. Defaults to `TRUE`.

- legend_txt:

  Character vector containing labels for legend.

- legend_pos:

  Legend position, used only when `mirror = FALSE`. The mirror layout
  always places the legends at the top and bottom right.

- legend_size:

  Legend size (`cex` argument). Default is `1`.

- mirror:

  Logical. Whether to plot as mirror or stacked plots. Defaults to
  `TRUE`.

- xlim:

  Numerical vector specifying limits for x axis.

- ylim:

  Numerical vector specifying limits for y axis.

- ...:

  Additional arguments to
  [`matplot`](https://rdrr.io/r/graphics/matplot.html).

## Value

No return value, called for side effects.

## Details

Compares the chromatograms in two groups, defined by the levels of
`var`. When `mirror` is `TRUE` both are drawn in one panel, with the
second inverted below the axis; when it is `FALSE` they are drawn in two
stacked panels.

## Side effects

Plots are rendered to the active graphics device.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`plot.chrom_list()`](https://ethanbass.github.io/chromatographR/reference/plot.chrom_list.md),
[`plot.matched_peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.matched_peak_list.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md),
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
data(pk_tab)
path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
meta <- read.csv(path)
pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
mirror_plot(pk_tab,lambdas = c("210","260"), var = "trt", mirror = TRUE, 
  col = c("green","blue"))
```
