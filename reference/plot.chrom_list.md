# Plot a list of chromatograms

Plots a `chrom_list`, as returned by
[`chromConverter::read_chroms()`](https://ethanbass.github.io/chromConverter/reference/read_chroms.html),
as traces
([`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md)),
a heatmap of samples
([`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md)),
or a contour plot of one sample
([`plot_chrom_contour()`](https://ethanbass.github.io/chromatographR/reference/plot_chrom_contour.md)).
If each sample holds several elements, the `"dad"` element is plotted,
or for `type = "traces"` the one named by `what`.

## Usage

``` r
# S3 method for class 'chrom_list'
plot(x, type = c("traces", "heatmap", "contour"), ...)
```

## Arguments

- x:

  A `chrom_list` object.

- type:

  Type of plot: `"traces"` (the default), `"heatmap"`, or `"contour"`.

- ...:

  Arguments to the plotting function selected by `type`. `lambdas` names
  the traces to plot for `"traces"`, and a single wavelength for
  `"heatmap"`; `idx` selects samples for `"traces"` and `"heatmap"`, and
  one sample for `"contour"`.

## Value

The value returned by the selected plotting function.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
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
data(Sa_pr)
class(Sa_pr) <- c("chrom_list", "list")
plot(Sa_pr, lambdas = 210)

plot(Sa_pr, type = "heatmap", lambdas = 210)

plot(Sa_pr, type = "contour", idx = 1)
```
