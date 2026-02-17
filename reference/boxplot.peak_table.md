# Make boxplot from peak table.

The function can take multiple response variables on the left hand side
of the formula (separated by `+`). In this case, a separate boxplot will
be produced for each response variable.

## Usage

``` r
# S3 method for class 'peak_table'
boxplot(x, formula, ...)
```

## Arguments

- x:

  A `peak_table` object.

- formula:

  A [`formula`](https://rdrr.io/r/stats/formula.html) object.

- ...:

  Additional arguments to
  [`boxplot`](https://rdrr.io/r/graphics/boxplot.html).

## See also

Other visualization functions:
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md)

## Examples

``` r
data(pk_tab)
path <- system.file("extdata", "Sa_metadata.csv", package = "chromatographR")
meta <- read.csv(path)
pk_tab <- attach_metadata(peak_table = pk_tab, metadata = meta, column="vial")
boxplot(pk_tab, formula=V11 ~ trt)
```
