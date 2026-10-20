# Plot matched peak pairs

Shows how peaks in each sample were matched to the reference by
[`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md).

## Usage

``` r
# S3 method for class 'matched_peak_list'
plot(x, idx = NULL, type = c("ladder", "scatter"), ...)

# S3 method for class 'matched_peak_list_set'
plot(x, ...)
```

## Arguments

- x:

  A `matched_peak_list` object from
  [`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md),
  or a `matched_peak_list_set`, whose wavelengths are plotted in turn.

- idx:

  Integer or character vector of sample indices or names to plot.
  Defaults to all non-reference samples.

- type:

  Plot type, either `"ladder"` (default) or `"scatter"`. See *Details*.

- ...:

  Additional arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) for
  `type = "ladder"`, such as `xlim` to zoom in on part of the run.

## Value

Invisibly returns `x`.

## Details

With `type = "ladder"` (the default), each sample gets one row, with
retention time on the x-axis. Dotted vertical lines mark the reference
peaks, and a short tick marks each one in every row. A matched sample
peak is drawn as a point filled by its spectral correlation with the
reference peak, joined to the reference tick by a horizontal bar, so the
length and direction of the bar show the retention time shift. A tick
without a bar is a reference peak left unmatched in that sample, and an
open circle is a sample peak left unmatched. Where a bar would overlap
another bar or run through an unmatched peak, it is lifted slightly
above its row.

With `type = "scatter"`, each sample gets its own panel, plotting every
candidate pair (a pair within the retention time window and above the
spectral cutoff) at its reference and sample retention times, so
retention time drift shows as distance from the \\y = x\\ diagonal.
Candidate pairs are colored by spectral correlation, matched pairs are
outlined in black, and unmatched peaks are marked on the axes.

## See also

[`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md),
[`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.chrom_list()`](https://ethanbass.github.io/chromatographR/reference/plot.chrom_list.md),
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
