# Plot spectrum from peak table

Plots the trace and/or spectrum for a given peak in peak table.

## Usage

``` r
# S3 method for class 'peak_table'
plot(
  x,
  loc,
  chrom_list = NULL,
  what = "peak",
  idx = "max",
  lambda = "max",
  plot_spectrum = TRUE,
  plot_trace = TRUE,
  box_plot = FALSE,
  vars = NULL,
  spectrum_labels = TRUE,
  scale_spectrum = FALSE,
  export_spectrum = FALSE,
  verbose = TRUE,
  engine = c("base", "plotly", "ggplot"),
  ...
)
```

## Arguments

- x:

  The peak table (output from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)).

- loc:

  A vector specifying the peak(s) or retention time(s) to plot. If
  missing when `what = "peak"`, you are prompted for a peak name. When a
  value is returned (see below), only the first element of `loc` is
  used, with a warning.

- chrom_list:

  A list of chromatograms in matrix format (timepoints × wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object using the pointer in the `peak_table`.

- what:

  What to look for. Either `peak` to extract spectral information for a
  certain peak, `rt` to scan by retention time, `idx` to scan by numeric
  index, or `click` to manually select retention time by clicking on the
  chromatogram. Defaults to `"peak"` mode. `"click"` requires
  `engine = "base"`.

- idx:

  Numerical index of chromatogram you wish to plot; `"max"` (default) to
  plot the chromatogram with the largest value for the peak in the peak
  table; or `"all"` to plot spectra for all chromatograms.

- lambda:

  The wavelength at which to plot the trace. `"max"` (default) selects
  the wavelength at which the spectrum at `loc` is highest.

- plot_spectrum:

  Logical. If `TRUE`, plots the spectrum of the chosen peak. Defaults to
  `TRUE`.

- plot_trace:

  Logical. If `TRUE`, plots the trace of the chosen peak at lambda.
  Defaults to `TRUE`.

- box_plot:

  Logical. If `TRUE`, plots box plot using categories defined by `vars`.

- vars:

  Character vector of column names in the sample metadata, used as the
  right-hand side of the boxplot formula.

- spectrum_labels:

  Logical. If `TRUE`, plots labels on maxima in spectral plot. Defaults
  to `TRUE`.

- scale_spectrum:

  Logical. If `TRUE`, scales spectrum to unit height. Defaults to
  `FALSE`.

- export_spectrum:

  Logical. If `TRUE`, invisibly returns the spectrum as a `data.frame`.
  Defaults to `FALSE`.

- verbose:

  Logical. If `TRUE` (default), prints the index and name of each
  chromatogram plotted, with its retention time, wavelength and
  absorbance.

- engine:

  Plotting backend to use. One of `"base"`,
  [`"ggplot"`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html),
  or [`"plotly"`](https://rdrr.io/pkg/plotly/man/plotly.html).
  `"ggplot2"` is also accepted here, for backward compatibility, but
  `"ggplot"` is the spelling the other plotting functions take.

- ...:

  Additional arguments to
  [boxplot](https://rdrr.io/r/graphics/boxplot.html).

## Value

- If `export_spectrum = FALSE` (default), a `plotly` or `ggplot` object,
  or nothing if `engine == "base"`.

- If `export_spectrum = TRUE`, the spectrum as a `data.frame` with
  wavelengths as rows and one column per sample giving absorbance at
  each wavelength (normalized if `scale_spectrum = TRUE`).

## Details

Use it to confirm the identity of a peak, or to check that a particular
column in the peak table represents a single compound. It also draws
simple box plots of the distribution of a peak against variables defined
in sample metadata.

When `plot_trace` is `TRUE`, plots the trace of chromatogram `idx` at
wavelength `lambda`, with a dotted red line at the retention time given
by `loc`. When `plot_spectrum` is `TRUE`, plots the spectrum of the same
chromatogram at that retention time. When `box_plot` is `TRUE`, draws a
[`boxplot`](https://rdrr.io/r/graphics/boxplot.html) of the peak with
groups defined by `vars`.

## Side effects

If `engine == "base"`, plots are rendered to the active graphics device.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.chrom_list()`](https://ethanbass.github.io/chromatographR/reference/plot.chrom_list.md),
[`plot.matched_peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.matched_peak_list.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md),
[`plot_chrom_contour()`](https://ethanbass.github.io/chromatographR/reference/plot_chrom_contour.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md)

## Author

Ethan Bass
