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

  A vector specifying the peak(s) or retention time(s) to plot.

- chrom_list:

  A list of chromatograms in matrix format (timepoints x wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object using the pointer in the `peak_table`.

- what:

  What to look for. Either `peak` to extract spectral information for a
  certain peak, `rt` to scan by retention time, `idx` to scan by numeric
  index, or `click` to manually select retention time by clicking on the
  chromatogram. Defaults to `"peak"` mode.

- idx:

  Numerical index of chromatogram you wish to plot; "max" to plot the
  chromatogram with the largest signal; or "all" to plot spectra for all
  chromatograms.

- lambda:

  The wavelength you wish to plot the trace at (if `plot_chrom` is
  `TRUE`. Otherwise, the wavelength to be used for the determination of
  signal abundance.

- plot_spectrum:

  Logical. If `TRUE`, plots the spectrum of the chosen peak. Defaults to
  `TRUE`.

- plot_trace:

  Logical. If `TRUE`, plots the trace of the chosen peak at lambda.
  Defaults to `TRUE`.

- box_plot:

  Logical. If `TRUE`, plots box plot using categories defined by `vars`.

- vars:

  Independent variables for boxplot. Righthand side of formula.

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

  Logical. If `TRUE` (default), prints verbose output to the console.

- engine:

  Which plotting engine to use: `base`, `ggplot2`, or `plotly`.

- ...:

  Additional arguments to
  [boxplot](https://rdrr.io/r/graphics/boxplot.html).

## Value

- If `export_spectrum = FALSE` (default), a `plotly` or `ggplot` object,
  or nothing if `engine == "base"`.

- If `export_spectrum = TRUE`, invisibly returns the spectrum as a
  `data.frame` with wavelengths as rows and a single column per sample
  encoding absorbance at each wavelength (normalized if
  `scale_spectrum = TRUE`). Otherwise, if `engine == "base"`, there is
  no return value. If `export_spectrum` is `TRUE`, returns the spectrum
  as a `data.frame` with wavelengths as rows and columns encoding the
  absorbance (or normalized absorbance, if `scale_spectrum` is `TRUE`)
  for the specified sample(s).

## Details

Can be used to confirm the identity of a peak or check that a particular
column in the peak table represents a single compound. Can also be used
to create simple box-plots to examine the distribution of a peak with
respect to variables defined in sample metadata.

When `plot_trace` is `TRUE`, plots the chromatographic trace of the
specified chromatogram (`idx`), at the specified wavelength (`lambda`)
with a dotted red line to indicate the retention time given by `loc`.
The trace is a single column from the chromatographic matrix. When
`plot_spectrum` is `TRUE`, plots the spectrum for the specified
chromatogram at the specified retention time. The spectrum represents a
single row from the chromatographic matrix. When `box_plot` is `TRUE`,
produces a [`boxplot`](https://rdrr.io/r/graphics/boxplot.html) from the
specified peak with groups defined by the `vars` argument.

## Side effects

If `engine == "base"`, plots are rendered to the active graphics device.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md),
[`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md)

## Author

Ethan Bass
