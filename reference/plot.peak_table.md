# Plot spectrum from peak table

Plots the trace and/or spectrum for a given peak in peak table.

## Usage

``` r
# S3 method for class 'peak_table'
plot(
  x,
  loc,
  chrom_list,
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
  chr = NULL,
  ...
)
```

## Arguments

- x:

  The peak table (output from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)
  function).

- loc:

  A vector specifying the peak(s) or retention time(s) that you wish to
  plot.

- chrom_list:

  A list of chromatograms in matrix format (timepoints x wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object used to create the `peak_table`.

- what:

  What to look for. Either `peak` to extract spectral information for a
  certain peak, `rt` to scan by retention time, or `click` to manually
  select retention time by clicking on the chromatogram. Defaults to
  `peak`.

- idx:

  Numerical index of chromatogram you wish to plot; "max" to plot the
  chromatogram with the largest signal; or "all" to plot spectra for all
  chromatograms.

- lambda:

  The wavelength you wish to plot the trace at (if `plot_chrom` is TRUE
  and/or the wavelength to be used for the determination of signal
  abundance.

- plot_spectrum:

  Logical. If TRUE, plots the spectrum of the chosen peak. Defaults to
  TRUE.

- plot_trace:

  Logical. If TRUE, plots the trace of the chosen peak at lambda.
  Defaults to TRUE.

- box_plot:

  Logical. If TRUE, plots box plot using categories defined by `vars`.

- vars:

  Independent variables for boxplot. Righthand side of formula.

- spectrum_labels:

  Logical. If TRUE, plots labels on maxima in spectral plot. Defaults to
  TRUE.

- scale_spectrum:

  Logical. If TRUE, scales spectrum to unit height. Defaults to FALSE.

- export_spectrum:

  Logical. If TRUE, exports spectrum to console. Defaults to FALSE.

- verbose:

  Logical. If TRUE, prints verbose output to console. Defaults to TRUE.

- engine:

  Which plotting engine to use: either `base` or `plotly`.

- chr:

  Deprecated. Please use `idx` instead.

- ...:

  Additional arguments to
  [`boxplot`](https://rdrr.io/r/graphics/boxplot.html).

## Value

If `export_spectrum` is TRUE, returns the spectrum as a ` data.frame`
with wavelengths as rows and columns encoding the absorbance (or
normalized absorbance, if `scale_spectrum` is TRUE) for the specified
sample(s). Otherwise, there is no return value.

## Details

Can be used to confirm the identity of a peak or check that a particular
column in the peak table represents a single compound. Can also be used
to create simple box-plots to examine the distribution of a peak with
respect to variables defined in sample metadata.

## Side effects

If `plot_trace` is `TRUE`, plots the chromatographic trace of the
specified chromatogram (`idx`), at the specified wavelength (`lambda`)
with a dotted red line to indicate the retention time given by `loc`.
The trace is a single column from the chromatographic matrix.

If `plot_spectrum` is TRUE, plots the spectrum for the specified
chromatogram at the specified retention time. The spectrum is a single
row from the chromatographic matrix.

If `box_plot` is TRUE, produces a
[`boxplot`](https://rdrr.io/r/graphics/boxplot.html) from the specified
peak with groups provided by `vars`.

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
