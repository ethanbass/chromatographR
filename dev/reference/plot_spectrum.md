# Plot spectrum from peak table

Plots the trace and/or spectrum for a given peak or retention time in a
`peak_table` object or a list of chromatograms.

## Usage

``` r
plot_spectrum(
  loc = NULL,
  peak_table,
  chrom_list,
  idx = "max",
  lambda = "max",
  plot_spectrum = TRUE,
  plot_trace = TRUE,
  spectrum_labels = TRUE,
  scale_spectrum = FALSE,
  export_spectrum = FALSE,
  verbose = TRUE,
  what = c("peak", "rt", "idx", "click"),
  engine = c("base", "plotly", "ggplot2"),
  chr = NULL,
  ...
)
```

## Arguments

- loc:

  The name of the peak or retention time for which you wish to extract
  spectral data.

- peak_table:

  The peak table (output from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/dev/reference/get_peaktable.md)).

- chrom_list:

  A list of chromatograms in matrix format (timepoints x wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object used to create the provided `peak_table`.

- idx:

  Numerical index of chromatogram you wish to plot, or "max" to
  automatically plot the chromatogram with the largest signal at the
  given peak or retention time.

- lambda:

  The wavelength you wish to plot the trace at if plot_trace == TRUE
  and/or the wavelength to be used for the determination of signal
  abundance.

- plot_spectrum:

  Logical. If `TRUE`, plots the spectrum of the chosen peak. Defaults to
  `TRUE`.

- plot_trace:

  Logical. If `TRUE`, plots the trace of the chosen peak at lambda.
  Defaults to `TRUE`.

- spectrum_labels:

  Logical. If `TRUE`, plots labels on maxima in spectral plot. Defaults
  to `TRUE`.

- scale_spectrum:

  Logical. If `TRUE`, scales spectrum to unit height. Defaults to
  `FALSE`.

- export_spectrum:

  Logical. If `TRUE`, exports spectrum to console. Defaults to `FALSE`.

- verbose:

  Logical. If `TRUE`, prints verbose output to console. Defaults to
  `TRUE`.

- what:

  What to look for. Either `peak` to extract spectral information for a
  certain peak, `rt` to scan by retention time, `idx` to scan by numeric
  index, or `click` to manually select retention time by clicking on the
  chromatogram. Defaults to "peak" mode.

- engine:

  Which plotting engine to use: `base`, `ggplot2`, or `plotly`.

- chr:

  Deprecated. Please use `idx` instead.

- ...:

  Additional arguments.

## Value

If `export_spectrum` is TRUE, returns the spectrum as a ` data.frame`
with wavelengths as rows and a single column encoding the absorbance (or
normalized absorbance, if `scale_spectrum` is TRUE) at each wavelength.
If `export_spectrum` is FALSE, the output depends on the plotting
`engine`. If `engine == "plotly"`, returns a `plotly` object containing
the specified plots. Otherwise, if `engine == "base"`, there is no
return value.

## Details

Can be used to confirm the identity of a peak or check that a particular
column in the peak table represents a single compound. Retention times
can also be selected by clicking on the plotted trace if
`what == 'click'`. Plots can be produced using either base R graphics,
[`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html),
or `plotly`, according to the value of the `engine` argument.

## Side effects

- If `plot_trace` is `TRUE`, plots the chromatographic trace of the
  specified chromatogram (`idx`), at the specified wavelength (`lambda`)
  with a dotted red line to indicate the retention time given by `loc`.
  The trace is a single column from the chromatographic matrix.

- If `plot_spectrum` is `TRUE`, plots the spectrum for the specified
  chromatogram at the specified retention time. The spectrum is a single
  row from the chromatographic matrix.

## See also

Other visualization functions:
[`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/dev/reference/boxplot.peak_table.md),
[`mirror_plot()`](https://ethanbass.github.io/chromatographR/dev/reference/mirror_plot.md),
[`plot.peak_list()`](https://ethanbass.github.io/chromatographR/dev/reference/plot.peak_list.md),
[`plot.peak_table()`](https://ethanbass.github.io/chromatographR/dev/reference/plot.peak_table.md),
[`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_all_spectra.md),
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_chroms.md),
[`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/dev/reference/plot_chroms_heatmap.md),
[`scan_chrom()`](https://ethanbass.github.io/chromatographR/dev/reference/scan_chrom.md)

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa)
pks <- get_peaks(Sa, lambda = "220.00000")
pk_tab <- get_peaktable(pks)
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(2, 1))
plot_spectrum(loc = "V10", peak_table = pk_tab, what = "peak")
par(oldpar)
}
```
