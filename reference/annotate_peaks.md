# Annotate peaks on a chromatogram

Adds text labels to peaks on a chromatogram plot produced by
[`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md).
Peak retention times and intensities are looked up automatically from
the peak table.

## Usage

``` r
annotate_peaks(
  p,
  loc,
  peak_table,
  chrom_list = NULL,
  label = NULL,
  lambda = NULL,
  idx = NULL,
  vjust = -0.5,
  ...
)
```

## Arguments

- p:

  A `ggplot` object produced by
  [`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md).

- loc:

  A character vector of peak names (e.g. `c("V8", "V15")`). If named,
  the names are used as labels (e.g.
  `c(Sinigrin = "V8", `1ME` = "V15")`).

- peak_table:

  A `peak_table` object created by
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- chrom_list:

  A list of chromatograms in matrix format (timepoints x wavelengths).
  If no argument is provided here, the function will try to find the
  `chrom_list` object using the pointer in the `peak_table`.

- label:

  Labels to display at each peak. Can be a character vector of labels
  (one per peak), `"rt"` to use retention times from the peak table. If
  labels are not specified here, the names of `loc` will be applied, or
  if the `loc` vector is not named, the peak names will be used
  directly.

- lambda:

  Wavelength(s) to use for locating the peak apex. Inherited from `p` if
  `NULL` (default). If multiple wavelengths are provided, the one with
  the highest absorbance is used.

- idx:

  Index of chromatogram(s) to use for locating the peak apex. Inherited
  from `p` if `NULL` (default).

- vjust:

  Vertical justification of the label relative to the peak apex.
  Defaults to `-0.5`.

- ...:

  Additional arguments passed to
  [`ggplot2::annotate()`](https://ggplot2.tidyverse.org/reference/annotate.html).

## Value

A `ggplot` object.

## See also

[`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md),
[`plot_spectrum_inset()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum_inset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot_chroms(dat, idx = 1:5, lambda = 228, engine = "ggplot") |>
  annotate_peaks(c(Sinigrin = "V8", `1ME` = "V15"), peak_table = pktab)
} # }
```
