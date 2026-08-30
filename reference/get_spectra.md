# Extract a spectrum from a chromatogram

Extracts the UV-Vis spectrum at a specified retention time from one or
more chromatograms in a `chrom_list`.

## Usage

``` r
get_spectra(
  loc,
  peak_table,
  chrom_list = NULL,
  idx = seq_along(chrom_list),
  lambda = "max",
  scale_spectrum = FALSE,
  what = c("peak", "rt", "idx"),
  format = c("wide", "long")
)
```

## Arguments

- loc:

  Peak or retention time to extract the spectrum at. Interpretation
  depends on `what`.

- peak_table:

  A `peak_table` object. Required if `what == "peak"`.

- chrom_list:

  A list of chromatographic matrices.

- idx:

  Indices of chromatograms to extract spectra from. Defaults to all
  chromatograms in `chrom_list`.

- lambda:

  Wavelength at which to extract the spectrum. Use `"max"` to select the
  wavelength with the highest absorbance.

- scale_spectrum:

  Logical. Whether to normalize absorbance values to the range \[0, 1\].
  Defaults to `FALSE`.

- what:

  How to interpret `loc`. One of `"peak"` (look up retention time from
  `peak_table`), `"rt"` (treat `loc` as a retention time directly), or
  `"idx"` (treat `loc` as a row index).

- format:

  One of `"wide"` or `"long"`. In wide format, wavelengths are rows and
  each column is a sample. In long format, columns are `wavelength`,
  `absorbance`, `sample`, and `rt`. Defaults to `"wide"`.

## Value

A `data.frame` in the format specified by `format`.
