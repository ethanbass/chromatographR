# Calculate peak purity

Estimates peak purity by assessing spectral dissimilarity across the
peak region, following Stahl (2003).

## Usage

``` r
get_purity(
  x,
  pos,
  noise_inflation = 1,
  cutoff = 0.05,
  noise_variance = NULL,
  noise_method = c("quantile", "threshold"),
  noise_quantile = 0.1,
  noise_window = NULL,
  noise_threshold = NULL,
  max_threshold = 0.99,
  background_correction = c("linear", "constant", "none"),
  reference = c("average", "apex"),
  metric = c("fraction", "max"),
  min_run = 3,
  lambdas,
  try = TRUE,
  weight = NULL
)
```

## Arguments

- x:

  A chromatographic matrix (timepoints × wavelengths).

- pos:

  A length-3 integer vector `c(apex, start, end)` giving the row indices
  of the peak apex, left boundary, and right boundary in `x`.

- noise_inflation:

  Multiplier applied to the noise variance, whether supplied via
  `noise_variance` or estimated. Values above `1` treat the noise as
  larger than measured, which loosens the threshold. Defaults to `1`.

- cutoff:

  Timepoints at or below `cutoff` times the background-corrected apex
  height are excluded, measured at the wavelength of maximum absorbance
  at the apex. Defaults to `0.05`.

- noise_variance:

  Variance of noise. If provided, it replaces the estimate, and the
  `noise_*` arguments that control estimation are ignored. Estimating
  the noise variance requires a pass over the whole chromatogram, so
  supplying it is much faster when many peaks are assessed.

- noise_method:

  How to select the timepoints used to estimate the noise variance.
  Either `"quantile"` (the default), which takes the timepoints in the
  lowest `noise_quantile` of maximum absorbance, or `"threshold"`, the
  legacy behavior, which takes timepoints falling below
  `noise_threshold` times the largest absorbance in the chromatogram.

- noise_quantile:

  Quantile of maximum absorbance below which timepoints are treated as
  noise. Defaults to `0.1`.

- noise_window:

  Optional range of row indices to use as the noise region. Takes
  precedence over `noise_method`.

- noise_threshold:

  Threshold to define noise, as a proportion of the largest absorbance
  in the chromatogram. Used only when `noise_method == "threshold"`,
  where `NULL` (the default) means `0.005`.

- max_threshold:

  Ceiling on the purity threshold. Without it, the threshold approaches
  one for high signal-to-noise peaks, so the tolerated dissimilarity
  approaches zero, and spectral variation that noise cannot explain,
  such as detector non-linearity or a changing solvent composition,
  makes such peaks impure. Corresponds to the user-defined similarity
  threshold in the Agilent ChemStation, where `0.99` is a similarity
  of 990. Defaults to `0.99`. Set to `1` to disable.

- background_correction:

  How to correct the peak spectra for background absorbance. Either
  `"linear"` (the default), `"constant"`, or `"none"`.

- reference:

  Which spectrum to compare the peak spectra against. Either `"average"`
  (the default), the mean of the assessed spectra, or `"apex"`, the
  spectrum at the peak apex. The apex spectrum correlates perfectly with
  itself, so it forces a purity ratio of exactly zero at the apex of
  every peak. Stahl recommends the average spectrum for unknown
  impurities.

- metric:

  Statistic to summarize the purity ratios. Either `"fraction"` (the
  default), the proportion of in-peak timepoints that are pure, or
  `"max"`, the largest purity ratio. The run-length rule applies to
  `"fraction"` only.

- min_run:

  Minimum number of consecutive failing timepoints required to count as
  impure. Defaults to `3`.

- lambdas:

  Wavelengths to include in calculations. Character values are matched
  to the column names of `x`; numeric values are taken as column
  indices. Defaults to all wavelengths.

- try:

  Logical. If `FALSE`, returns `NA` without estimating purity. Defaults
  to `TRUE`.

- weight:

  The `weight` argument is deprecated. Please use `noise_inflation`
  instead.

## Value

Returns the purity of the peak specified by `pos`: with
`metric = "fraction"`, the proportion of timepoints above `cutoff` not
counted as impure; with `metric = "max"`, the largest purity ratio.
Returns `NA_real_` if the peak cannot be assessed, including when fewer
than three timepoints lie above `cutoff`. An error during the
calculation is converted to a warning and also returns `NA_real_`.

## Details

At each timepoint in the peak, the spectrum is compared to a reference
spectrum, by default the mean of the assessed spectra. The resulting
dissimilarity (\\1 - r^2\\) is divided by the dissimilarity that noise
alone could explain, one minus the threshold returned by
[`get_agilent_threshold`](https://ethanbass.github.io/chromatographR/reference/get_agilent_threshold.md).
A ratio below one means the spectrum is as similar to the reference as
the noise permits. With the default `metric = "fraction"`, a timepoint
only counts as impure if it belongs to a run of at least `min_run`
consecutive failing timepoints, so isolated failures are ignored. Stahl
requires a run of three.

Spectra are background-corrected before comparison by subtracting the
spectra at the peak boundaries, interpolated across the peak
(`background_correction = "linear"`) or averaged
(`background_correction = "constant"`). Anything that changes linearly
across the peak is removed with the background, including a co-eluting
impurity much broader than the peak. If
[`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md)
has removed the baseline completely, `background_correction = "none"`
avoids this, but any remaining baseline then enters the comparison as
part of each spectrum.

With the default `metric = "fraction"`, pure peaks do not always score
exactly one, so treat purity below about 0.9, rather than any value
below one, as a sign of impurity. Wide, tailing and weak peaks fall
below 0.9 more often than others, even when pure. So do peaks above
about 2 AU, whose spectra near the apex are distorted by the detector's
non-linear response. The purity calculation is still being refined, and
its values may change in future releases.

## References

Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array
Technology.” Agilent Technologies, April 1, 2003, 16.
[5988-8647EN.pdf](https://web.archive.org/web/20220615145615/https://www.agilent.com/cs/library/applications/5988-8647EN.pdf)

## Author

Ethan Bass

## Examples

``` r
data(Sa_warp)
x <- Sa_warp[[1]]
pos <- c(37, 22, 51)
get_purity(x, pos)
#> [1] 0.4545455
```
