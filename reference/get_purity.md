# Calculate mean peak purity

Estimates peak purity by assessing the dissimilarity of the spectra
comprising the peak, using the method described in Stahl 2003.

## Usage

``` r
get_purity(
  x,
  pos,
  weight = 1,
  cutoff = 0.05,
  noise_variance = NULL,
  noise_threshold = 0.01,
  lambdas,
  try = TRUE
)
```

## Arguments

- x:

  A chromatogram in matrix format.

- pos:

  A vector containing the center, lower and upper bounds of a peak as
  numeric indices.

- weight:

  Weight provided to
  [`get_agilent_threshold`](https://ethanbass.github.io/chromatographR/reference/get_agilent_threshold.md).

- cutoff:

  Proportion of maximum absorbance to use as cutoff. Argument to
  [`trim_peak`](https://ethanbass.github.io/chromatographR/reference/trim_peak.md).
  Defaults to `.05`.

- noise_variance:

  Variance of noise. Argument to
  [`get_agilent_threshold`](https://ethanbass.github.io/chromatographR/reference/get_agilent_threshold.md).

- noise_threshold:

  Threshold to define noise. Highest proportion of maximum absorbance.
  Defaults to `.01`.

- lambdas:

  Wavelengths to include in calculations.

- try:

  Logical. Whether to estimate the purity or not. Defaults to TRUE.

## Value

Returns the mean purity of the peak specified by `pos`, defined as the
proportion of timepoints with purity values below 1.

## References

Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array
Technology.” Agilent Technologies, April 1, 2003, 16.
<https://www.agilent.com/cs/library/applications/5988-8647EN.pdf>

## Author

Ethan Bass
