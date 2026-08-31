# Calculate peak purity

Estimates peak purity by assessing spectral dissimilarity across the
peak region, following Stahl (2003).

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

  A chromatographic matrix (timepoints × wavelengths).

- pos:

  A length-3 integer vector `c(apex, start, end)` giving the row indices
  of the peak apex, left boundary, and right boundary in `x`.

- weight:

  Weight parameter passed to
  [`get_agilent_threshold`](https://ethanbass.github.io/chromatographR/reference/get_agilent_threshold.md).

- cutoff:

  Proportion of maximum absorbance to use as cutoff. Passed to
  [`trim_peak`](https://ethanbass.github.io/chromatographR/reference/trim_peak.md).
  Defaults to `0.05`.

- noise_variance:

  Variance of noise. Passed to `get_agilent_threshold`.

- noise_threshold:

  Threshold to define noise. Highest proportion of maximum absorbance.
  Defaults to `0.01`.

- lambdas:

  Wavelengths to include in calculations.

- try:

  Logical. Whether to estimate the purity or not. Defaults to `TRUE`.

## Value

Returns the mean purity of the peak specified by `pos`, defined as the
proportion of timepoints with purity values below 1.

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
#> [1] 0.6666667
```
