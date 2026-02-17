# Calculate purity thresholds

Calculate purity thresholds

## Usage

``` r
get_agilent_threshold(
  x,
  pos,
  weight = 1,
  noise_variance = NULL,
  noise_threshold = 0.005,
  lambdas
)
```

## Arguments

- x:

  A chromatogram in matrix format.

- pos:

  A vector containing peak information.

- weight:

  Scaling parameter affecting stringency of threshold. Defaults to `1`.

- noise_variance:

  Variance of noise.

- noise_threshold:

  Threshold to define noise. Highest proportion of maximum absorbance.
  Defaults to `.005`.

- lambdas:

  Wavelengths to include.

## Value

Returns a vector of purity thresholds at each retention time index
within the peak specified by `pos`.

## References

Stahl, Mark. “Peak Purity Analysis in HPLC and CE Using Diode-Array
Technology.” Agilent Technologies, April 1, 2003, 16.
<https://www.agilent.com/cs/library/applications/5988-8647EN.pdf>

## Author

Ethan Bass
