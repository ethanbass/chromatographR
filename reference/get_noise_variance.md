# Calculate variance of noise regions

Estimates the across-wavelength variance of spectral noise regions used
in peak purity calculations.

## Usage

``` r
get_noise_variance(
  x,
  noise_method = c("quantile", "threshold"),
  noise_quantile = 0.1,
  noise_window = NULL,
  noise_threshold = 0.005,
  lambdas
)
```

## Arguments

- x:

  A chromatogram in matrix format

- noise_method:

  Either `"quantile"` or `"threshold"`.

- noise_quantile:

  Quantile of maximum absorbance below which timepoints are treated as
  noise.

- noise_window:

  Optional range of row indices to use as the noise region. Takes
  precedence over `noise_method`.

- noise_threshold:

  Threshold to define noise. Highest proportion of maximum absorbance.
  Used only when `noise_method == "threshold"`.

- lambdas:

  Wavelengths to include.

## Value

Returns the average across-wavelength variance of the signal over the
retention times treated as noise, or `NA_real_` if no timepoints
qualify.

## Note

The variance is measured across wavelengths, not across time, since that
is the quantity the purity threshold compares against the variance of
each spectrum. It is therefore not directly comparable to the temporal
noise estimate used by
[`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md)
and
[`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md),
and the two diverge on data that has been smoothed along the wavelength
axis by
[`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md).

## Author

Ethan Bass
