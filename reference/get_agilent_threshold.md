# Calculate purity thresholds

Computes the threshold used by
[`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md):
the similarity expected between a spectrum and the reference if they
differed only by noise, the product of their reliabilities \\(1 -
\sigma^2\_{noise} / \sigma^2_i)(1 - \sigma^2\_{noise} /
\sigma^2\_{ref})\\, capped at `max_threshold`.

## Usage

``` r
get_agilent_threshold(
  x,
  pos,
  noise_inflation = 1,
  noise_variance = NULL,
  noise_method = c("quantile", "threshold"),
  noise_quantile = 0.1,
  noise_window = NULL,
  noise_threshold = NULL,
  max_threshold = 0.99,
  background_correction = c("linear", "constant", "none"),
  reference = c("average", "apex"),
  cutoff = 0.05,
  lambdas,
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

  Multiplier applied to the noise variance. Defaults to `1`.

- noise_variance:

  Variance of noise.

- noise_method:

  How to select the timepoints used to estimate the noise variance.
  Either `"quantile"` or `"threshold"`.

- noise_quantile:

  Quantile of maximum absorbance below which timepoints are treated as
  noise.

- noise_window:

  Optional range of row indices to use as the noise region.

- noise_threshold:

  Threshold to define noise. Used only when
  `noise_method == "threshold"`, where `NULL` (the default) means
  `0.005`.

- max_threshold:

  Ceiling on the threshold. Defaults to `0.99`.

- background_correction:

  How to correct the peak spectra for background absorbance. Either
  `"linear"`, `"constant"`, or `"none"`.

- reference:

  Which spectrum to compare against. Either `"average"` (the default) or
  `"apex"`. See
  [`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md).

- cutoff:

  Proportion of the corrected apex height below which timepoints are
  excluded from the average reference spectrum. Defaults to `0.05`.

- lambdas:

  Wavelengths to include. Character values are matched to the column
  names of `x`; numeric values are taken as column indices.

- weight:

  The `weight` argument is deprecated. Please use `noise_inflation`
  instead.

## Value

Returns a vector of purity thresholds, one for each timepoint from the
start to the end of the peak specified by `pos`, or `NULL` if the peak
cannot be assessed.

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
get_agilent_threshold(x, pos)
#>  [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.6221288 0.8321277 0.9184135
#>  [8] 0.9579810 0.9755804 0.9831036 0.9863842 0.9878527 0.9885007 0.9887190
#> [15] 0.9886271 0.9881708 0.9870904 0.9847495 0.9795133 0.9667497 0.9359573
#> [22] 0.8829043 0.8173311 0.7015591 0.4184964 0.0000000 0.0000000 0.0000000
#> [29] 0.0000000 0.0000000
```
