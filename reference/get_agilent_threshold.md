# Calculate purity thresholds

Intermediate function used in
[`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md)
to compute variance-based thresholds.

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

  A chromatographic matrix (timepoints × wavelengths).

- pos:

  A length-3 integer vector `c(apex, start, end)` giving the row indices
  of the peak apex, left boundary, and right boundary in `x`.

- weight:

  Scaling parameter affecting stringency of threshold. Defaults to `1`.

- noise_variance:

  Variance of noise.

- noise_threshold:

  Threshold to define noise. Highest proportion of maximum absorbance.
  Defaults to `0.005`.

- lambdas:

  Wavelengths to include.

## Value

Returns a vector of purity thresholds at each retention time index
within the peak specified by `pos`.

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
#>  [1] 0.08919129 0.04889951 0.02983187 0.00000000 0.00000000 0.00000000
#>  [7] 0.18706175 0.64143327 0.85443767 0.93649520 0.96786675 0.98059458
#> [13] 0.98603512 0.98835261 0.98907636 0.98863759 0.98678229 0.98243533
#> [19] 0.97297769 0.95181639 0.90252245 0.80466241 0.68708839 0.56780858
#> [25] 0.38645915 0.12893417 0.00000000 0.00000000 0.00000000 0.00000000
```
