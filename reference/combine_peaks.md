# Combine duplicate peaks

Identifies groups of duplicate peaks in a peak table and retains a
single representative peak from each group. Peaks are considered
duplicates when their retention times differ by less than `tol` and
their spectral correlation exceeds `min.cor`. This is useful for
collapsing peaks that were integrated at more than one wavelength or
chromatographic component.

## Usage

``` r
combine_peaks(
  peak_table,
  tol = 0.01,
  min.cor = 0.9,
  choose = c("max", "least_sparse", "lambda"),
  lambda = NULL,
  verbose = getOption("verbose")
)
```

## Arguments

- peak_table:

  Peak table from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- tol:

  Tolerance for matching retention times (maximum retention time
  difference). Defaults to `0.01`.

- min.cor:

  Minimum spectral correlation to confirm a match. Defaults to `0.9`.

- choose:

  Method used to select a representative peak from each group of
  duplicate peaks. `"max"` retains the peak with the greatest total
  intensity, `"least_sparse"` retains the peak detected in the greatest
  number of samples, and `"lambda"` retains peaks matching the preferred
  wavelength(s).

- lambda:

  Character vector of preferred wavelength(s) to retain when
  `choose = "lambda"`. The function keeps the duplicate peak whose
  integration wavelength matches one of the supplied values. An error is
  thrown if none of the supplied values match any peak in the peak
  table.

- verbose:

  Logical. Whether to print status to the console.

## Value

A peak table derived from the original, but with columns corresponding
to duplicate peaks combined according to the specified criteria.

## See also

[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)

Other utility functions:
[`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md),
[`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md),
[`merge_peaks()`](https://ethanbass.github.io/chromatographR/reference/merge_peaks.md),
[`reshape_chroms()`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md),
[`reshape_peaktable()`](https://ethanbass.github.io/chromatographR/reference/reshape_peaktable.md)

## Author

Ethan Bass

## Examples

``` r
data(pk_tab)
data(Sa_warp)
pk_tab <- attach_ref_spectra(pk_tab)
combine_peaks(pk_tab, tol = 0.02, min.cor = 0.9)
#>           V1       V2        V3        V4       V5        V6       V7        V8
#> 119 0.000000 5.111190 0.0000000 0.5864823 1.858693 15.079233 0.000000 0.4135257
#> 121 0.000000 4.033081 0.9163944 0.5221266 0.000000  7.705485 0.000000 0.0000000
#> 122 3.368245 8.228498 0.9735932 0.0000000 1.963956 25.888893 0.000000 1.8394472
#> 458 0.000000 5.655876 1.8309616 0.0000000 1.735467 14.650605 1.459654 0.3488550
#>            V9      V10      V11      V12       V13      V14      V15       V16
#> 119  79.71546 2.637016 41.60604 1.724824 1.3743838 2.273449 36.58533 15.189655
#> 121  57.27332 1.042586 18.33105 0.000000 0.8614246 7.504094 20.56684  8.879872
#> 122  59.62042 1.417466 35.83855 0.000000 1.2422142 4.241067 31.60616 12.910084
#> 458 108.93550 4.157266 33.22952 3.405194 3.1656321 5.231752 42.99198 23.968663
#>          V17      V18      V19       V20       V21      V22       V23
#> 119 33.71629 4.186043 3.890605 1.7999191 12.435065 15.09572 2.5155590
#> 121 22.59679 1.169449 1.750044 0.1143296  6.619292 12.54686 0.7921663
#> 122 28.40685 2.232165 3.037413 0.4855545  7.976247 14.88744 2.7807659
#> 458 25.33267 5.011990 5.348641 0.2741214 12.648599 15.42859 2.3459283
#>            V24        V25      V26       V27       V28       V29       V30
#> 119 0.00000000 0.11261310 17.80024 2.0024282 0.0000000 0.7981487 0.8774407
#> 121 0.01374328 0.43874747 14.41614 0.6184847 0.5505384 0.4542899 0.5315334
#> 122 0.01900430 0.11078837 16.60867 0.9699587 1.2199894 1.2251103 0.8634827
#> 458 0.00512731 0.02601551 24.77821 1.0894581 1.0049019 1.5687582 1.7486409
#>          V31      V32      V33       V34       V35       V36       V39
#> 119 1.819662 1.023753 1.541590 0.0000000 0.2105941 0.5041651 0.1046104
#> 121 1.085835 0.000000 1.811022 0.0000000 0.5440478 0.2934669 0.0000000
#> 122 1.288960 0.000000 4.117218 0.1280478 0.5797748 0.3649663 0.3733894
#> 458 1.459723 1.647317 3.730485 0.2736762 0.7127212 0.6514810 0.0000000
#>            V40       V42       V43        V44       V46      V49       V50
#> 119 0.09932078 0.2558299 14.841595 0.00000000 0.8899389 16.64689 0.0000000
#> 121 0.18110972 0.0000000  6.941289 0.03897730 1.8268337 11.45662 0.1379431
#> 122 0.00000000 0.3773357 14.120471 0.03016515 1.7424893 14.23811 0.7083038
#> 458 0.07382435 0.7167379 11.020983 0.53962131 0.7439920 12.98416 0.6251158
#>           V52      V54      V55         V57         V58       V60       V62
#> 119 1.1409371 5.373707 6.566692 0.000000000 0.059937260 0.1983225 0.3076778
#> 121 0.4724857 2.819725 5.966240 0.002892882 0.023415605 0.1066713 0.1992947
#> 122 1.1684493 3.438104 7.081215 0.004187678 0.041697073 0.2845064 0.0000000
#> 458 1.6551325 5.200684 7.292224 0.023957666 0.005351955 0.2704478 0.3774084
#>           V63       V64
#> 119 0.6844486 0.5443093
#> 121 0.9549370 0.0000000
#> 122 0.9325107 0.0000000
#> 458 0.5895916 0.6156869
```
