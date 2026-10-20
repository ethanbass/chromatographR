# Combine duplicate peaks

Identifies groups of duplicate peaks in a peak table and retains a
single representative peak from each group. Peaks are considered
duplicates when their retention times differ by less than `tol` and
their spectral correlation exceeds `min_cor`. Use it to collapse a
compound that was integrated at more than one wavelength into a single
feature.

## Usage

``` r
combine_peaks(
  peak_table,
  tol = 0.01,
  min_cor = 0.9,
  choose = c("max", "least_sparse", "lambda"),
  lambda = NULL,
  verbose = getOption("verbose"),
  min.cor = NULL
)
```

## Arguments

- peak_table:

  Peak table from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- tol:

  Tolerance for matching retention times (maximum retention time
  difference), in minutes. Defaults to `0.01`.

- min_cor:

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
  integration wavelength matches one of the supplied values; if several
  match, it keeps the first in the peak table, whatever the order of
  `lambda`. A group with no match falls back to `choose = "max"`. An
  error is thrown if none of the supplied values match any peak in the
  peak table.

- verbose:

  Logical. Whether to print status to the console.

- min.cor:

  The `min.cor` argument is deprecated. Please use `min_cor` instead.

## Value

A peak table derived from the original, but with columns corresponding
to duplicate peaks combined according to the specified criteria.

## See also

[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)

Other utility functions:
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
combine_peaks(pk_tab, tol = 0.02, min_cor = 0.9)
#>           V1       V2        V3        V4        V5        V6        V7
#> 119 0.000000 5.104960 0.0000000 0.6009758 15.228827 0.4598398  80.38146
#> 121 0.000000 3.989016 0.8661555 0.5518728  7.734209 0.0000000  57.49868
#> 122 3.360145 8.139282 0.9361464 2.0773758 25.864168 2.2679018  59.76982
#> 458 0.000000 5.696499 1.7732372 1.8204290 14.969887 0.2980591 110.89038
#>           V8       V9      V10      V11      V12      V13       V14      V15
#> 119 2.909701 42.02479 1.425470 1.561397 2.446028 37.29032 14.016502 33.91009
#> 121 1.080224 18.42549 0.000000 1.021819 8.088705 20.28485  8.652057 22.63074
#> 122 1.164881 36.05080 0.000000 1.552761 4.585314 32.38382 12.627866 28.68807
#> 458 4.284673 33.86360 4.100628 2.969487 5.969776 43.45102 24.337351 25.50186
#>          V16      V17       V18       V19      V20      V21         V22
#> 119 4.463381 5.263536 1.9309993 12.848093 14.94647 2.406808 0.000000000
#> 121 1.219043 1.823485 0.1217485  6.788176 12.21597 0.740373 0.007342364
#> 122 2.452972 3.209959 0.5920201  8.346982 14.60876 2.760261 0.000000000
#> 458 5.262720 5.705804 0.4092139 12.684912 15.18446 2.299992 0.000000000
#>             V23       V24      V25       V26       V27       V28       V29
#> 119 0.088458294 0.0000000 17.56996 1.7764505 0.0000000 1.2744367 0.8476725
#> 121 0.027751603 0.4036467 14.27335 0.5863850 0.5314740 0.4319262 0.5374353
#> 122 0.080516724 0.0000000 16.53137 0.8943947 1.0992387 0.0000000 1.1170351
#> 458 0.003952009 0.0000000 24.12003 0.9896030 0.9328859 1.4820802 1.5240373
#>          V30       V31      V32        V33       V34       V35       V38
#> 119 1.568035 0.9404432 1.450265 0.00000000 0.1798730 0.4618064 0.1179606
#> 121 1.087723 0.0000000 1.850937 0.00000000 0.5016980 0.2508542 0.0000000
#> 122 1.206059 0.0000000 3.916146 0.09581451 0.5626279 0.3415882 0.3624550
#> 458 0.000000 3.3193129 3.663237 0.17479033 0.6311858 0.5666003 0.0000000
#>            V39       V49       V50         V55       V58       V60       V61
#> 119 0.10407311 2.3029954 1.2029546 0.000000000 0.1794781 0.3020769 0.6494688
#> 121 0.20955663 0.8001079 0.5176482 0.003438554 0.1411927 0.1795083 0.9209173
#> 122 0.00000000 1.2558211 1.2144860 0.004641440 0.2659416 0.0000000 0.9074791
#> 458 0.08872707 3.0731105 1.8331265 0.021059462 0.2412048 0.2994864 0.9530536
```
