# Trim peak

Trim peak

## Usage

``` r
trim_peak(x, pos, cutoff = 0.05)
```

## Arguments

- x:

  A numeric vector of signal intensities (a single wavelength trace).

- pos:

  A length-3 integer vector `c(apex, start, end)` giving the row indices
  of the peak apex, left boundary, and right boundary respectively.

- cutoff:

  Proportion of maximum absorbance to use as cutoff. Defaults to `0.05`.

## Value

Returns indices within the peak region (`start:end`) where the signal
exceeds `cutoff` times the apex intensity.

## Author

Ethan Bass

## Examples

``` r
data(Sa_warp)
x <- Sa_warp[[1]]
signal <- x[, "210"]
pos <- c(37, 22, 51)
trim_peak(signal, pos)
#> 10.52 10.54 10.56 10.58  10.6 10.62 10.64 10.66 10.68  10.7 10.72 10.74 10.76 
#>     6     7     8     9    10    11    12    13    14    15    16    17    18 
#> 10.78  10.8 10.82 10.84 10.86 10.88  10.9 10.92 10.94 10.96 10.98    11 
#>    19    20    21    22    23    24    25    26    27    28    29    30 
```
