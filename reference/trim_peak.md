# Trim peak

Trim peak

## Usage

``` r
trim_peak(x, pos, cutoff = 0.05)
```

## Arguments

- x:

  A chromatogram in matrix format.

- pos:

  A vector containing peak information.

- cutoff:

  Proportion of maximum absorbance to use as cutoff. Defaults to `0.05`.

## Value

Returns indices within the peak specified by `pos` with a higher signal
intensity than the specified cutoff.

## Author

Ethan Bass
