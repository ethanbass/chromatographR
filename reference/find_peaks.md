# Find peaks

Detects peaks in a chromatographic profile by identifying zero-crossings
in the smoothed first derivative of the signal (`y`) that exceed a
specified slope threshold (`slope_thresh`). Smoothing reduces spurious
local extrema that do not represent true features. Peaks can be filtered
using a minimum amplitude threshold (`amp_thresh`), which removes
low-intensity peaks.

## Usage

``` r
find_peaks(
  y,
  smooth_type = c("gaussian", "box", "savgol", "mva", "tmva", "none"),
  smooth_window = 0.001,
  slope_thresh = 0,
  amp_thresh = 0,
  bounds = TRUE
)
```

## Arguments

- y:

  A chromatographic signal (as a numeric vector).

- smooth_type:

  Type of smoothing. Either gaussian kernel (`"gaussian"`), box kernel
  (`"box"`), Savitzky-Golay smoothing (`"savgol"`), moving average
  (`"mva"`), triangular moving average (`"tmva"`), or no smoothing
  (`"none"`).

- smooth_window:

  Smoothing window. Larger values of this parameter will exclude sharp,
  narrow features. If the supplied value is between 0 and 1, the window
  will be interpreted as a proportion of points to include. Otherwise,
  the window will be interpreted as the absolute number of points to
  include in the window. Defaults to `0.001`.

- slope_thresh:

  Minimum threshold for slope of the smoothed first derivative. This
  parameter filters peaks on the basis of their width, such that larger
  values will exclude broad peaks from the peak list. Defaults to `0`.

- amp_thresh:

  Minimum threshold for peak amplitude. This parameter filters on the
  basis of peak height, such that larger values will exclude small peaks
  from the peak list. Defaults to `0`.

- bounds:

  Logical. If `TRUE` (default), includes peak boundaries in
  `data.frame`.

## Value

If `bounds == TRUE`, returns a `data.frame` containing the center,
start, and end of each identified peak. Otherwise, returns a numeric
vector of peak centers. All locations are expressed as indices.

## Details

Available smoothing methods include Gaussian kernel (`"gaussian"`), box
kernel (`"box"`), Savitzky–Golay (`"savgol"`), moving average (`"mva"`),
triangular moving average (`"tmva"`), and no smoothing (`"none"`).

Preprocessing with
[`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md)
is recommended prior to peak detection. Overly high chromatographic
resolution may sometimes cause peak splitting. In such cases, it is
recommended to either increase the `smooth_window` or reduce the
time-axis resolution by adjusting the `dim1` argument during
preprocessing.

## Note

The `find_peaks` function is adapted from MATLAB code included in Prof.
Tom O'Haver's [Pragmatic Introduction to Signal
Processing](http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm).

## References

O'Haver, Tom. Pragmatic Introduction to Signal Processing: Applications
in scientific measurement. <https://terpconnect.umd.edu/~toh/spectrum/>
(Accessed January, 2022).

## See also

[`fit_peaks`](https://ethanbass.github.io/chromatographR/reference/fit_peaks.md),
[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_pr)
find_peaks(Sa_pr[[1]][,"220"])
#>    pos lower upper
#> 1   28    14    46
#> 2   49    46    52
#> 3   60    53    61
#> 4   69    61    78
#> 5   81    78    85
#> 6   89    85    93
#> 7  104    93   114
#> 8  115   114   118
#> 9  120   118   123
#> 10 131   123   139
#> 11 141   139   147
#> 12 155   149   161
#> 13 167   161   171
#> 14 182   171   188
#> 15 193   188   200
#> 16 215   200   228
#> 17 233   228   238
#> 18 243   238   254
#> 19 261   254   264
#> 20 272   264   279
#> 21 283   279   293
#> 22 297   293   303
#> 23 305   303   309
#> 24 319   315   323
#> 25 336   323   347
#> 26 351   347   363
#> 27 378   363   383
#> 28 388   383   390
#> 29 395   390   401
#> 30 408   401   411
#> 31 417   411   421
```
