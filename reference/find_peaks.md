# Find peaks

Detects peaks in a chromatographic profile by finding zero-crossings in
the smoothed first derivative of the signal (`y`). Smoothing suppresses
spurious local extrema. Peaks can be filtered by minimum curvature
(`min_curvature`), minimum amplitude (`min_amplitude`), or minimum
prominence relative to the baseline noise (`min_prominence_snr`).

## Usage

``` r
find_peaks(
  y,
  smooth_type = c("gaussian", "box", "savgol", "mva", "tmva", "none"),
  smooth_window = 3,
  min_curvature = 0,
  min_amplitude = 0,
  min_prominence_snr = 0,
  bounds = TRUE,
  slope_thresh = NULL,
  amp_thresh = NULL
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

  Width of the smoothing window. A value between 0 and 1 is interpreted
  as a proportion of the length of `y`; a value of 1 or greater as a
  number of timepoints. Defaults to `3`. See *Choosing the smoothing
  window*.

- min_curvature:

  Minimum curvature at the apex, as the magnitude of the second
  derivative of the smoothed signal. Curvature scales as height over
  width squared, so larger values exclude broader peaks. Defaults to
  `0`.

- min_amplitude:

  Minimum peak amplitude, as the intensity of the unsmoothed signal at
  the apex. Apices at or below this value are removed. Defaults to `0`,
  which removes apices at or below zero.

- min_prominence_snr:

  Minimum peak prominence, as a multiple of the noise standard
  deviation. Requires `bounds = TRUE`. Only meaningful on
  baseline-corrected data (see *Filtering peaks*). Defaults to `0`,
  which applies no filtering.

- bounds:

  Logical. If `TRUE` (default), includes peak boundaries in a
  `data.frame`.

- slope_thresh:

  The `slope_thresh` argument is deprecated. Please use `min_curvature`
  instead.

- amp_thresh:

  The `amp_thresh` argument is deprecated. Please use `min_amplitude`
  instead.

## Value

If `bounds == TRUE`, a `data.frame` with one row per peak and columns
`pos` (apex), `lower` and `upper` (boundaries), and `prominence`.
Otherwise, a numeric vector of apex positions. All locations are
expressed as indices.

## Details

How many peaks are detected depends mostly on `smooth_window` (see
*Choosing the smoothing window*). Three filters, `min_curvature`,
`min_amplitude` and `min_prominence_snr`, can then remove small or
spurious peaks (see *Filtering peaks*).

## Note

The `find_peaks` function is adapted from MATLAB code included in Prof.
Tom O'Haver's [Pragmatic Introduction to Signal
Processing](https://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm).

## Choosing the smoothing window

Peaks are detected where the first derivative of the signal changes from
positive to negative. Noise makes the derivative change sign many times
on its own, so the signal is first smoothed over a window set by the
`smooth_window` argument. A wider window suppresses more noise, but can
also blur neighboring peaks together if the window is too wide.

Set `smooth_window` to about a third of the full width at half maximum
of a typical peak, in timepoints. The default, `3`, suits peaks about 10
timepoints wide. Wider windows may merge neighboring peaks. This may be
easy to miss through casual inspection of peak tables — a merged pair is
reported as one peak, with a single retention time and their combined
area. Plot the peak list over the trace to check. Narrower windows add
spurious apices. At fine time resolution, you can also reduce the
resolution with the `dim1` argument to
[`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md).

A `smooth_window` between 0 and 1 is interpreted as a proportion of the
run. A proportional smoothing window is invariant to sampling rate, but
it scales with run length: `0.001` is only 3 timepoints for a
3,000-point trace but 100 for a 100,000-point one. A width supplied in
timepoints has the opposite trade-off: it doesn't depend on run length,
but it scales with the sampling rate, so if you change the time grid
with `preprocess(dim1 = )`, the `smooth_window` should be scaled
accordingly.

Under-smoothing also undermines `min_prominence_snr`. It can split a
small peak into two apices, each measured from the shallow valley
between them, so the filter discards both as noise and the peak is lost.
Set the smoothing window wide enough to avoid splitting peaks before
applying prominence filtering. `find_peaks` warns when `smooth_window`
is under a third of the typical peak width.

The [peak-finding
article](https://ethanbass.github.io/chromatographR/articles/peak-finding.html)
tests these recommendations on simulated chromatograms with known peaks.

**Method-specific caveats.** At the same `smooth_window`, `"savgol"`
smooths about a third as much as the other methods, so it needs a wider
window for the same noise reduction. It errors below 7 timepoints.
`"savgol"`, `"mva"`, and `"tmva"` round the window to the nearest
integer, then up to the next odd one if it is even. The moving averages
also use at least three points, since a one-point window does no
smoothing.

## Filtering peaks

The curvature and amplitude filters (`min_curvature` and
`min_amplitude`) are in absolute instrument units, so a usable value
must be chosen for each dataset. The prominence filter
(`min_prominence_snr`) instead measures each peak against the noise. A
peak's prominence is its height above the higher of the two valleys on
either side of it. For an isolated peak this is its full height; when a
neighboring peak raises one of the valleys, it is less. The noise is
estimated iteratively across the whole chromatogram. Starting with a
rough estimate from the change in signal between adjacent timepoints,
`find_peaks` excludes every peak whose apex is more than three times the
current estimate above the median of the remaining timepoints, then
re-estimates the noise as the scaled median absolute deviation of the
remaining timepoints. This process is repeated until the estimate
changes by less than 0.5%. `min_prominence_snr = 4` keeps peaks with
prominence exceeding 4x the estimated noise standard deviation.

Because the remaining timepoints include baseline drift, the estimated
standard deviation can reflect drift rather than noise, so
`min_prominence_snr` is only meaningful on baseline-corrected data.

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
#>    pos lower upper  prominence
#> 1   29    14    48   9.0254038
#> 2   49    48    54   0.1546514
#> 3   70    54    79 125.1386438
#> 4   82    79    86   2.0772941
#> 5   90    86    93   2.0177936
#> 6  104    93   124 876.5683456
#> 7  132   124   149 471.3078495
#> 8  156   149   161   7.3412163
#> 9  168   161   173  27.5675529
#> 10 183   173   189 386.7094028
#> 11 194   189   201  98.8385660
#> 12 216   201   229 421.3882160
#> 13 234   229   240  40.4113073
#> 14 244   240   255   8.5081663
#> 15 262   255   265  13.6675782
#> 16 273   265   279  81.8281072
#> 17 284   279   294 139.1889566
#> 18 298   294   311  13.9836910
#> 19 320   311   323   0.8434209
#> 20 337   323   348 174.3813353
#> 21 352   348   364   4.7501747
#> 22 379   364   384   0.4880988
#> 23 396   384   403   1.4539185
#> 24 409   403   412   0.1748946
#> 25 418   412   422   1.0162429
```
