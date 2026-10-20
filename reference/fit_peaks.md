# Fit chromatographic peaks

Peak parameters are estimated by fitting one of several models to the
data using nonlinear least squares
([`nlsLM`](https://rdrr.io/pkg/minpack.lm/man/nlsLM.html)). Supported
models include Gaussian, bidirectional exponentially modified Gaussian,
or exponential-Gaussian hybrid functions. The area under each fitted
curve is estimated using trapezoidal approximation.

## Usage

``` r
fit_peaks(
  x,
  lambda,
  pos = NULL,
  sd_max = NULL,
  fit = c("egh", "bemg", "gaussian", "raw"),
  max_iter = 1000,
  estimate_purity = FALSE,
  purity_args = list(),
  noise_variance = NULL,
  baseline = c("none", "flat", "sloped"),
  ...
)
```

## Arguments

- x:

  A single chromatogram in matrix format.

- lambda:

  Wavelength at which to fit peaks. Must correspond to a wavelength
  present in `x`.

- pos:

  A `data.frame` of peak apices and boundaries, as returned by
  [`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md),
  with columns `pos`, `lower` and `upper` giving row indices of `x`. If
  `NULL` (the default),
  [`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md)
  is run on `x` at `lambda`.

- sd_max:

  Deprecated. Filter peaks by width with
  [`filter_peaks`](https://ethanbass.github.io/chromatographR/reference/filter_peaks.md)
  instead. If supplied, it is applied as a cutoff on `FWHM` of
  `2.355 * sd_max` timepoints.

- fit:

  Peak model. Options are bidirectional exponentially modified Gaussian
  (`"bemg"`), exponential-Gaussian hybrid (`"egh"`), `"gaussian"`, or
  `"raw"`. `"raw"` integrates the signal directly by trapezoidal
  approximation without fitting a model. Defaults to `"egh"`. See *Peak
  models* in
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md).

- max_iter:

  Maximum number of iterations to use in nonlinear least squares
  peak-fitting. Defaults to `1000`.

- estimate_purity:

  Logical. Whether to estimate purity. Defaults to `FALSE`, and is
  always `FALSE` when `x` holds a single wavelength.

- purity_args:

  A named list of additional arguments passed to
  [`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md),
  e.g. `list(background_correction = "none", metric = "max")`. A
  `lambdas` element here selects the wavelengths used to compare
  spectra, not the wavelengths at which peaks are integrated.

- noise_variance:

  Variance of noise, passed to
  [`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md).
  If not provided, it is estimated once for the whole chromatogram
  rather than once per peak.

- baseline:

  Whether to measure each peak above a local baseline: `"none"` (the
  default), `"flat"` or `"sloped"`. See *Local baselines* in
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md).

- ...:

  Additional arguments passed to
  [`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md),
  such as `smooth_window`, `smooth_type`, `min_prominence_snr`,
  `min_curvature` and `min_amplitude`. See *Choosing the smoothing
  window* in
  [`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md).

## Value

A `data.frame` with one row per peak and the following columns:

- `rt`: Apex of the fitted peak model, or of the detected peak when
  `fit = "raw"`.

- `start`: Peak start.

- `end`: Peak end.

- `sd`: Standard deviation of the Gaussian component of the fitted
  model, which is not the width of the peak unless `fit = "gaussian"`.
  `NA` when `fit = "raw"`.

- `tau`: Exponential decay constant (only when `fit = "egh"`).

- `tau_right`: Time constant of the exponential tail on the right side
  of the peak, the reciprocal of its fitted rate (only when
  `fit = "bemg"`).

- `tau_left`: Time constant of the exponential tail on the left side of
  the peak, the reciprocal of its fitted rate (only when
  `fit = "bemg"`).

- `h`: The amplitude parameter of the fitted BEMG shape function. It
  scales the shape but is neither its value at `center` nor the peak
  `height` (only when `fit = "bemg"`).

- `center`: The center parameter of the fitted BEMG shape function, the
  location around which the exponential tails are applied. May differ
  from `rt` (only when `fit = "bemg"`).

- `FWHM`: The full width at half maximum of the peak, comparable across
  values of `fit`. `NA` for a `"raw"` peak whose integration window does
  not fall below half its height on both sides.

- `height`: Peak height above the fitted baseline.

- `area`: Peak area above the fitted baseline.

- `r.squared`: Coefficient of determination (\\R^2\\) of the fitted
  curve over the peak window, `1 - SS_res / SS_tot`. Negative when the
  fit is worse than a flat line at the mean. Not returned when
  `fit = "raw"`.

- `purity`: The spectral purity of peak as assessed by
  [`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md).

- `window_limited`: `1` if a fitted shape parameter hit a limit set by
  the integration window, or a `"gaussian"` or `"egh"` fit narrowed to
  one sampling interval at half height, meaning the window or the
  sampling rather than the data determined the reported shape, and `0`
  otherwise. Not returned when `fit = "raw"`.

- `converged`: `1` if the optimizer reported convergence for the final
  step of the fit, and `0` otherwise, including when part or all of the
  fit was skipped or failed (see *Fitting procedure* in
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)).
  Not returned when `fit = "raw"`.

- `floor`: Constant baseline offset (only when `baseline = "flat"`).

- `floor_start`: Baseline level at the left edge of the peak window
  (only when `baseline = "sloped"`).

- `floor_end`: Baseline level at the right edge of the peak window (only
  when `baseline = "sloped"`). The table is in index units rather than
  retention times: `rt`, `start`, `end` and `center` are index
  positions; `sd`, `FWHM`, `tau`, `tau_right` and `tau_left` are widths
  on the index scale; and `area` is integrated over the index axis.
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)
  converts them to the chromatogram's time axis.

## Details

Peaks that can't be fitted are dropped if they've been detected
automatically with `find_peaks`. Otherwise, if the peak positions are
provided explicitly through `pos`, `height` and `area` are set to `NA`
and they are marked `0` in the `converged` column.

The peak models, local baselines and fitting procedure are described in
[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md),
under *Peak models*, *Local baselines* and *Fitting procedure*.

## Note

The loop over peaks in `fit_peaks` is adapted from the `fitpeaks`
function in Ron Wehrens' [alsace](https://github.com/rwehrens/alsace)
package. The Gaussian peak shape and the interface of the internal
fitting functions are adapted from Robert Morrison's
[DuffyTools](https://github.com/robertdouglasmorrison/DuffyTools)
package.

## References

- Lan, K. & Jorgenson, J. W. 2001. A hybrid of exponential and gaussian
  functions as a simple model of asymmetric chromatographic peaks.
  *Journal of Chromatography A* 915:1-13.
  [doi:10.1016/S0021-9673(01)00594-5](https://doi.org/10.1016/S0021-9673%2801%2900594-5)
  .

- Naish, P. J. & Hartwell, S. 1988. Exponentially Modified Gaussian
  functions - A good model for chromatographic peaks in isocratic HPLC?
  *Chromatographia*, 26: 285-296.
  [doi:10.1007/BF02268168](https://doi.org/10.1007/BF02268168) .

## See also

[`find_peaks`](https://ethanbass.github.io/chromatographR/reference/find_peaks.md),
[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_pr)
fit_peaks(Sa_pr[[1]], lambda = 220)
#>           rt start end        sd          tau      FWHM      height
#> 1   30.66114    14  48  6.101575   1.44697364 14.403076   8.9748818
#> 2   49.94050    48  54  1.181971  -1.10373785  2.886559   0.3915597
#> 3   70.01920    54  79  2.729871  -1.17486581  6.479732 123.7545250
#> 4   82.33054    79  86  1.482960  -0.43913100  3.505344   2.5630020
#> 5   89.65612    86  93  1.725862   0.74138982  4.096455   2.8793770
#> 6  104.49705    93 124  2.183908  -0.19404617  5.144468 894.3629580
#> 7  131.85646   124 149  1.833434   0.11605405  4.318157 468.3178181
#> 8  156.03382   149 161  1.858660   0.17723271  4.378535   7.4302224
#> 9  167.74077   161 173  1.681790  -0.15232333  3.961719  28.6002885
#> 10 183.20931   173 189  2.094857  -0.16099380  4.934272 414.2717933
#> 11 193.90505   189 201  2.400072  -0.87407104  5.684120 128.4663742
#> 12 216.44076   201 229  1.842958  -0.13155394  4.340793 433.8264796
#> 13 234.37947   229 240  2.087220   0.30929995  4.919702  47.6597945
#> 14 243.28110   240 255  3.721672   3.76517021  9.144208  16.0564058
#> 15 261.90773   255 265  2.102433  -0.57181758  4.966692  18.9817970
#> 16 272.86037   265 279  3.095311   0.09729706  7.289213  94.6165502
#> 17 284.46078   279 294  1.877091  -0.35653644  4.427114 161.9127876
#> 18 298.16377   294 311  1.928674  -0.35565627  4.548366  17.6025746
#> 19 320.10828   311 323  1.046951   0.27961372  2.472988   0.7954206
#> 20 336.89987   323 348  2.228935  -0.81884813  5.279340 171.3788902
#> 21 351.76261   348 364  3.699194   3.37989773  9.020474   6.9747887
#> 22 378.91468   364 384 10.262991  -7.88194072 24.777329   3.4392483
#> 23 396.45231   384 403  9.840007  -2.54701213 23.238604   5.8916446
#> 24 411.06601   403 412  7.643896 -97.75200512 70.106682   5.3320058
#> 25 419.24718   412 422  6.088367  -6.96651196 15.128362   7.2786470
#>            area r.squared purity window_limited converged
#> 1   136.6168187 0.9834201     NA              0         1
#> 2     0.9895294 0.1434314     NA              0         1
#> 3   861.2327696 0.9910190     NA              0         1
#> 4     9.2206621 0.9982737     NA              0         1
#> 5    11.6265191 0.9937378     NA              0         1
#> 6  4899.5751670 0.9985326     NA              0         1
#> 7  2153.0681287 0.9981330     NA              0         1
#> 8    34.3806776 0.9355800     NA              0         1
#> 9   120.5844854 0.9995746     NA              0         1
#> 10 2172.4932812 0.9960212     NA              0         1
#> 11  735.7342180 0.9970448     NA              0         1
#> 12 2005.0667012 0.9957840     NA              0         1
#> 13  247.0319456 0.9661289     NA              0         1
#> 14  145.4116996 0.9588278     NA              0         1
#> 15   96.0671652 0.9711060     NA              0         1
#> 16  710.5446283 0.8959060     NA              0         1
#> 17  758.8991868 0.9748398     NA              0         1
#> 18   82.6349173 0.9781477     NA              0         1
#> 19    2.0678877 0.1506531     NA              0         1
#> 20  969.4820324 0.9942008     NA              0         1
#> 21   64.0434095 0.8251971     NA              0         1
#> 22   57.4508700 0.9935989     NA              0         1
#> 23   94.7135320 0.8801906     NA              0         1
#> 24   46.5196604 0.8663922     NA              1         1
#> 25   64.5439100 0.9002292     NA              0         1
```
