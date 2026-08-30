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
  sd_max = 50,
  fit = c("egh", "bemg", "gaussian", "raw"),
  max_iter = 1000,
  estimate_purity = TRUE,
  noise_threshold = 0.001,
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

  Locations of peaks in vector `y`. If `NULL`, `find_peaks` will run
  automatically to find peak positions.

- sd_max:

  Maximum width (standard deviation) for peaks. Defaults to `50`.

- fit:

  Peak model to use. Options include `"bemg"` (bidirectional
  exponentially modified Gaussian), `"egh"` (exponential-Gaussian
  hybrid), `"gaussian"`, and `"raw"`. If `raw` is selected, peaks are
  integrated using trapezoidal integration without model fitting.
  Defaults to `"egh"`.

- max_iter:

  Maximum number of iterations to use in nonlinear least squares
  peak-fitting. Defaults to `1000`.

- estimate_purity:

  Logical. Whether to estimate purity or not. Defaults to `TRUE`.

- noise_threshold:

  Noise threshold. Input to `get_purity`.

- baseline:

  Whether to fit a baseline offset beneath each peak model. Options are
  `"none"` (no baseline, default), `"flat"` (constant offset), or
  `"sloped"` (linearly varying offset). Not available when
  `fit = "raw"`.

- ...:

  Additional arguments to `find_peaks`.

## Value

A matrix with one row per peak and the following columns:

- `rt`: Peak maximum location.

- `start`: Peak start (only included in table if `bounds = TRUE`).

- `end`: Peak end (only included in table if `bounds = TRUE`).

- `sd`: Peak standard deviation.

- `tau`: Exponential decay constant (only included in table if
  `fit = "egh"`).

- `tau_right`: Exponential decay constant for right side of peak (when
  `bemg` is selected).

- `tau_left`: Exponential decay constant for left side of peak (when
  `bemg` is selected).

- `FWHM`: The full width at half maximum.

- `height`: Peak height above the fitted baseline.

- `area`: Peak area above the fitted baseline.

- `r.squared`: The R² value for linear fit of the model to the data.

- `purity`: The spectral purity of peak as assessed by
  [`get_purity`](https://ethanbass.github.io/chromatographR/reference/get_purity.md).

- `floor`: Constant baseline offset (only when `baseline = "flat"`).

- `floor_start`: Baseline level at the left edge of the peak window
  (only when `baseline = "sloped"`).

- `floor_end`: Baseline level at the right edge of the peak window (only
  when `baseline = "sloped"`). The first five elements (`rt`, `start`,
  `end`, `sd` and `FWHM`) are expressed as index positions rather than
  absolute retention times. The transformation to real time is done in
  `get_peaks`.

## Note

The `fit_peaks` function is adapted from Dr. Robert Morrison's
[DuffyTools](https://github.com/robertdouglasmorrison/DuffyTools)
package as well as code published in Ron Wehrens'
[alsace](https://github.com/rwehrens/alsace) package.

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
#>           rt start end        sd           tau      FWHM      height
#> 1   30.69332    14  46 6.2865644    1.20769992 14.773426   9.1034627
#> 2   50.94330    46  52 0.5711967  -11.41366848  1.342312   0.4022556
#> 3   61.44532    53  61 3.2042222    2.37308275  7.529922  14.3385995
#> 4   69.97526    61  78 2.7067311   -1.01879444  6.360818 124.1300652
#> 5   82.42709    78  85 1.3744994   -0.73109961  3.230074   2.3051668
#> 6   89.64284    85  93 1.5296385    0.61202488  3.594650   2.6256571
#> 7  104.49707    93 114 2.1826134   -0.19394678  5.129142 894.1248943
#> 8   81.84177   114 118 2.4690722   22.61770113  5.802320  42.5947951
#> 9  119.85327   118 123 5.2832737  -13.86056689 12.415693   8.6545653
#> 10 131.85675   123 139 1.8315236    0.11513594  4.304081 468.1212567
#> 11 141.10828   139 147 4.1299082   -4.37157201  9.705284   8.4526758
#> 12 156.09501   149 161 1.7963783   -0.01487998  4.221489   7.1019209
#> 13 167.73853   161 171 1.6517982   -0.14556475  3.881726  28.2608430
#> 14 183.22815   171 188 2.0807233   -0.22145455  4.889700 414.6928515
#> 15 193.98654   188 200 2.4682084   -1.29245512  5.800290 126.0967365
#> 16 216.44045   200 228 1.8375563   -0.13007408  4.318257 433.2796803
#> 17 234.50440   228 238 1.9619252   -0.17648212  4.610524  47.8442003
#> 18 243.58113   238 254 4.9483859    0.66568143 11.628707  14.0416866
#> 19 262.00099   254 264 1.8302984   -0.98597515  4.301201  19.1541885
#> 20 272.87432   264 279 3.0837104    0.05107119  7.246719  94.1813171
#> 21 284.45995   279 293 1.8712197   -0.34947351  4.397366 161.7599016
#> 22 298.28530   293 303 2.0217365   -0.90317976  4.751081  17.1078695
#> 24 320.15708   315 323 1.2587749    0.34042184  2.958121   0.9339542
#> 25 336.89949   323 347 2.2339189   -0.82019092  5.249709 171.5472176
#> 26 351.82167   347 363 4.3854006    2.73641711 10.305691   6.9880680
#> 27 379.01651   363 383 9.6193330  -13.36671566 22.605432   3.8160109
#> 28 389.96461   383 390 0.4505398  -11.19848606  1.058768   5.7371410
#> 29 396.32753   390 401 7.3086188   -2.14825864 17.175254   6.4654057
#> 30 410.94947   401 411 1.7808715 -123.44165218  4.185048   5.7088185
#> 31 419.42472   411 421 4.8293966  -14.83514636 11.349082   7.7061414
#>           area r.squared     purity
#> 1   141.237161 0.9869708 0.63636364
#> 2     1.664893 0.7267476 1.00000000
#> 3    37.107139 0.9997068 1.00000000
#> 4   846.567716 0.9954198 0.17647059
#> 5     7.964953 0.9529150 1.00000000
#> 6     9.667361 0.9974965 1.00000000
#> 7  4895.364758 0.9991512 0.08333333
#> 8    38.579138 0.8292605 0.60000000
#> 9    38.911532 0.9637970 1.00000000
#> 10 2149.558127 0.9991802 0.15384615
#> 11   46.615211 0.9910774 1.00000000
#> 12   31.861289 0.9960339 1.00000000
#> 13  114.667078 0.9997472 0.90909091
#> 14 2149.492159 0.9987345 0.15384615
#> 15  759.917479 0.9870224 0.18181818
#> 16 1996.653868 0.9988400 0.50000000
#> 17  227.707869 0.9889112 0.72727273
#> 18  150.579316 0.8325459 0.52941176
#> 19   84.009309 0.9912112 0.63636364
#> 20  707.817033 0.9118992 0.46666667
#> 21  755.914103 0.9843055 0.26666667
#> 22   85.445828 0.9455491 0.72727273
#> 24    2.857647 0.9897500 1.00000000
#> 25  972.589353 0.9950354 0.33333333
#> 26   71.516308 0.7701393 1.00000000
#> 27   63.064242 0.9955775 1.00000000
#> 28   29.703398 0.9392470 1.00000000
#> 29   64.730873 0.8411157 1.00000000
#> 30   54.811944 0.8662737 1.00000000
#> 31   66.589250 0.8950062 1.00000000
```
