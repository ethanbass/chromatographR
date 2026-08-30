# Preprocess time/wavelength data

Standard preprocessing of time × wavelength response matrices (e.g.
HPLC-DAD/UV data), consisting of: (i) baseline correction along the time
axis, (ii) smoothing along the spectral axis, and (iii) optional
interpolation in either dimension for dimensionality reduction. For
densely sampled data, (e.g., UV-VIS spectra, the size of the matrix can
be reduced by interpolation. By default, the data are baseline-corrected
along the time axis using
[`baseline.corr`](https://rdrr.io/pkg/ptw/man/baseline.corr.html), and
smoothed along the spectral axis using cubic smoothing splines
([`smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html)).

## Usage

``` r
preprocess(
  X,
  dim1,
  dim2,
  remove.time.baseline = TRUE,
  spec.smooth = TRUE,
  maxI = NULL,
  interpolate_rows = TRUE,
  interpolate_cols = TRUE,
  cl = 2,
  show_progress = NULL,
  outlier_cutoff = 5/60,
  ...
)
```

## Arguments

- X:

  A numeric matrix (time × wavelength) or a list of such matrices. Row
  names must correspond to time points and column names to wavelengths.

- dim1:

  Numeric vector specifying the target time grid. The data will be
  interpolated onto these time points. The range of the new values
  should not exceed the range of the original time points.

- dim2:

  Numeric vector specifying the target wavelength grid. The data will be
  interpolated onto these wavelengths. The range of the new values
  should not exceed the range of the original wavelengths.

- remove.time.baseline:

  Logical, indicating whether baseline correction should be done along
  the time axis, according to
  [`baseline.corr`](https://rdrr.io/pkg/ptw/man/baseline.corr.html).
  Defaults to `TRUE`.

- spec.smooth:

  Logical, indicating whether smoothing should be done along the
  spectral axis, according to
  [`smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html).
  Defaults to `TRUE`.

- maxI:

  If supplied, all values are rescaled so that the maximum intensity
  equals `maxI`.

- interpolate_rows:

  Logical. Whether to interpolate along the time axis (`dim1`). Defaults
  to `TRUE`.

- interpolate_cols:

  Logical. Whether to interpolate along the spectral axis (`dim2`).
  Defaults to `TRUE`.

- cl:

  Either an integer specifying the number of cores to use for parallel
  processing or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to `2`. On Windows integer values will be ignored.

- show_progress:

  Logical. Whether to show progress bar. Defaults to `TRUE` if
  [`pbapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  is installed.

- outlier_cutoff:

  Threshold (in seconds) for excluding chromatograms that end earlier
  than expected. Samples ending more than this value before the median
  end time are removed. Defaults to `5` seconds. Only applies when
  `dim1` is not specified.

- ...:

  Additional arguments to
  [`baseline.corr`](https://rdrr.io/pkg/ptw/man/baseline.corr.html).

## Value

The function returns the preprocessed data matrix (or list of matrices),
with row names and column names indicating the time points and
wavelengths, respectively.

## Note

Adapted from the `preprocess` function in the alsace package by Ron
Wehrens:
<https://github.com/rwehrens/alsace/blob/master/R/preprocess.R>.

## References

- Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast parametric
  time warping of peak lists. *Bioinformatics* 31:3063-3065.
  [doi:10.1093/bioinformatics/btv299](https://doi.org/10.1093/bioinformatics/btv299)
  .

- Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
  LC–DAD using multivariate curve resolution: the alsace package for R.
  *Metabolomics* 11:1:143-154.
  [doi:10.1007/s11306-014-0683-5](https://doi.org/10.1007/s11306-014-0683-5)
  .

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa)
new.ts <- seq(10,18.66,by=.01) # choose time-points
new.lambdas <- seq(200, 318, by = 2) # choose wavelengths
Sa_pr <- preprocess(Sa[[1]], dim1 = new.ts, dim2 = new.lambdas)
}
```
