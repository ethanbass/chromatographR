# Preprocess time/wavelength data

Standard pre-processing of response matrices, consisting of a time axis
and a spectral axis (e.g. HPLC-DAD/UV data). For smooth data, like
UV-VIS data, the size of the matrix can be reduced by interpolation. By
default, the data are baseline-corrected in the time direction
([`baseline.corr`](https://rdrr.io/pkg/ptw/man/baseline.corr.html)) and
smoothed in the spectral dimension using cubic smoothing splines
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
  ...
)
```

## Arguments

- X:

  A numerical data matrix, or list of data matrices. Missing values are
  not allowed. If rownames or colnames attributes are used, they should
  be numerical and signify time points and wavelengths, respectively.

- dim1:

  A new, usually shorter, set of time points (numerical). The range of
  these should not exceed the range of the original time points.

- dim2:

  A new, usually shorter, set of wavelengths (numerical). The range of
  these should not exceed the range of the original wavelengths.

- remove.time.baseline:

  Logical, indicating whether baseline correction should be done in the
  time direction, according to
  [`baseline.corr`](https://rdrr.io/pkg/ptw/man/baseline.corr.html).
  Default is `TRUE`.

- spec.smooth:

  Logical, indicating whether smoothing should be done in the spectral
  direction, according to
  [`smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html). Default
  is `TRUE`.

- maxI:

  if given, the maximum intensity in the matrix is set to this value.

- interpolate_rows:

  Logical. Whether to interpolate along the time axis (`dim1`). Defaults
  to `TRUE`.

- interpolate_cols:

  Logical. Whether to interpolate along the spectral axis (`dim2`).
  Defaults to `TRUE`.

- cl:

  Argument to
  [`pblapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  or [`mclapply`](https://rdrr.io/r/parallel/mclapply.html). Either an
  integer specifying the number of clusters to use for parallel
  processing or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to 2. On Windows integer values will be ignored.

- show_progress:

  Logical. Whether to show progress bar. Defaults to `TRUE` if
  [`pbapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  is installed.

- ...:

  Further optional arguments to
  [`baseline.corr`](https://rdrr.io/pkg/ptw/man/baseline.corr.html).

## Value

The function returns the preprocessed data matrix (or list of matrices),
with row names and column names indicating the time points and
wavelengths, respectively.

## Note

Adapted from the
[preprocess](https://github.com/rwehrens/alsace/blob/master/R/preprocess.R)
function in the alsace package by Ron Wehrens.

## References

- Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast parametric
  time warping of peak lists. *Bioinformatics* **31**:3063-3065.
  [doi:10.1093/bioinformatics/btv299](https://doi.org/10.1093/bioinformatics/btv299)
  .

- Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
  LC–DAD using multivariate curve resolution: the alsace package for R.
  *Metabolomics* **11:1**:143-154.
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
