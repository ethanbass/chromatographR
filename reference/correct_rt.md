# Correct retention time

Aligns chromatograms using one of two algorithms, according to the value
of `alg`: either parametric time warping, as implemented in
[`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html), or variable penalty
dynamic time warping, as implemented in
[`VPdtw`](https://rdrr.io/pkg/VPdtw/man/VPdtw.html). The `init.coef` and
`n.traces` arguments apply only to `ptw` warping, while `penalty` and
`maxshift` apply only to `vpdtw` warping.

## Usage

``` r
correct_rt(
  chrom_list,
  lambdas,
  models = NULL,
  reference = "best",
  alg = c("ptw", "vpdtw"),
  what = c("corrected.values", "models"),
  init.coef = c(0, 1, 0),
  n.traces = NULL,
  n.zeros = 0,
  scale = FALSE,
  trwdth = 200,
  plot_it = FALSE,
  penalty = 5,
  maxshift = 50,
  verbose = getOption("verbose"),
  show_progress = NULL,
  cl = 2,
  ...
)
```

## Arguments

- chrom_list:

  List of chromatograms in matrix format.

- lambdas:

  Select wavelengths to use by name.

- models:

  List of models to warp by. The models provided here (if any) must
  match the algorithm selected in `alg`.

- reference:

  Index of the sample that is to be considered the reference sample.

- alg:

  algorithm to use: parametric time warping (`ptw`) or variable penalty
  dynamic time warping (`vpdtw`).

- what:

  What to return: either the 'corrected.values' (useful for visual
  inspection) or the warping 'models' (for further programmatic use).

- init.coef:

  Starting values for the optimization.

- n.traces:

  Number of traces to use.

- n.zeros:

  Number of zeros to add.

- scale:

  Logical. If true, scale chromatograms before warping.

- trwdth:

  width of the triangle in the WCC criterion.

- plot_it:

  Logical. Whether to plot alignment.

- penalty:

  The divisor used to calculate the penalty for
  [`VPdtw`](https://rdrr.io/pkg/VPdtw/man/VPdtw.html). The warping
  penalty is calculated by dividing the
  [`dilation`](https://rdrr.io/pkg/VPdtw/man/dilation.html) by this
  number. Thus, a higher number will produce a lower penalty and be more
  permissive, while a lower number will produce a higher penalty and
  allow less warping. Defaults to 5.

- maxshift:

  Integer. Maximum allowable shift for
  [`VPdtw`](https://rdrr.io/pkg/VPdtw/man/VPdtw.html). Defaults to 50.

- verbose:

  Whether to print verbose output.

- show_progress:

  Logical. Whether to show progress bar. Defaults to `TRUE` if
  [`pbapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  is installed. Currently works only for `ptw` alignments.

- cl:

  Argument to
  [`pblapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  or [`mclapply`](https://rdrr.io/r/parallel/mclapply.html). Either an
  integer specifying the number of clusters to use for parallel
  processing or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to 2. On Windows integer values will be ignored.

- ...:

  Optional arguments for the
  [`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html) function. The only
  argument that cannot be changed is `warp.type`: this is always equal
  to `"global"`.

## Value

A list of warping models or a list of warped absorbance profiles,
depending on the value of the `what` argument.

## Note

Adapted from
[correctRT](https://github.com/rwehrens/alsace/blob/master/R/correctRT.R)
function in the alsace package by Ron Wehrens.

## References

- Clifford, D., Stone, G., Montoliu, I., Rezzi, S., Martin, F. P., Guy,
  P., Bruce, S., & Kochhar, S. 2009. Alignment using variable penalty
  dynamic time warping. *Analytical chemistry*, **81(3)**:1000-1007.
  [doi:10.1021/ac802041e](https://doi.org/10.1021/ac802041e) .

- Clifford, D., & Stone, G. 2012. Variable Penalty Dynamic Time Warping
  Code for Aligning Mass Spectrometry Chromatograms in R. *Journal of
  Statistical Software*, **47(8)**:1-17.
  [doi:10.18637/jss.v047.i08](https://doi.org/10.18637/jss.v047.i08) .

- Eilers, P.H.C. 2004. Parametric Time Warping. *Anal. Chem.*,
  **76**:404-411.
  [doi:10.1021/ac034800e](https://doi.org/10.1021/ac034800e) .

- Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast parametric
  time warping of peak lists. *Bioinformatics*, **31**:3063-3065.
  [doi:10.1093/bioinformatics/btv299](https://doi.org/10.1093/bioinformatics/btv299)
  .

- Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
  LC–DAD using multivariate curve resolution: the alsace package for R.
  *Metabolomics*, **11**:143-154.
  [doi:10.1007/s11306-014-0683-5](https://doi.org/10.1007/s11306-014-0683-5)
  .

## See also

[`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html),
[`correct_peaks`](https://ethanbass.github.io/chromatographR/reference/correct_peaks.md),
[`VPdtw`](https://rdrr.io/pkg/VPdtw/man/VPdtw.html)

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa_pr)
warp <- correct_rt(chrom_list = Sa_pr, lambdas=210)
}
```
