# Retention time alignment via time warping

Chromatograms are aligned using either parametric time warping
(`"ptw"`), as implemented in
[`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html) or variable penalty
dynamic time warping (`"vpdtw"`), as implemented in
[`VPdtw`](https://rdrr.io/pkg/VPdtw/man/VPdtw.html).

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
  fill_zeros = FALSE,
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

  A character or numeric vector specifying the wavelengths to use for
  alignment. Only one wavelength should be specified for VPdtw warping.
  For one-dimensional chromatograms, this argument can be ignored.

- models:

  List of models to warp by. The models provided here (if any) must
  match the algorithm selected in `alg`.

- reference:

  Index of the sample to be used as the reference. If no reference is
  specified, the reference will be chosen algorithmically from a
  similarity matrix of the supplied chromatograms using the
  [`bestref`](https://rdrr.io/pkg/ptw/man/bestref.html) function from
  `ptw`.

- alg:

  Alignment algorithm to use: parametric time warping (`"ptw"`), or
  variable penalty dynamic time warping (`"vpdtw"`).

- what:

  Output type: either the `"corrected.values"` (useful for visual
  inspection and downstream analysis) or the warping `"models"` (for
  further programmatic use).

- init.coef:

  Starting values for the optimization.

- n.traces:

  Number of traces to use.

- fill_zeros:

  Logical. If `TRUE`, out-of-bounds regions produced by warping are
  filled with zeros; otherwise they are returned as `NA` (default).

- n.zeros:

  Number of zeros to add for padding chromatograms at the edges.

- scale:

  Logical. If `TRUE`, scale chromatograms before warping.

- trwdth:

  Argument to [`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html). Width of
  the triangle in the WCC criterion. Defaults to `200`.

- plot_it:

  Logical. Whether to plot alignment. Defaults to `FALSE`.

- penalty:

  The divisor used to calculate the penalty for
  [`VPdtw`](https://rdrr.io/pkg/VPdtw/man/VPdtw.html). The warping
  penalty is calculated by dividing the
  [`dilation`](https://rdrr.io/pkg/VPdtw/man/dilation.html) by this
  number. Thus, a higher number will produce a lower penalty and be more
  permissive, while a lower number will produce a higher penalty and
  allow less warping. Defaults to `5`.

- maxshift:

  Integer. Maximum allowable shift for `VPdtw` warping. Defaults to
  `50`.

- verbose:

  Logical. Whether to print verbose output.

- show_progress:

  Logical. Whether to show progress bar. Defaults to `TRUE` if
  [`pbapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  is installed. Currently works only for `ptw` alignments.

- cl:

  Either an integer specifying the number of cores to use for parallel
  processing or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to `2`. On Windows systems, integer values will be ignored.

- ...:

  Optional additional arguments to `ptw`. The only argument that cannot
  be changed is `warp.type` which is hard-coded to `"global"` to permit
  warping on multiple wavelengths.

## Value

A list of warping models or a list of warped absorbance profiles,
according to the value of the `what` argument.

## Details

Some arguments are specific to particular warping functions. For example
the `init.coef` and `n.traces` arguments apply only to `"ptw"` warping,
while `penalty` and `maxshift` apply only to `"vpdtw"` warping.

## Note

Adapted from the `correctRT` function in the alsace package by Ron
Wehrens
(<https://github.com/rwehrens/alsace/blob/master/R/correctRT.R>).

## References

- Clifford, D., Stone, G., Montoliu, I., Rezzi, S., Martin, F. P., Guy,
  P., Bruce, S., & Kochhar, S. 2009. Alignment using variable penalty
  dynamic time warping. *Analytical chemistry*, 81(3):1000-1007.
  [doi:10.1021/ac802041e](https://doi.org/10.1021/ac802041e) .

- Clifford, D., & Stone, G. 2012. Variable Penalty Dynamic Time Warping
  Code for Aligning Mass Spectrometry Chromatograms in R. *Journal of
  Statistical Software*, 47(8):1-17.
  [doi:10.18637/jss.v047.i08](https://doi.org/10.18637/jss.v047.i08) .

- Eilers, P.H.C. 2004. Parametric Time Warping. *Analytical Chemistry*,
  76:404-411. [doi:10.1021/ac034800e](https://doi.org/10.1021/ac034800e)
  .

- Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast parametric
  time warping of peak lists. *Bioinformatics*, 31:3063-3065.
  [doi:10.1093/bioinformatics/btv299](https://doi.org/10.1093/bioinformatics/btv299)
  .

- Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
  LC–DAD using multivariate curve resolution: the alsace package for R.
  *Metabolomics*, 11:143-154.
  [doi:10.1007/s11306-014-0683-5](https://doi.org/10.1007/s11306-014-0683-5)
  .

## See also

[`correct_rt_group`](https://ethanbass.github.io/chromatographR/reference/correct_rt_group.md),
[`correct_peaks`](https://ethanbass.github.io/chromatographR/reference/correct_peaks.md),
[`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html),
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
