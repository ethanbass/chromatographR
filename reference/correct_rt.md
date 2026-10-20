# Retention time alignment via time warping

Chromatograms are aligned using one of three algorithms, according to
the value of `alg`: parametric time warping (`"ptw"`), as implemented in
[`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html); variable penalty dynamic
time warping (`"vpdtw"`), as implemented in
[`VPdtw`](https://rdrr.io/pkg/VPdtw/man/VPdtw.html); or landmark
parametric time warping (`"lptw"`), which fits the warping polynomial
directly to matched peak positions rather than full chromatographic
signals.

## Usage

``` r
correct_rt(
  chrom_list,
  lambdas,
  models = NULL,
  reference = "best",
  alg = c("ptw", "lptw", "vpdtw"),
  peaks = NULL,
  what = c("corrected.values", "models"),
  init.coef = c(0, 1, 0),
  n_traces = NULL,
  fill_zeros = FALSE,
  n_zeros = 0,
  scale = FALSE,
  trwdth = 0.2,
  plot_it = FALSE,
  penalty = 0.2,
  max_shift = 1,
  dilation_span = 0.3,
  rt_window = 0.3,
  spectral_cutoff = 0.95,
  regularization = 0,
  verbose = getOption("verbose"),
  show_progress = NULL,
  cl = 2,
  ...,
  n.traces = NULL,
  n.zeros = NULL,
  maxshift = NULL
)
```

## Arguments

- chrom_list:

  List of chromatograms in matrix format, sharing a common time axis
  (see
  [`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md)).

- lambdas:

  A character or numeric vector specifying the wavelengths to use for
  alignment. `"vpdtw"` and `"lptw"` accept only a single wavelength and
  raise an error otherwise. For one-dimensional chromatograms, this
  argument can be omitted.

- models:

  List of models to warp by. The models provided here (if any) must
  match the algorithm selected in `alg`, and are matched to the samples
  of `chrom_list` by name.

- reference:

  Index or name of the sample to be used as the reference. Fitted models
  record it by name. If not specified, chosen by
  [`bestref`](https://rdrr.io/pkg/ptw/man/bestref.html) from a
  similarity matrix of the chromatograms.

- alg:

  Alignment algorithm to use: parametric time warping (`"ptw"`),
  landmark parametric time warping (`"lptw"`), or variable penalty
  dynamic time warping (`"vpdtw"`).

- peaks:

  A peak list object from
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)
  or a `matched_peak_list` object from
  [`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md).
  If a `matched_peak_list` is provided, peak matching will be skipped.
  Only used when `alg = "lptw"`. If `NULL`, peaks are found with
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)
  and matched with
  [`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md),
  both at their default settings. To change either, pass their output
  instead.

- what:

  Output type: either the `"corrected.values"` (useful for visual
  inspection and downstream analysis) or the warping `"models"` (for
  further programmatic use).

- init.coef:

  Starting values for the optimization. Named to match the argument of
  the same name in [`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html), which
  receives it. For `alg = "lptw"`, only its length is used, setting the
  degree of the warping polynomial.

- n_traces:

  Number of wavelengths to warp by, chosen from `lambdas` (or from all
  wavelengths if `lambdas` is missing) by the CODA criterion in
  [`select.traces`](https://rdrr.io/pkg/ptw/man/select.traces.html).
  Defaults to `NULL`, which uses every wavelength in `lambdas`.

- fill_zeros:

  Logical. If `TRUE`, out-of-bounds regions produced by warping are
  filled with zeros; otherwise they are returned as `NA` (default).

- n_zeros:

  Number of zeros to pad onto each end of the chromatograms before
  fitting. Only used when `alg = "ptw"`.

- scale:

  Logical. If `TRUE`, divide each chromatogram by its largest absolute
  value before fitting the warp, so that samples of very different
  concentration are matched on the shape of their signal. The corrected
  values are returned in their original units; to rescale the data
  themselves, use the `maxI` argument of
  [`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md).
  `"ptw"` is unaffected under its default criterion, weighted
  cross-correlation, which does not change when a sample is multiplied
  by a constant.

- trwdth:

  Width of the triangle in the weighted cross-correlation (WCC)
  criterion, passed to [`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html).
  Defaults to `0.2`, which matches the default of `20` timepoints in
  [`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html) itself on a chromatogram
  sampled every `0.01` minutes. Roughly a peak width is a reasonable
  choice.

- plot_it:

  Logical. Whether to plot alignment. Defaults to `FALSE`.

- penalty:

  Numeric \>= 0. Strength of the penalty against warping for
  `alg = "vpdtw"`: the constant of proportionality in \\\lambda
  \propto\\
  [`dilation`](https://rdrr.io/pkg/VPdtw/man/dilation.html)`(reference)`
  over a `dilation_span`-wide window, as recommended by Clifford and
  Stone (2012). A higher number allows less warping; `0` applies no
  penalty at all. Defaults to `0.2`. In chromatographR \<= 0.8.0
  `penalty` was a divisor, so use the reciprocal of an old value. Values
  above `1` raise a warning, since they are more likely to be old
  divisors than deliberate choices.

- max_shift:

  Maximum allowable shift for `VPdtw` warping, in minutes. Set it from
  the retention time drift you actually expect; the default of `1`
  allows a minute of drift. Although `max_shift` is a time, `VPdtw`
  itself works on the index axis, so the `shift` element of the returned
  models is reported in timepoints.

- dilation_span:

  Width of the moving-maximum window used to build the `"vpdtw"`
  penalty. The penalty at each point is the largest reference intensity
  within `dilation_span` either side of it, computed by
  [`dilation`](https://rdrr.io/pkg/VPdtw/man/dilation.html) after
  conversion to timepoints. Smaller values track individual peaks more
  closely, allowing more warping in the baseline between them. Larger
  values approach a constant penalty, which becomes exactly constant
  once `dilation_span` exceeds the distance between the tallest peak of
  the reference and the furthest end of the chromatogram. Defaults to
  `0.3`, around a couple of peak widths on the example data.

- rt_window:

  Numeric. Argument to
  [`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md).
  Maximum allowed difference in retention time, in minutes, between a
  reference and sample peak to be considered a candidate match. Defaults
  to `0.3`.

- spectral_cutoff:

  Numeric. Argument to
  [`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md).
  Minimum spectral correlation for a candidate peak pair to be retained.
  Defaults to `0.95`.

- regularization:

  Numeric \>= 0. Strength of the penalty against warping for
  `alg = "lptw"`: a dimensionless weight pulling the warp toward the
  identity warp, passed as `penalty` to `fit_penalized_warp`. A higher
  number allows less warping; `0` (the default) applies no penalty at
  all. `1` weights a unit of squared warp displacement equally against a
  unit of squared peak-position mismatch, and around `2.4` halves the
  applied warp, so useful values run roughly `0`-`10`. Unbounded above,
  approaching the identity warp asymptotically.

- verbose:

  Logical. Whether to print a progress message as the warping models are
  fit and applied. Defaults to `getOption("verbose")`.

- show_progress:

  Logical. Whether to show progress bar. Defaults to `TRUE` if
  [`pbapply`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  is installed. Has no effect when `alg = "vpdtw"`.

- cl:

  Either an integer specifying the number of cores to use for parallel
  processing or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to `2`. On Windows systems, integers are ignored. Has no effect when
  `alg = "vpdtw"`.

- ...:

  For `alg = "ptw"`, optional additional arguments passed to
  [`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html). The only argument that
  cannot be changed is `warp.type`, which is hard-coded to `"global"` to
  permit warping on multiple wavelengths. For `alg = "lptw"`, optional
  arguments passed to `fit_penalized_warp`: `penalty_scale`
  (per-coefficient penalty multipliers) and `min_peaks` (minimum matched
  peaks required to fit the polynomial). Its `penalty` argument is set
  by `regularization` rather than through `...`.

- n.traces:

  The `n.traces` argument is deprecated. Please use `n_traces` instead.

- n.zeros:

  The `n.zeros` argument is deprecated. Please use `n_zeros` instead.

- maxshift:

  The `maxshift` argument is deprecated. Please use `max_shift`, which
  is in minutes. `maxshift` is read as a count of timepoints, as in
  chromatographR 0.8.0, and converted.

## Value

A list of warping models or a list of warped absorbance profiles,
according to the value of the `what` argument. With `alg = "vpdtw"`, the
warped chromatograms share a new time axis, which can be longer than the
input.

## Details

Some arguments are specific to particular warping functions. For example
`n_traces` applies only to `"ptw"` warping; `penalty`, `max_shift`, and
`dilation_span` apply only to `"vpdtw"` warping; and `peaks`,
`rt_window`, `spectral_cutoff`, and `regularization` apply only to
`"lptw"` warping.

`trwdth`, `max_shift` and `dilation_span` are widths along the time
axis, given in minutes, regardless of the time unit of `chrom_list`, and
rounded internally to a whole number of timepoints. `trwdth` was a count
of timepoints up to chromatographR 0.8.0; multiply an old value by the
sampling interval in minutes to convert it. A width wider than the
chromatogram itself is rejected with an error, since it is almost
certainly a stale timepoint count. `max_shift` replaces `maxshift`,
which is still accepted as a count of timepoints.

`regularization` (`"lptw"`) and `penalty` (`"vpdtw"`) are not
comparable: both are dimensionless and higher means less warping, but
they are on unrelated scales.

## Note

The PTW branch of `correct_rt` is adapted from the `correctRT` function
in Ron Wehrens' [alsace](https://github.com/rwehrens/alsace) package.

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
data(Sa_pr)
warp <- correct_rt(chrom_list = Sa_pr, lambdas=210)
```
