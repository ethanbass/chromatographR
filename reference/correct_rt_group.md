# Correct retention time with group-based warping

Aligns chromatograms using parametric time warping (`"ptw"`), as
implemented in [`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html), or
landmark parametric time warping (`"lptw"`), with warping functions
estimated from within-group averages. This is useful when samples fall
into batches with shared retention time shifts, or when individual
samples contain peaks absent in other groups that would otherwise
confound alignment. One warp is fit per group rather than one per
sample, so the cost of warping falls roughly in proportion to the number
of groups.

## Usage

``` r
correct_rt_group(
  chrom_list,
  lambdas,
  groups,
  reference = "best",
  reference_group = NULL,
  alg = c("ptw", "lptw"),
  init.coef = c(0, 1, 0),
  n_traces = NULL,
  fill_zeros = FALSE,
  n_zeros = 0,
  scale = FALSE,
  trwdth = 0.2,
  plot_it = FALSE,
  penalty = NULL,
  max_shift = 1,
  rt_window = 0.3,
  spectral_cutoff = 0.95,
  regularization = 0,
  verbose = getOption("verbose"),
  show_progress = NULL,
  cl = 2,
  what = c("corrected.values", "models"),
  ...,
  n.traces = NULL,
  n.zeros = NULL,
  maxshift = NULL
)
```

## Arguments

- chrom_list:

  A list of chromatograms in matrix format, sharing a common time axis
  (see
  [`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md)).

- lambdas:

  A character or numeric vector specifying the wavelengths to use for
  alignment.

- groups:

  A vector of group assignments for each chromatogram, or a single
  string naming a metadata attribute from which to extract the group
  assignments (e.g. `"batch"`). If a vector is provided, it must be the
  same length as `chrom_list` and may be named (matching
  `names(chrom_list)`) or positional. All samples must have a group
  assignment; `NA` values will trigger an error.

- reference:

  Index or name of a sample in `chrom_list` whose group will be used as
  the alignment reference. Defaults to `"best"`, which selects the
  reference group automatically using
  [`bestref`](https://rdrr.io/pkg/ptw/man/bestref.html) applied to the
  group averages.

- reference_group:

  Name of the group to use as the reference group. If supplied,
  overrides `reference`. Defaults to `NULL`.

- alg:

  Alignment algorithm to use: parametric time warping (`"ptw"`, the
  default) or landmark parametric time warping (`"lptw"`).

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

  Unused. `correct_rt_group` does not support `"vpdtw"`, so the
  `"vpdtw"` penalty has no effect here.

- max_shift:

  Unused. `correct_rt_group` does not support `"vpdtw"`, so the
  `"vpdtw"` shift limit has no effect here.

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

- what:

  What to return: `"corrected.values"` (the default) for the warped
  chromatograms, or `"models"` for a `ptw_list` or `lptw_list` holding
  each sample's group warp, which can be passed to
  [`correct_peaks`](https://ethanbass.github.io/chromatographR/reference/correct_peaks.md)
  or to
  [`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md)
  via `models`.

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

  Deprecated and unused. See `max_shift`.

## Value

If `what = "corrected.values"`, a list of warped chromatogram matrices
in the same order as `chrom_list`, with each sample warped using the
warp coefficients estimated from its group average. If
`what = "models"`, a `ptw_list` (`alg = "ptw"`) or `lptw_list`
(`alg = "lptw"`) with one model per sample, in the same order.

## Details

With `alg = "lptw"`, peaks are found in each group average by
[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)
at its default settings and matched to the peaks of the reference
group's average, so landmarks come from the averaged signal rather than
from individual samples.

## See also

[`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md),
[`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html)

## Author

Ethan Bass

## Examples

``` r
data(Sa_pr)
groups <- c("a", "a", "b", "b")
warp <- correct_rt_group(chrom_list = Sa_pr, lambdas = 210, groups = groups)
```
