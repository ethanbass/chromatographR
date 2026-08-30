# Correct retention time with group-based warping

Aligns chromatograms using parametric time warping (`"ptw"`), as
implemented in [`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html), with
warping functions estimated from within-group averages. This is useful
when samples fall into batches with shared retention time shifts, or
when individual samples contain peaks absent in other groups that would
otherwise confound alignment. In addition to potentially being more
accurate, this should also be faster than computing warping functions
individually on each sample.

## Usage

``` r
correct_rt_group(
  chrom_list,
  lambdas,
  groups,
  reference = "best",
  reference_group = NULL,
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

  A list of chromatograms in matrix format.

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

  Index or name of the group average to use as the alignment reference.
  Defaults to `"best"`, which selects the reference automatically using
  [`bestref`](https://rdrr.io/pkg/ptw/man/bestref.html) applied to the
  group averages.

- reference_group:

  Name of the group to use as the reference group. If supplied,
  overrides `reference`. Defaults to `NULL`.

- init.coef:

  Starting values for the optimization.

- n.traces:

  Number of traces to use.

- fill_zeros:

  Logical. If `TRUE`, out-of-bounds regions produced by warping are
  filled with zeros. If `FALSE` (default), these regions are returned as
  `NA`.

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

  Argument to `pbapply` or
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html). Either an
  integer specifying the number of clusters to use for parallel
  processing or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). Defaults
  to `2`. On Windows systems, integer values will be ignored.

- ...:

  Optional additional arguments to `ptw`. The only argument that cannot
  be changed is `warp.type` which is hard-coded to `"global"` to permit
  warping on multiple wavelengths.

## Value

A list of warped chromatogram matrices in the same order as
`chrom_list`, with each sample warped using the warp coefficients
estimated from its group average.

## See also

[`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md),
[`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html)

## Author

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
if (FALSE) { # \dontrun{
data(Sa_pr)
groups <- c("a", "a", "b", "b")
warp <- correct_rt_group(chrom_list = Sa_pr, lambdas = 210, groups = groups)
} # }
}
```
