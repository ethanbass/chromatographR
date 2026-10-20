# Convert peak lists into an aligned peak table

Converts a `peak_list` object containing peaks from multiple samples
into an aligned peak table by clustering peaks across samples according
to retention time and (optionally) spectral similarity. The tree is cut
at height `hmax`, whose meaning depends on `clust` (see `hmax`).

## Usage

``` r
get_peaktable(peak_list, ...)

# Default S3 method
get_peaktable(
  peak_list,
  chrom_list,
  response = c("area", "height"),
  use_corrected_rt = NULL,
  hmax = 0.2,
  summarize_by = c("weighted.mean", "median", "mean", "max"),
  plot_it = FALSE,
  ask = plot_it,
  clust = c("rt", "sp.rt"),
  sigma_t = NULL,
  sigma_r = 0.5,
  deepSplit = FALSE,
  verbose = FALSE,
  ...,
  sigma.t = NULL,
  sigma.r = NULL,
  use.cor = NULL
)

# S3 method for class 'matched_peak_list'
get_peaktable(
  peak_list,
  pks,
  response = c("area", "height"),
  rt_tol = 0.01,
  na.fill = 0,
  verbose = FALSE,
  ...
)

# S3 method for class 'matched_peak_list_set'
get_peaktable(
  peak_list,
  pks,
  response = c("area", "height"),
  rt_tol = 0.01,
  na.fill = 0,
  verbose = FALSE,
  ...
)
```

## Arguments

- peak_list:

  A `peak_list` object created by
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md),
  or a `matched_peak_list` created by
  [`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md).

- ...:

  Additional arguments (currently unused).

- chrom_list:

  A list of chromatographic matrices. Only required when
  `clust = "sp.rt"`.

- response:

  Peak intensity measure: either `"area"` (default) or `"height"`.

- use_corrected_rt:

  Logical. Indicates whether to use corrected retention times (`rt.cor`
  column) or raw retention times (`rt` column). Unless otherwise
  specified, the `rt.cor` column will be used by default if it exists in
  the provided `peak_list`. Either way, they are stored in the `rt`,
  `start` and `end` rows of `pk_meta`, and `args$use_corrected_rt`
  records which were used.

- hmax:

  Height at which the clustering tree is cut. With `clust = "rt"`, peaks
  are clustered by complete linkage on retention time, so `hmax` is the
  largest retention time difference allowed between any two peaks in the
  same cluster, in minutes. With `clust = "sp.rt"`, peaks are clustered
  by average linkage on the dissimilarity \\1 - S\_{ij}\\, which runs
  from 0 to 1, and `hmax` is the maximum tree height passed to
  [`cutreeDynamicTree`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamicTree.html);
  clusters may be split below it (see `deepSplit`). There, the default
  of `0.2` lets peaks with identical spectra sit about `0.67 * sigma_t`
  apart, a third of the mean `end - start` at the default `sigma_t`.
  Ignored for peak lists produced by `get_peaks(windows = )` or
  `get_peaks(targets = )`, since those peaks are grouped by integration
  window or target instead of being clustered (see *Windows and
  targets*).

- summarize_by:

  How to summarize peak metadata within each cluster. Options are
  `"weighted.mean"` (default, weights peaks by their intensity as given
  by `response`), `"mean"`, `"median"` (which aggregate metadata across
  all peaks in the cluster) or `"max"` (which selects the most intense
  peak in the cluster and uses its metadata directly).

- plot_it:

  Logical. If `TRUE`, for every component a strip plot will be shown
  indicating the clustering.

- ask:

  Logical. Ask before showing new plot? Defaults to `plot_it`.

- clust:

  Specify whether to perform hierarchical clustering based on spectral
  similarity and retention time (`sp.rt`) or retention time alone
  (`rt`). Defaults to `rt`. Ignored for peak lists produced by
  `get_peaks(windows = )` or `get_peaks(targets = )`, which are grouped
  by window or target (see *Windows and targets*).

- sigma_t:

  Width of Gaussian kernel for retention time similarity, in minutes.
  Controls the weight given to retention time when `clust = "sp.rt"`.
  Defaults to half the mean `end - start` across all peaks.

- sigma_r:

  Width of Gaussian kernel for spectral similarity. Controls the weight
  given to spectral correlation when `clust = "sp.rt"`. Defaults to
  `0.5`.

- deepSplit:

  Logical. Controls sensitivity to cluster splitting when
  `clust = "sp.rt"`. If `TRUE`, returns a larger number of smaller
  clusters. Defaults to `FALSE`. For additional information, see
  [`cutreeDynamicTree`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamicTree.html).

- verbose:

  Logical. Whether to warn when two peaks from the same sample fall in
  the same cluster. Defaults to `FALSE`.

- sigma.t:

  The `sigma.t` argument is deprecated. Please use `sigma_t`, which is
  in minutes. `sigma.t` is read in the units of the retention times, as
  before.

- sigma.r:

  The `sigma.r` argument is deprecated. Please use `sigma_r` instead.

- use.cor:

  The `use.cor` argument is deprecated. Please use `use_corrected_rt`
  instead.

- pks:

  A `peak_list` object from
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md).
  Required when `peak_list` is a `matched_peak_list`. Used to retrieve
  per-peak intensities and metadata.

- rt_tol:

  Numeric. Maximum allowed deviation (in retention-time units) when
  joining matched peak RTs back to the corresponding row in `pks`.
  Defaults to `0.01`.

- na.fill:

  Value used to fill cells where no peak was matched in a given sample.
  Defaults to `0`. Only used by the `matched_peak_list` method.

## Value

An S3
[`peak_table`](https://ethanbass.github.io/chromatographR/reference/peak_table-class.md)
object, containing the following elements:

- `tab`: The peak table itself — a `data.frame` of intensities in a
  sample × peak configuration. A peak that was not found in a sample, or
  whose response is negative, is recorded as `0` (`na.fill` for a
  `matched_peak_list`). A peak whose model failed, which
  `get_peaks(targets = )` keeps with a missing `height` and `area`, is
  recorded as `NA`.

- `pk_meta`: A `data.frame` of peak metadata, one column per peak. Its
  rows are the wavelength (`lambda`), the peak number (`peak`), and
  whichever of the peak list's `rt`, `start`, `end`, `sd`, `tau`,
  `FWHM`, `r.squared`, `purity`, `window_limited` and `converged`
  columns it has. See the `Value` section of
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)
  for their meanings. The rows taken from the peak list are summarized
  across samples according to `summarize_by`; by default, peaks are
  weighted by intensity, which reduces the influence of small peaks.
  With `"weighted.mean"` or `"mean"`, the `converged` and
  `window_limited` rows become the share of the column's peaks with that
  flag, so `0.25` means a quarter of them. In a table built from a
  `matched_peak_list`, the rows hold the reference sample's values
  instead.

- `sample_meta`: A `data.frame` of sample metadata: the fields added to
  the chromatograms with
  [`chromConverter::add_metadata()`](https://ethanbass.github.io/chromConverter/reference/add_metadata.html)
  (chromConverter 0.11.0 or later), and any attached with
  [`attach_metadata`](https://ethanbass.github.io/chromatographR/reference/attach_metadata.md).
  `NA` if there are none.

- `instrument_meta`: A `data.frame` of instrumental metadata transferred
  automatically from the chromatograms.

- `ref_spectra`: A matrix of reference spectra (in a wavelength × peak
  configuration). `NA` until added with
  [`attach_ref_spectra`](https://ethanbass.github.io/chromatographR/reference/attach_ref_spectra.md).

- `args`: A list of the arguments given to `get_peaktable` to generate
  the peak table, with the `fit`, `baseline` and `max_iter` the peaks
  were fitted with in
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md).
  Its `filters` element holds any
  [`filter_peaks`](https://ethanbass.github.io/chromatographR/reference/filter_peaks.md)
  and
  [`filter_peaktable`](https://ethanbass.github.io/chromatographR/reference/filter_peaktable.md)
  calls applied, as code, separated by `"; "`.

## Note

This function is adapted from the `getPeakTable` function in the alsace
package by Ron Wehrens:
<https://github.com/rwehrens/alsace/blob/master/R/getPeakTable.R>.

## Clustering

By default, clustering is performed on retention times alone (when
`clust = "rt"`). Alternatively, clustering can also incorporate
information about spectral similarity (when `clust = "sp.rt"`) using a
similarity function adapted from Broeckling et al., 2014: \$\$ S\_{ij} =
\exp\left(-\frac{(t_i - t_j)^2}{2\sigma_t^2}\right) \cdot
\exp\left(-\frac{(1 - c\_{ij})^2}{2\sigma_r^2}\right) \$\$ where
\\c\_{ij}\\ is the spectral correlation between spectra \\i\\ and \\j\\,
and \\\sigma_t\\ and \\\sigma_r\\ control the relative contributions of
retention time and spectral similarity, respectively (see `sigma_t` and
`sigma_r`).

When two peaks from the same sample fall in the same cluster, only the
more intense one enters the peak table, with a warning if
`verbose = TRUE`. Such warnings are often benign, but they can also mean
`hmax` is too large. Lowering `hmax` merges fewer peaks, at the cost of
splitting some across multiple clusters; filtering low-intensity peaks
before clustering reduces the warnings without that trade-off.

## Windows and targets

Peak lists produced by `get_peaks(windows = )` are grouped by
integration window rather than clustered by retention time, including
after
[`correct_peaks`](https://ethanbass.github.io/chromatographR/reference/correct_peaks.md).
This gives one peak table column per window, ordered by retention time.

Peak lists produced by `get_peaks(targets = )` are grouped by target
instead of being clustered, each peak going to the matching target.
Targets found in no sample are omitted from the peak table. If `targets`
has a `name` column, the columns are named after the targets.

Peak tables built from `get_peaks(windows = )` or
`get_peaks(targets = )` keep those `windows` or `targets` as an
attribute of the same name, as a record of how the table was built. The
`targets` attribute lists every target that was searched for, including
any found in no sample, which have no column. Neither attribute is
updated when columns are removed.

## Matched peak lists

[`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md)
pairs the peaks in each sample with the peaks of a reference sample, by
spectrum as well as retention time. Each column then follows one
compound across the samples, even where retention times drift too far
for clustering to group the peaks reliably. In crowded stretches of the
chromatogram, this works only if neighboring compounds have distinct
spectra. This suits targeted analysis (see
[`targeted_analysis`](https://ethanbass.github.io/chromatographR/reference/targeted_analysis.md))
and studies in which one sample, such as a standard mix or a pooled
sample, contains every compound of interest.

The peak table has one column for each peak in the reference sample, at
each of the specified wavelengths. A peak that's missing in the
reference sample gets no column, even if it's present in other samples,
so choose the reference with care (see
[`match_peaks`](https://ethanbass.github.io/chromatographR/reference/match_peaks.md)).

## References

- Broeckling, C. D., Afsar F.A., Neumann S., Ben-Hur A., and Prenni
  J.E. 2014. RAMClust: A Novel Feature Clustering Method Enables
  Spectral-Matching-Based Annotation for Metabolomics Data. *Anal.
  Chem.* 86:6812-6817.
  [doi:10.1021/ac501530d](https://doi.org/10.1021/ac501530d) .

- Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
  LC–DAD using multivariate curve resolution: the alsace package for R.
  *Metabolomics* 11:143-154.
  [doi:10.1007/s11306-014-0683-5](https://doi.org/10.1007/s11306-014-0683-5)
  .

## See also

[`peak_table`](https://ethanbass.github.io/chromatographR/reference/peak_table-class.md),
[`attach_ref_spectra`](https://ethanbass.github.io/chromatographR/reference/attach_ref_spectra.md),
[`attach_metadata`](https://ethanbass.github.io/chromatographR/reference/attach_metadata.md)

## Author

Ron Wehrens, Ethan Bass

## Examples

``` r
data(Sa_pr)
pks <- get_peaks(Sa_pr, lambdas = c('210'))
get_peaktable(pks, response = "area")
#>           V1       V2        V3         V4         V5        V6         V7
#> 119 0.000000 5.160717 0.6216202 16.6408509  0.4577718 0.4599784 81.1954044
#> 121 0.000000 4.012988 0.8692032  0.5517739  7.7387426 0.0000000  0.0000000
#> 122 0.000000 8.168948 1.3471271  2.0773758 25.8641682 2.2710192 59.7697574
#> 458 5.640076 5.769393 1.7765152  1.7162467 18.6509054 0.0000000  0.2636266
#>             V8        V9      V10      V11      V12       V13      V14
#> 119   2.725664 42.309860  1.55284 1.537438 2.509823 37.504493 14.03672
#> 121  57.367316  0.000000 18.42121 0.000000 1.022751  8.082183 20.06353
#> 122   1.164881  0.000000 36.05080 1.536124 4.585314 32.383822 12.62787
#> 458 110.608512  3.767623 33.80274 4.325797 1.389959  5.667755 43.02514
#>           V15      V16       V17      V18      V19      V20        V21      V22
#> 119  0.000000 34.13395  0.000000 4.567092 0.000000 5.128388 1.95140616 12.59202
#> 121  8.359219  0.00000 22.311973 0.000000 1.171736 1.812437 0.09512587  0.00000
#> 122  0.000000 28.69373  1.215313 2.452972 0.000000 3.390005 0.59202006  0.00000
#> 458 23.415509  0.00000 25.148310 0.000000 5.164612 5.576050 0.34027136  0.00000
#>           V23      V24       V25         V26        V27      V28       V29
#> 119 14.918412  0.00000 2.3572466 0.092200823 17.6381605  0.00000 1.7600985
#> 121  6.833451 11.80269 0.7327165 0.007687478  0.3798519 13.93831 0.5340067
#> 122  8.923256 15.44319 2.7602610 0.080516725  0.0000000 16.53251 0.8198640
#> 458 12.455040 14.88247 2.2405245 0.004549712  0.0000000 23.49422 0.9877410
#>           V30       V31       V32      V33       V34      V35      V36
#> 119 0.0000000 1.2108244 0.9140629 1.683199 0.8260853 1.447895 1.980108
#> 121 0.2088148 0.4814268 0.4277401 0.000000 2.0490665 0.000000 1.746997
#> 122 0.0000000 1.0992387 1.1170352 0.000000 0.0000000 0.000000 5.719170
#> 458 0.5647096 0.9305848 1.3541306 1.668973 2.9601134 0.000000 3.644116
```
