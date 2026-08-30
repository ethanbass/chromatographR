# Convert peak lists into an aligned peak table.

Converts a `peak_list` object containing peaks from multiple samples
into an aligned peak table by clustering peaks across samples according
to retention time and (optionally) spectral similarity. Clusters are
defined by cutting a complete-linkage dendrogram at height `hmax`, which
can be understood as the maximal allowed clustering distance.

## Usage

``` r
get_peaktable(
  peak_list,
  chrom_list,
  response = c("area", "height"),
  use.cor = NULL,
  hmax = 0.2,
  summarize_by = c("weighted.mean", "median", "mean", "max"),
  plot_it = FALSE,
  ask = plot_it,
  clust = c("rt", "sp.rt"),
  sigma.t = NULL,
  sigma.r = 0.5,
  deepSplit = FALSE,
  verbose = FALSE,
  out = c("data.frame", "matrix")
)
```

## Arguments

- peak_list:

  A `peak_list` object created by
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md),
  containing a nested list of features, where the first level is the
  sample, and the second level is the spectral wavelength. Every sample
  x wavelength component is described by a `data.frame` with a row for
  each peak and columns containing information on various peak
  parameters.

- chrom_list:

  A list of chromatographic matrices. Only required when
  `clust = "sp.rt"`.

- response:

  Peak intensity measure: either `"area"` (default) or `"height"`.

- use.cor:

  Logical. Indicates whether to use corrected retention times (`rt.cor`
  column) or raw retention times `rt` column). Unless otherwise
  specified, the `rt.cor` column will be used by default if it exists in
  the provided `peak_list`.

- hmax:

  Height at which the complete linkage dendrogram will be cut. Can be
  interpreted as the maximal intercluster retention time difference.

- summarize_by:

  How to select the representative peak from each cluster. Options are
  `"weighted.mean"` (default, weights peaks by their intensity as given
  by `response`), `"mean"`, `"median"` (which aggregate metadata across
  all peaks in the cluster) or `"max"` (which selects the most intense
  peak in the cluster and uses its metadata directly).

- plot_it:

  Logical. If `TRUE`, for every component a strip plot will be shown
  indicating the clustering.

- ask:

  Logical. Ask before showing new plot? Defaults to `TRUE`.

- clust:

  Specify whether to perform hierarchical clustering based on spectral
  similarity and retention time (`sp.rt`) or retention time alone
  (`rt`). Defaults to `rt`.

- sigma.t:

  Width of Gaussian kernel for retention time similarity. Controls
  weight given to retention time when `"sp.rt"` is selected.

- sigma.r:

  Width of Gaussian kernel for spectral similarity. Controls weight
  given to spectral correlation if `"sp.rt"` is selected.

- deepSplit:

  Logical. Controls sensitivity to cluster splitting. If `TRUE`,
  function will return a larger number of smaller clusters. For
  additional information, see
  [`cutreeDynamic`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html).

- verbose:

  Logical. Whether to print warning when combining peaks into single
  time window. Defaults to `FALSE`.

- out:

  Specify `data.frame` (default) or `matrix` as output.

## Value

The function returns an S3
[`peak_table`](https://ethanbass.github.io/chromatographR/reference/peak_table-class.md)
object, containing the following elements:

- `tab`: The peak table itself – a `data.frame` of intensities in a
  sample × peak configuration.

- `pk_meta`: A `data.frame` containing peak metadata (e.g., the spectral
  component, peak number, and average retention time), summarized across
  all peaks in each cluster according to `summarize_by`.

- `sample_meta`: A `data.frame` of sample metadata. Must be added using
  [`attach_metadata`](https://ethanbass.github.io/chromatographR/reference/attach_metadata.md).

- `instrument_meta`: A `data.frame` of instrumental metadata transferred
  automatically from the chromatograms.

- `ref_spectra`: A `data.frame` of reference spectra (in a wavelength ×
  peak configuration). Must be added using
  [`attach_ref_spectra`](https://ethanbass.github.io/chromatographR/reference/attach_ref_spectra.md).

- `args`: A vector of arguments given to `get_peaktable` to generate the
  peak table.

## Details

By default, clustering is performed on retention times alone (when
`clust = "rt"`). Alternatively, clustering can also incorporate
information about spectral similarity (when `clust = "sp.rt"`) using a
similarity function adapted from Broeckling et al., 2014: \$\$ S\_{ij} =
\exp\left(-\frac{(t_i - t_j)^2}{2\sigma_t^2}\right) \cdot
\exp\left(-\frac{(1 - c\_{ij})^2}{2\sigma_r^2}\right) \$\$ where
\\c\_{ij}\\ is the spectral correlation between spectra \\i\\ and \\j\\,
and \\\sigma_t\\ and \\\sigma_r\\ control the relative contributions of
retention time and spectral similarity, respectively (see `sigma.t` and
`sigma.r`).

If two peaks from the same sample are assigned to the same cluster, a
warning message is printed to the console. These warnings are often
benign but may also indicate that `hmax` is too large. Reducing `hmax`
may prevent excessive clustering, although overly small values may
instead lead to splitting of peaks across multiple clusters. Filtering
low-intensity peaks prior to clustering may also help reduce these
warnings.

Once clusters are formed, peak metadata are summarized within each
cluster according to `summarize_by`. By default (`"weighted.mean"`),
peak metadata are averaged using peak intensity as weights, reducing the
influence of small peaks which are more likely to represent noise.
Alternatively, the `"max"` argument selects the metadata associated with
the most intense peak in each cluster.

## Note

This function is adapted from the `getPeakTable` function in the alsace
package by Ron Wehrens:
<https://github.com/rwehrens/alsace/blob/master/R/getPeakTable.R>.

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

Ethan Bass

## Examples

``` r
if (FALSE) { # interactive()
data(Sa_pr)
pks <- get_peaks(Sa_pr, lambdas = c('210'))
get_peaktable(pks, response = "area")
}
```
