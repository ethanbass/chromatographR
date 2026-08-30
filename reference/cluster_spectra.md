# Cluster spectra

Cluster peaks by spectral similarity.

## Usage

``` r
cluster_spectra(
  peak_table,
  alpha = 0.05,
  min_size = 5,
  max_size = NULL,
  nboot = 1000,
  plot_dend = TRUE,
  plot_spectra = TRUE,
  verbose = getOption("verbose"),
  outfile = NULL,
  save = NULL,
  parallel = TRUE,
  max.only = FALSE,
  output = c("pvclust", "clusters"),
  ...
)
```

## Arguments

- peak_table:

  A `peak_table` object from
  [`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md).

- alpha:

  Significance threshold used to select clusters based on bootstrap AU
  p-values.

- min_size:

  Minimum number of peaks a cluster may have.

- max_size:

  Maximum number of peaks a cluster may have.

- nboot:

  Number of bootstrap replicates for
  [`pvclust`](https://rdrr.io/pkg/pvclust/man/pvclust.html).

- plot_dend:

  Logical. If `TRUE`, plots dendrogram with bootstrap values.

- plot_spectra:

  Logical. If `TRUE`, plots overlapping spectra for each cluster.

- verbose:

  Logical. If `TRUE`, prints progress report to console.

- outfile:

  Optional path to an `.RDS` file. If not `NULL`, the `pvclust` object
  is saved file using
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html). Existing files are
  overwritten.

- save:

  Logical. If `TRUE`, saves pvclust object to current directory.

- parallel:

  Logical. If `TRUE`, use parallel processing for
  [`pvclust`](https://rdrr.io/pkg/pvclust/man/pvclust.html).

- max.only:

  Logical. If `TRUE`, returns only highest level cluster from each
  nested set of clusters.

- output:

  What to return. Either `"clusters"` to return list of clusters,
  `"pvclust"` to return pvclust object, or `"both"` to return both
  items.

- ...:

  Additional arguments to
  [`pvclust`](https://rdrr.io/pkg/pvclust/man/pvclust.html).

## Value

Returns clusters and/or `pvclust` object according to the value of the
`output` argument.

- If `output = clusters`, returns a list of S4 `cluster` objects.

- If `output = pvclust`, returns a `pvclust` object.

- If `output = both`, returns a nested list containing `[[1]]` the
  `pvclust` object, and `[[2]]` the list of S4 `cluster` objects.

The `cluster` objects consist of the following components:

- `peaks`: a character vector containing the names of all peaks
  contained in the given cluster.

- `pval`: a numeric scalar containing the AU bootstrap p-value for the
  given cluster.

## Details

Before using this function, reference spectra must be attached to the
`peak_table` using the
[`attach_ref_spectra`](https://ethanbass.github.io/chromatographR/reference/attach_ref_spectra.md)
function. These spectra are then used to construct a distance matrix
derived from pairwise Pearson correlations among peaks. Hierarchical
clustering with bootstrap resampling is performed on the resulting
correlation matrix to classify peaks by their spectral similarity as
implemented in
[`pvclust`](https://rdrr.io/pkg/pvclust/man/pvclust.html). Finally,
bootstrap values can be used to select clusters that exceed a certain
confidence specified by `alpha`.

Clusters can be filtered by the minimum and maximum size of the cluster
using the `min_size` and `max_size` arguments respectively. Users should
be aware that the clustering algorithm will often return nested
clusters. Thus, an individual peak could appear in more than one
cluster. If `max.only` is `TRUE`, only the largest cluster in a nested
tree of clusters meeting the specified confidence threshold will be
returned.

It is strongly recommended to use more than 100 bootstraps if you run
the clustering algorithm on real data even though we use `nboot = 100`
in the example to reduce run time. The authors of `pvclust` suggest
using at least 10,000 bootstrap replicates.

## References

R. Suzuki & H. Shimodaira. 2006. Pvclust: an R package for assessing the
uncertainty in hierarchical clustering. *Bioinformatics*,
22(12):1540-1542.
[doi:10.1093/bioinformatics/btl117](https://doi.org/10.1093/bioinformatics/btl117)
.

## Author

Ethan Bass

## Examples

``` r
# \donttest{
if (requireNamespace("pvclust", quietly = TRUE)) {
data(pk_tab)
data(Sa_warp)
pk_tab <- attach_ref_spectra(pk_tab, Sa_warp, ref = "max.int")
cl <- cluster_spectra(pk_tab, nboot = 100, max.only = FALSE, 
save = FALSE, alpha = 0.03)
}










# }
```
