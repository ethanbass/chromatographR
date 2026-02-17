# Converts peak list into an ordered peak table.

The function performs a complete linkage clustering of retention times
across all samples, and cuts at a height given by the user (which can be
understood as the maximal inter-cluster retention time difference) in
the simple case based on retention times. Clustering can also
incorporate information about spectral similarity using a distance
function adapted from Broeckling et al., 2014:
\$\$e^{-\frac{(1-c\_{ij})^2}{2\sigma_r^2}} \cdot
e^{-\frac{(1-(t_i-t_j)^2)}{2\sigma_t^2}}\$\$ If two peaks from the same
sample are assigned to the same cluster, a warning message is printed to
the console. These warnings can usually be ignored, but one could also
consider reducing the `hmax` variable. However, this may lead to
splitting of peaks across multiple clusters. Another option is to filter
the peaks by intensity to remove small features.

## Usage

``` r
get_peaktable(
  peak_list,
  chrom_list,
  response = c("area", "height"),
  use.cor = NULL,
  hmax = 0.2,
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
  [`get_peaks`](https://ethanbass.github.io/chromatographR/dev/reference/get_peaks.md),
  containing a nested list of peak tables: the first level is the
  sample, and the second level is the spectral wavelength. Every
  component is described by a `data.frame` with a row for each peak and
  columns containing information on various peak parameters.

- chrom_list:

  A list of chromatographic matrices.

- response:

  Indicates whether peak area or peak height is to be used as intensity
  measure. Defaults to `area` setting.

- use.cor:

  Logical. Indicates whether to use corrected retention times (`rt.cor`
  column) or raw retention times (`rt` column). Unless otherwise
  specified, the `rt.cor` column will be used by default if it exists in
  the provided `peak_list`.

- hmax:

  Height at which the complete linkage dendrogram will be cut. Can be
  interpreted as the maximal intercluster retention time difference.

- plot_it:

  Logical. If `TRUE`, for every component a strip plot will be shown
  indicating the clustering.

- ask:

  Logical. Ask before showing new plot? Defaults to `TRUE`.

- clust:

  Specify whether to perform hierarchical clustering based on spectral
  similarity and retention time (`sp.rt`) or retention time alone
  (`rt`). Defaults to `rt`. The `sp.rt` option is experimental and
  should be used with caution.

- sigma.t:

  Width of gaussian in retention time distance function. Controls weight
  given to retention time if `sp.rt` is selected.

- sigma.r:

  Width of gaussian in spectral similarity function. Controls weight
  given to spectral correlation if `sp.rt` is selected.

- deepSplit:

  Logical. Controls sensitivity to cluster splitting. If `TRUE`,
  function will return more smaller clusters. See documentation for
  [`cutreeDynamic`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html)
  for additional information.

- verbose:

  Logical. Whether to print warning when combining peaks into single
  time window. Defaults to `FALSE`.

- out:

  Specify `data.frame` or `matrix` as output. Defaults to `data.frame`.

## Value

The function returns an S3
[`peak_table`](https://ethanbass.github.io/chromatographR/dev/reference/peak_table-class.md)
object, containing the following elements:

- `tab`: The peak table itself – a `data.frame` of intensities in a
  sample x peak configuration.

- `pk_meta`: A `data.frame` containing peak meta-data (e.g., the
  spectral component, peak number, and average retention time).

- `sample_meta`: A `data.frame` of sample meta-data. Must be added using
  [`attach_metadata`](https://ethanbass.github.io/chromatographR/dev/reference/attach_metadata.md).

- `ref_spectra`: A `data.frame` of reference spectra (in a wavelength ×
  peak configuration). Must be added using
  [`attach_ref_spectra`](https://ethanbass.github.io/chromatographR/dev/reference/attach_ref_spectra.md).

- `args`: A vector of arguments given to `get_peaktable` to generate the
  peak table.

## Note

This function is adapted from
[getPeakTable](https://github.com/rwehrens/alsace/blob/master/R/getPeakTable.R)
function in the alsace package by Ron Wehrens.

## References

- Broeckling, C. D., Afsar F.A., Neumann S., Ben-Hur A., and Prenni
  J.E. 2014. RAMClust: A Novel Feature Clustering Method Enables
  Spectral-Matching-Based Annotation for Metabolomics Data. *Anal.
  Chem.* **86**:6812-6817.
  [doi:10.1021/ac501530d](https://doi.org/10.1021/ac501530d) .

- Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
  LC–DAD using multivariate curve resolution: the alsace package for R.
  *Metabolomics* **11**:143-154.
  [doi:10.1007/s11306-014-0683-5](https://doi.org/10.1007/s11306-014-0683-5)
  .

## See also

[`attach_ref_spectra`](https://ethanbass.github.io/chromatographR/dev/reference/attach_ref_spectra.md)
[`attach_metadata`](https://ethanbass.github.io/chromatographR/dev/reference/attach_metadata.md)

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
