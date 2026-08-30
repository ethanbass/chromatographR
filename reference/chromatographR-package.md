# chromatographR: Chromatographic Data Analysis Toolset

Tools for high-throughput analysis of HPLC-DAD/UV chromatographic data.
The package provides functionality for preprocessing, alignment, peak
detection and fitting, peak table construction, and visualization of
chromatographic data.

## Details

A typical workflow includes signal preprocessing and peak table
construction following general approaches described in Wehrens et al.
(2015).

Chromatogram alignment is implemented using parametric time warping
(PTW) or variable penalty dynamic time warping (VPdtw).

Peak detection is based on zero-crossings in smoothed derivatives. Peaks
are then modeled using Gaussian, exponential-Gaussian hybrid, or
bidirectional exponentially modified Gaussian peak functions via
nonlinear least squares.

Further details and example workflows are provided in the package
vignettes.

## Core analysis functions

- [`read_chroms()`](https://ethanbass.github.io/chromConverter/reference/read_chroms.html):
  Import chromatograms from a variety of vendor formats.

- [`preprocess()`](https://ethanbass.github.io/chromatographR/reference/preprocess.md):
  Preprocess chromatographic matrices.

- [`correct_rt()`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md):
  Align chromatograms.

- [`get_peaks()`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md):
  Find and integrate peaks.

- [`get_peaktable()`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md):
  Assemble peak table.

- [`attach_metadata()`](https://ethanbass.github.io/chromatographR/reference/attach_metadata.md):
  Attach metadata to peak table.

- [`attach_ref_spectra()`](https://ethanbass.github.io/chromatographR/reference/attach_ref_spectra.md):
  Attach reference spectra to peak table.

- [`normalize_data()`](https://ethanbass.github.io/chromatographR/reference/normalize_data.md):
  Normalize `peak_table` or `chrom_list`.

## Peak refinement

- [`combine_peaks()`](https://ethanbass.github.io/chromatographR/reference/combine_peaks.md):
  Combine duplicate peaks based on retention time and spectral
  similarity.

- [`merge_peaks()`](https://ethanbass.github.io/chromatographR/reference/merge_peaks.md):
  Merge split peaks into a single feature in a peak table.

- [`filter_peaktable()`](https://ethanbass.github.io/chromatographR/reference/filter_peaktable.md):
  Filter peak features in a peak table.

- [`filter_peaks()`](https://ethanbass.github.io/chromatographR/reference/filter_peaks.md):
  Filter peaks in peak lists.

## Visualization functions

### Chromatogram plots

- [`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md):
  Plot chromatograms as traces.

- [`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md):
  Plot chromatograms as heatmap.

- [`annotate_peaks()`](https://ethanbass.github.io/chromatographR/reference/annotate_peaks.md):
  Add peak labels to chromatogram plots.

- [`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md):
  Plot chromatograms as mirror plots.

- [`plot.ptw_list()`](https://ethanbass.github.io/chromatographR/reference/plot.ptw_list.md):
  Plot parametric time warping (PTW) alignment object.

### Spectral plots

- [`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md):
  Plot spectrum and/or trace of specified peak.

- [`plot_spectrum_inset()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum_inset.md):
  Plot spectrum over chromatogram.

- [`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md):
  Plot all spectra for specified peak.

- [`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md):
  Interactively extract spectra from a chromatogram.

### Peak plots

- [`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md):
  Plot fitted peaks over chromatographic trace.

- [`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md):
  Plot chromatograms and/or spectra from a peak table.

- [`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md):
  Create boxplot from a peak table object.

## Data utilities and I/O

- [`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md):
  Return retention times from a peak table or a list of chromatograms.

- [`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md):
  Return wavelengths from a peak table or a list of chromatograms.

- [`write_peaktable()`](https://ethanbass.github.io/chromatographR/reference/write_peaktable.md):
  Export peak table to CSV or XLSX format.

- [`reshape_chroms()`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md):
  Reshape a list of chromatograms to long format.

- [`reshape_peaktable()`](https://ethanbass.github.io/chromatographR/reference/reshape_peaktable.md):
  Reshape a peak table to long format.

## Example data

- `Sa`: A list of four goldenrod root chromatograms.

- `Sa_pr`: Preprocessed goldenrod root chromatograms.

- `Sa_warp`: Preprocessed and aligned goldenrod root chromatograms.

- `pk_tab`: Peak table from aligned goldenrod root chromatograms.

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

Useful links:

- <https://ethanbass.github.io/chromatographR/>

- <https://github.com/ethanbass/chromatographR/>

- Report bugs at <https://github.com/ethanbass/chromatographR/issues/>

## Author

Ethan Bass
