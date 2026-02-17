# chromatographR

Chromatographic Data Analysis Toolset

## Details

Tools for high-throughput analysis of HPLC-DAD/UV chromatograms (or
similar data). Includes functions for preprocessing, alignment,
peak-finding and fitting, peak-table construction, data-visualization,
etc. Preprocessing and peak-table construction follow the rough formula
laid out in alsace (Wehrens, R., Bloemberg, T.G., and Eilers P.H.C.,
2015.
[doi:10.1093/bioinformatics/btv299](https://doi.org/10.1093/bioinformatics/btv299)
). Alignment of chromatograms is available using parametric time warping
(PTW) (Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015.
[doi:10.1093/bioinformatics/btv299](https://doi.org/10.1093/bioinformatics/btv299)
) or variable penalty dynamic time warping (VPdtw) (Clifford, D., &
Stone, G. 2012.
[doi:10.18637/jss.v047.i08](https://doi.org/10.18637/jss.v047.i08) ).
Peak-finding relies on the algorithm suggested by Tom O'Haver in his
[Pragmatic Introduction to Signal
Processing](https://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm).
Peaks are then fitted to a gaussian or exponential-gaussian hybrid peak
shape using non-linear least squares (Lan, K. & Jorgenson, J. W. 2001.
[doi:10.1016/S0021-9673(01)00594-5](https://doi.org/10.1016/S0021-9673%2801%2900594-5)
). More details on package usage and a suggested workflow can be found
in the vignette.

## Analysis functions

- - [`read_chroms()`](https://ethanbass.github.io/chromConverter/reference/read_chroms.html):
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
    Normalize `peaktable` or `chrom_list`.

## Visualization functions

- - [`boxplot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md):
    Create boxplot from peaktable object.

  - [`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md):
    Plot chromatograms as traces.

  - [`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md):
    Plot chromatograms as heatmap.

  - [`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md):
    Plot spectrum and/or trace of specified peak.

  - [`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md):
    Plot all spectra for specified peak.

  - [`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md):
    Plot chromatograms as mirror plot.

  - [`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md):
    Plot spectrum at wavelength specified by clicking on a chromatogram.

  - [`plot.peak_list()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md):
    Plot fitted peaks over chromatographic trace.

  - [`plot.peak_table()`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md):
    Plot the trace and/or spectrum of a specified peak from the peak
    table.

  - [`plot.ptw_list()`](https://ethanbass.github.io/chromatographR/reference/plot.ptw_list.md):
    Plot PTW alignment object.

## Utility functions

- - [`combine_peaks()`](https://ethanbass.github.io/chromatographR/reference/combine_peaks.md):
    Combine duplicate peaks in peak table based on retention time and
    spectral similarity.

  - [`merge_peaks()`](https://ethanbass.github.io/chromatographR/reference/merge_peaks.md):
    Merge split peaks into a single column of a peak table.

  - [`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md):
    Return retention times from a peak table or a list of chromatograms.

  - [`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md):
    Return wavelengths from a peak table or a list of chromatograms.

  - [`reshape_chroms()`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md):
    Reshape a list of chromatograms to long format.

  - [`reshape_peaktable()`](https://ethanbass.github.io/chromatographR/reference/reshape_peaktable.md):
    Reshape a `peak_table` object to long format.

  - [`write_peaktable()`](https://ethanbass.github.io/chromatographR/reference/write_peaktable.md):
    Export peak table in `csv` or `xlsx` format.

## Example data

- - [Sa](https://ethanbass.github.io/chromatographR/reference/Sa.md): A
    list of four goldenrod root chromatograms.

  - [Sa_pr](https://ethanbass.github.io/chromatographR/reference/Sa_pr.md):
    Preprocessed goldenrod root chromatograms.

  - [Sa_warp](https://ethanbass.github.io/chromatographR/reference/Sa_warp.md):
    Preprocessed and aligned goldenrod root chromatograms.

  - [pk_tab](https://ethanbass.github.io/chromatographR/reference/pk_tab.md):
    Peak table from aligned goldenrod root chromatograms.

## See also

Useful links:

- <https://ethanbass.github.io/chromatographR/>

- <https://github.com/ethanbass/chromatographR/>

- Report bugs at <https://github.com/ethanbass/chromatographR/issues/>

## Author

Ethan Bass
