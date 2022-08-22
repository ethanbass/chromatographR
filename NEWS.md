# chromatograph 0.4.3

* Minor changes to documentation.
* Added additional check of chrom_list dimensions and names.

# chromatograph 0.4.2

* Added option to select `time.units` for peak area in `get_peaks` function
facilitating better comparison with vendor software.
* Now allow preservation of instrumental metadata through preprocessing and alignment steps.
* Added `filter_peaktable` function.
* Changed default behavior in `correct_rt` to `corrected_values` rather than `models`.
* Fixed bug in `cluster_spectra` affecting peaks with 0 standard deviation.
* Fixed bug affecting `peak_list` metadata.
* Deprecated `load_chroms` function. Please use `read_chroms` from chromConverter
to import files instead.

### Minor changes:

* Added more informative warnings and error messages to various functions.
* Now recommend installation of `VPdtw` from CRAN instead of https://ethanbass.github.io/drat/
* Fixed typos in vignette

# chromatographR 0.4.1

* Extended package DESCRIPTION and added citations to relevant references.
* Added `\value` and `\section{Side effects}` fields to docs for the various plot
functions. 
* Fixed bug in `mirror_plot` so legend can be fully hidden by setting `plot_legend`
to FALSE.
* Fixed bug in `scan_chrom` so additional arguments are passed to `plot_spectrum`.
* Added `color` argument to customize color of fitted peaks in `plot.peak_list`.
* Fixed plot functions and examples so they don't change par settings.
* Other small updates to documents (mostly formatting or small clarifications).
* Fixed `plot.peak_table` so it can return spectra if `export_spectrum` is TRUE.
* Added error for `box_plot` option in `plot.peak_table` if metadata is not attached.
* Added error for `box_plot` option in `plot.peak_table` if peak is not provided to `loc`.
* Added error in `plot_spectrum` function for user-supplied retention times beyond
the edges of the chromatogram.
* Added error in `plot_spectrum` function for unspecified `lambda` (if `what=="click`).
* Released on CRAN

# chromatographR 0.4.0

* Added support for variable penalty dynamic time warping (VPdtw) through
`correct_rt` function.
* Fixed bug in get_peaks function.
* Allow preprocessing without interpolation.
* Fixed bug so preprocess can work on Windows (without parallel processing).
* Allow use of raw data for peak integration in `get_peaks`.
* Added `verbose` option to `correct_rt` to print reference chromatogram.

# chromatographR 0.3

* Added a `NEWS.md` file to track changes to the package.
