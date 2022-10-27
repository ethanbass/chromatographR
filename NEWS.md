# chromatographR 0.4.5

#### New features
* Added `reshape_chroms` function for converting chromatograms to "long" format.
* Added `export.peaktable` function to easily write peak_table to `csv` or `xlsx`.
* Added functions for assessing peak purity: `get_purity_values` and `get_mean_purity`.
* Added `fill_gaps` function for filling gaps in `peak_table`.
* Allow multiple peaks as arguments to `plot.peaktable`.
* Added functions for plotting trace and spectra with [plotly](https://plotly.com/r/):
`plotly_trace` and `plotly_spec`.
* Fixed `preprocess` so it will no longer try to interpolate along columns for 2D data.
* Added stand-alone boxplot function for peak_table object.
* Fixed bug in `attach_metadata` that could result in disordered rows.

#### Changes to *fit_peaks* function:
* Simplified logic in `fit_peaks` function.
* Modified `fit_peaks` syntax so it now takes a matrix (`x`) and a wavelength
(`lambda`) instead of a numeric vector (`y`).
* Incorporated assessment of peak purity during peak fitting.
* Added wavelength (`lambda`) to `peak_list` and `peak_table` metadata.
* Fixed bug to allow fitting of a single peak with `fit_peaks`.

# chromatographR 0.4.4

* Fixed issue with tests when run on certain machines (MKL).

# chromatographR 0.4.3

* Minor changes to documentation.
* Added additional check of chrom_list dimensions and names.

# chromatographR 0.4.2

#### New features
* Added option to select `time.units` for peak area in `get_peaks` function
facilitating better comparison with vendor software.
* Now allow preservation of instrumental metadata through pre-processing and alignment steps.
* Added `filter_peaktable` function.

#### Minor changes:
* Deprecated `load_chroms` function. Please use `read_chroms` from chromConverter
to import files instead.
* Changed default behavior in `correct_rt` to `corrected_values` rather than `models`.
* Added more informative warnings and error messages to various functions.
* Now recommend installation of `VPdtw` from CRAN instead of https://ethanbass.github.io/drat/
* Fixed typos in vignette

#### Bug fixes
* Fixed bug in `cluster_spectra` affecting peaks with 0 standard deviation.
* Fixed bug affecting `peak_list` metadata.

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

# chromatographR 0.3.0

* Added a `NEWS.md` file to track changes to the package.
