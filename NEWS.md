# chromatographR 0.7.3

* Updated `reshape_peaktable` function to include wavelength and retention time data.
* Updated `pk_tab` data file to latest peak table format.
* Fixed minor issues with roxygen formatting.
* Fixed bug in `get_peaks` causing occasional errors due to edge cases.

# chromatographR 0.7.2

* Export `get_times` and `get_lambdas` functions.
* Small revisions to documentation.
* Use `tryCatch` to allow for missing spectra in `plot_all_spectra`.
* Fixed rounding of p-values in `cluster_spectra`.
* Changed `alpha` in `cluster_peaks` to match common usage, such that the 
`alpha` parameter now specifies the significance level rather than the 
confidence level (1-alpha).
* Deprecated `peak_no` argument in `cluster_peaks` in favor of new `min_size` and `max_size` arguments.

# chromatographR 0.7.1

* Fixed bug in `get_peaks` causing peaks to erroneously filtered out in some cases.
* Made small updates to documentation (in `preprocess` and `fit_peaks` functions) to better describe arguments.
* Added warning in `mirror_plot` when `var` contains more than two levels and levels aren't specified.
* Fixed bug to allow `mirror_plot` to work properly with 2D data.
* Return dimnames for 1D ptw model objects returned by `correct_rt`.
* Minor updates to vignette.

# chromatographR 0.7.0

* Updated `correct_peaks` function so it works properly for correcting retention times in peak lists.
* Added `fixed_levels` argument to `reshape_peaktable` so features can be plotted in the order they're provided by the user.
* Added option for summing split peaks using the `merge_peaks` function by selecting `method = "sum"`.
* Updated `get_peaktable` so that the `use.cor` argument works correctly (to use corrected retention times stored in a separate column).
* Fixed `mirror_plot` so it can take numeric input for lambdas.
* Changed default setting of `verbose` argument in `correct_rt` from `FALSE`
to default setting.
* Removed `load_chroms` function. Use `read_chroms` instead.
* Eliminated spurious warning from `attach_ref_spectra` function.
* Changed name of `index` argument in `plot.peak_list` to `idx`. The original argument is now deprecated.
* Fixed bug affecting `plot_purity` argument in `plot.peak_list`.
* Fixed bug in `reshape_chroms` so empty metadata column no longer appears.
* The `plot_spectrum` function now includes the peak names when plotting spectra.
* Fixed `correct_rt` so it no longer requires user-provided `lambdas` for 1D chromatograms.
* Added `subset.peak_table` function for easily subsetting peak_tables (e.g. to exclude specific
peaks or samples).
* Added `what` argument for `plot_all_spectra` (e.g. to plot multiple spectra at a particular retention time).

#### Refactoring of `cluster_spectra` function:

* For simplicity, `cluster_spectra` now requires reference spectra to be attached to peak table.
* Accordingly, the `chrom_list` argument is no longer needed.
* Saving to RDS is now turned off by default.
* The `pvclust` package is now suggested instead of being required.

#### Updates to vignette and documentation

* Suggest numeric input to lambdas instead of character input to reduce unnecessary confusion.
* Made other minor changes to text of vignette to (hopefully) improve clarity.
* Added a short section on the attachment of reference spectra.

# chromatographR 0.6.1

* Fixed bug in plot functions (e.g. `plot_chroms` and `plot_spectrum`) causing error when retention times are inconsistent between chromatograms.
* Eliminated spurious warning in preprocess function.
* Updated `read_chroms` syntax in vignette.

# chromatographR 0.6.0

#### New features

* Enabled use of `parallel` package for parallel processing (in addition to current options using `mcapply`). (These options require the installation of suggested package `pbapply`).
* Updated `get_peaktable` for greater flexibility (e.g. for usage of 'ChemStation' peak lists as input).

#### Other changes
* Made some minor changes to vignette to improve clarity (e.g. using single wavelength for integration, etc.)

#### Bug fixes

* Fixed error in `attach_metadata` when there are NA values in merge column.

# chromatographR 0.5.6

* Fixed bug in preprocess function causing fatal error due to misrecognition of matrices.
* Fixed behavior of `plot_chroms` and `correct_rt` to allow automatic detection of `lambda` for 1D chromatograms.
* Fixed bug in `combine_peaks` (due to misplaced parenthesis).
* Added new option to filter by maximum peak area or height in `filter_peaktable` (`what = "max"`), as suggested by Katherine Holmes.

# chromatographR 0.5.5

* Fixed bug in `get_peaktable` causing failure to print strip plot when `plot_it == TRUE`.

# chromatographR 0.5.4

* Added `.zenodo.json` file.

# chromatographR 0.5.3

* Fixed bug in `plot_chroms` causing mismatched legend labels in base R plot.
* Added additional arguments to `plot_chroms` function: `xlim`,`ylim`, and `legend_position`. 
* Added additional information about arguments available in `get_peaks` for fine-tuning the peak-finding algorithm (in response to [#27](https://github.com/ethanbass/chromatographR/issues/27)).

# chromatographR 0.5.2

* Added `metadata` argument to `reshape_peaktable` for filtering metadata fields.
* Added option for renaming peaks via `reshape_peaktable` by providing a named character vector.

# chromatographR 0.5.1

* In `plot_chroms`, `show_legend` now defaults to FALSE to prevent overloading of the plot.

#### Bug fixes

* Fixed syntactical bug in `get_peaktable` when applied to gaussian peak list.
* Fixed bug caused by improper transfer of `time.units` metadata by `filter_peaks` function.
* Added default for missing `time.units` in `plot.peak_list`.

# chromatographR 0.5.0

#### New features

* Added `ggplot2` option to `plot_spectrum`, `plot.peak_table` and `plot_all_spectra` functions.
* Reworked `write_chroms` for more sensible handling of paths and added `filename` argument.
* Updated `get_purity` function to improve speed.
* Added additional argument to `reshape_chroms` function for subsetting data by 
retention times (`rts`).
* Added parallel processing through the `pbapply` package for the `correct_rt`,
`get_peaks`, and `preprocess` functions by setting the `cl` argument.

#### Other changes

* Changed behavior of `preprocess` when inferring retention times so chromatograms are no longer rounded down to the largest integer.
* In `preprocess`, spectral smoothing is no longer applied on 2D chromatograms, removing error message when preprocess is used with default settings.
* Moved position of `...` argument to end in `plot.peak_table`.
* Changed `progress_bar` argument to `show_progress` in `correct_rt`, `preprocess`
and `get_peaks` to fix strange `pmatch` behavior with additional arguments to
preprocess.
* Changed orientation of "plotly" plots generated by `plot_spectrum` to match other
plotting engines.
* Deprecated the `mc.cores` argument in `correct_rt` is now deprecated in favor of the new
`cl` argument.
* Deprecated the `parallel` argument in `preprocess` in favor of just using `cl`.
* Changed name of first argument in `mirror_plot` from `peak_table` to `x`. Otherwise the function has not changed.
* Added additional tests, improving test coverage to 80%.
* Updated `get_chrom_list` (internal) to allow parsing of subsetted lists.

# chromatographR 0.4.8

* Fixed bug in `merge_peaks` function so it works properly (to combine 2 or more
peaks in a peak table).
* Fixed bugs in `plot_chroms` preventing plotting with `ggplot2` and plotting wrong chromatograms in base R. 
* Added additional tests for `plot_chroms` and reshape functions.

# chromatographR 0.4.7

* Added `reshape_peaktable` function for conversion of peak tables to long format.
* Turned off `estimate_purity` in `get_peaks` function by default.
* Added option to filter by wavelength in `reshape_chroms`, speeding up `plot_chroms`.

# chromatographR 0.4.6

### New Features

* Added `plot_chroms` function for easily plotting multiple traces from a list of chromatograms.
* Minor changes to internal syntax of `correct_rt` to give more informative error messages.
* Added `estimate_purity` argument in `get_peaks` to toggle peak purity estimation.
* Changed default setting for `progress_bar` in `correct_rt` and `get_peaks`. Now defaults to `TRUE` if `pbapply` is installed.
* Added additional tests of utility functions and new `plot_chroms` function.
* Minor changes to vignette.
* Minor changes to documentation.

### Bug fixes

* Fixed bug causing mismatched time axes and alignment issues after VPdtw warping
(again), so that it returns matrices with a consistent time axis.
* Fixed y unit label in `boxplot.peak_table` function.
* Fixed behavior of `plot_spectrum` so spectrum is exported properly when `engine == plotly`.
* Fixed bug in `write_peaktable` when writing to `xlsx`.

# chromatographR 0.4.5

### New Features

* Added `reshape_chroms` function for converting chromatograms to "long" format.
* Added `write_peaktable` function to easily write peak_table to `csv` or `xlsx`.
* Added `get_purity` function for assessing peak purity.
* Allow multiple peaks as arguments to `plot.peaktable`.
* Added functions for plotting traces and spectra with [plotly](https://plotly.com/r/):
`plotly_trace` and `plotly_spec`.
* Fixed `preprocess` so it will no longer try to interpolate along columns for 2D data.
* Added stand-alone `boxplot` function for `peak_table` objects.
* Added a new class (`ptw_list`) and plotting function for lists of `ptw` alignment objects.
* Added `plot_it` argument in `correct_rt` for plotting alignments.
* Added [VPdtw](https://github.com/ethanbass/VPdtw/) as a dependency (instead of being only suggested).
* Added `progress_bar` option to `get_peaks` and `correct_rt`.
* Improved error handling in `plot.peaklist`.
* Updated find_peaks function with more and better smoothing options to improve peak-finding. Now defaults to gaussian smoothing.
* Changed `fit_peaks` function and syntax (see below).
* Minor updates to vignette.

#### Changes to *fit_peaks* function:

* Simplified logic in `fit_peaks` function.
* Modified `fit_peaks` syntax so it now takes a matrix (`x`) and a wavelength
(`lambda`) instead of a numeric vector (`y`).
* Incorporated assessment of peak purity during peak fitting.
* Added wavelength (`lambda`) to `peak_list` and `peak_table` metadata.
* Fixed bug to allow fitting of a single peak with `fit_peaks`.

### Bug fixes

* Fixed bug causing mismatched time axes (and improper alignment of chromatograms) after variable penalty dynamic time warping (VPdtw). 
* Fixed bug in `attach_metadata` that could result in disordered rows.
* Fixed occasional test failure on MKL server by skipping `cluster_spectra` test on CRAN.
* Adjusted `cluster_spectra` and `combine_peaks` functions so messages can be suppressed with `verbose == FALSE`.

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
