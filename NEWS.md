# chromatographR 0.8.0

### Breaking changes

#### Preprocessing

* **`preprocess` no longer clamps negative values to zero by default.** Set the
new `zero_floor = TRUE` for bit-identical previous output.
* **Deprecated `preprocess`'s `remove.time.baseline` argument in favor of
`baseline = "none"`.**
* **`maxI` now follows `...` in `preprocess`'s signature, so it must be named in
full.** Partial matching otherwise captured `m`, which belongs to the baseline
estimators. `maxI` can no longer be given positionally.
* **`preprocess` now validates `dim1` and `dim2`, which must be strictly
increasing.** A named `by`, `length.out` or `along.with` element is reported
specifically, since `c(from, to, by = x)` is easily written for `seq(from, to,
by = x)`.

#### Retention time correction

* Added `fill_zeros` to `correct_rt`, controlling whether out-of-bounds regions
produced by warping are filled with zeros or returned as `NA`. Defaults to
`FALSE`, which changes VPdtw behavior.

#### Peak detection and integration

* Corrected an off-by-one error in `find_peaks` that reported every apex one
timepoint early. Retention times from `get_peaks` and `get_peaktable` shift
accordingly, and peaks whose preceding timepoint had been clipped to zero by
`preprocess` are now detected.
* `plot_spectrum` and `plot_all_spectra` now label the correct spectral
maxima, having inherited the same off-by-one error. Reported lambda-max values
shift by one grid step.
* `find_peaks` now rounds the `"mva"` and `"tmva"` smoothing windows up to an
odd width of at least three points, which changes which apices are detected
for these methods.
* `find_peaks` now errors when the provided `smooth_window` is narrower 
than seven timepoints when Savitzky-Golay (`"savgol"`) smoothing is
selected, rather than silently returning the signal unsmoothed.
* Renamed `find_peaks`'s `slope_thresh` and `amp_thresh` arguments to
`min_curvature` and `min_amplitude`. `slope_thresh` was a misnomer
because the argument thresholded the second derivative, not the slope. 
For now, the old names still work with a deprecation warning.
* `get_peaks` and `fit_peaks` no longer filter on peak width by default.
`sd_max` used to drop peaks with sd > 50 timepoints. Its replacement,
`max_fwhm`, defaults to no filtering; set it to restore the filter.

#### Peak purity

* `get_purity` now decides which timepoints belong to a peak from the
wavelength of maximum absorbance at the apex, rather than the first wavelength
in the chromatogram.

#### Peak tables and plotting

* Corrected `get_peaks` to return model-fitted peak heights rather than raw peak
heights. This changes numerical output for workflows relying on the previous
behavior.
* Corrected `merge_peaks` so the result no longer depends on the order in which
`peaks` are supplied.
* Added `summarize_by` argument to `get_peaktable` controlling how peak metadata
are summarized within each cluster. The new default, `"weighted.mean"`, weights
peaks by intensity; `"mean"`, `"median"` and `"max"` are also available. Set
`summarize_by = "mean"` to reproduce previous behavior.
* Retention times in dimnames are now formatted to a consistent number of
decimal places, so an interpolated time axis is no longer ragged. Lookups by the
old names need updating (`x["10.10", ]` rather than `x["10.1", ]`).
* Added automatic thinning of chromatograms in `plot_chroms` to speed up
plotting of raw data, adjustable with the new `time_resolution` argument.
* Changed behavior of `plot_chroms` so the legend is no longer displayed by
default.
* When `export_spectrum == TRUE`, `plot_spectrum` and `plot_all_spectra` now
return spectra invisibly.
* Changed the name of the `a` parameter in `plot.peak_list` to `alpha`.
* Renamed the `chrom_list` argument to `x` in `plot_chroms_heatmap`, which now
also accepts a peak table in the first position.
* Removed `index` and `chr` arguments across plotting functions, previously
deprecated in v0.7.0. Use `idx` instead.

#### Other renamed and deprecated arguments

* **Renamed the remaining dot-separated argument names to snake_case**:
`spec.smooth` to `smooth_spectra`, `min.cor` to `min_cor`, and
`n.traces`/`n.zeros` to `n_traces`/`n_zeros`. The old names still work, with a
warning. `init.coef` and `max.only` are retained, since they are passed through
to functions that define them.
* Deprecated arguments in `get_peaks` still using dot notation (e.g., `sd.max`).
Use the snake_case equivalents instead, or `max_fwhm` where the snake_case name
is itself deprecated below.
* Deprecated `sd_max` in `get_peaks` and `fit_peaks`, and `min_sd`/`max_sd` in
`filter_peaks`, in favor of `min_fwhm` and `max_fwhm`. A standard deviation is
equivalent to peak width only for `fit = "gaussian"`. See `?filter_peaks`.
* Renamed the Whittaker smoothing penalty passed to `preprocess` and
`correct_baseline` from `lambda` to `stiffness`, reserving `lambda` for
the common usage to specify a wavelength used throughout the package.
* Deprecated `peak` argument in `plot_all_spectra` in favor of `loc`.
* Deprecated the `save` argument in `cluster_spectra` in favor of `outfile`.
* Removed deprecated `peak_no` argument from `cluster_spectra`.

### New features

#### New functions

* Added `correct_baseline` function for estimating and subtracting a baseline
from a `chrom_list` along the time axis. Unlike `preprocess` it does not
interpolate, smooth or rescale, so it can be applied to data that are already
preprocessed or aligned.
* Added `get_baseline` and `plot_baseline` for retrieving and plotting the
baseline retained by `preprocess` or `correct_baseline` when
`attach_baseline = TRUE`.
* Added `crop_chroms` function for restricting chromatograms to a range of
retention times (`rt_range`) and/or wavelengths (`lambda_range`). Cropping is
pure subsetting, so values and dimnames pass through unchanged, unlike
`preprocess(dim1 = )`, which interpolates onto a new axis.
* Added `subtract_blanks` function for subtracting blank chromatograms from a
`chrom_list`.
* Added `correct_rt_group` function for group-wise warping function estimation
using within-group averages. This improves alignment when samples share
batch-level shifts or contain group-specific peaks, and is generally faster than
per-sample alignment.
* Added `plot_spectrum_inset` function for plotting UV spectra overlaid onto a
chromatogram.
* Added `annotate_peaks` function for adding text labels to peaks in plots
created by `plot_chroms`.
* Added `get_spectra` function for extracting spectra from chromatograms. The
`export_spectrum` argument in `plot_spectrum` and `plot_all_spectra` is
deprecated in favor of `get_spectra`.
* Added `get_time_resolution` for reading a chromatogram's sampling interval,
in the units of its time axis by default or in minutes via `unit`.

#### New generic functions

* Added `summary` method for `peak_table` objects, reporting dimensions, retention time range, intensity statistics, and metadata availability.
* Added `print.peak_table` function to return just the peak data.
* Added `[` method for `peak_table` objects for simultaneous subsetting of `tab`, `pk_meta`, `sample_meta`, and `ref_spectra` elements.
* Added `[` method for `peak_list` objects to preserve metadata attributes when subsetting.
* Added `print` method for `peak_list` objects, summarizing samples, wavelengths, fit method, and peak counts.

#### Preprocessing

* Preprocess can now take nested lists containing a `dad` element in each list. This change allows better compatibility with chromConverter.

#### Baseline correction

* Added a `baseline` argument to `preprocess` selecting the estimator used
along the time axis: `"als"` (still the default), `"arpls"`, `"aspls"` or
`"snip"`. See `?baseline_estimators` for what each does and where each fails.
`"snip"` requires `m`, the filter half-width in timepoints, and has no default.
* The smoothing parameter each estimator used is recorded on the result as a
`"baseline_param"` attribute, carrying the estimator's name.
* Added an `attach_baseline` argument to `preprocess`, retaining the estimated
baseline as a `"baseline"` attribute so a poor fit can be inspected or plotted.
It is dropped by any transform that invalidates it (`correct_rt`,
`subtract_blanks`), since it is indexed by the time axis.
* `correct_baseline` now accepts a single trace as a plain numeric vector, and
returns one, carrying names through to the result and its baseline.
* Baseline estimators now require a complete trace, raising an informative error
when the signal contains missing or non-finite values rather than failing
obscurely or silently splicing the two sides of a gap together.

#### Retention time correction

* Refactored `ptw` warping in `correct_rt` using whole-matrix interpolation resulting in substantial speedup.

#### Peak detection and integration

* Added support for bidirectional exponentially modified gaussian (BEMG) peak
shape (`fit = "bemg"` in `get_peaks`/`fit_peaks`). The default peak shape
remains `"egh"` for now but this may change in a future release.
* Added support for fitting a local baseline beneath each peak through the
`baseline` argument in `get_peaks`: `"none"` (default), `"flat"` or `"sloped"`.
Peak heights and areas are computed above the fitted baseline.
* Added a `min_prominence_snr` argument to `find_peaks`, filtering apices by
topographic prominence in multiples of the estimated noise. Unlike
`min_curvature` and `min_amplitude` this threshold is scale-free, so it
transfers between datasets. `find_peaks` now also returns a `prominence`
column.
* Added manual integration of fixed retention-time windows to `get_peaks`
and `fit_peaks` via a `windows` argument, quantifying a compound even where
no peak is detected. `get_peaktable` gives one column per window. See
`?get_peaks`.

#### Fit diagnostics

* Two fit-diagnostic columns are now reported by `get_peaks`: `converged`
(flagging peaks that didn't converge) and `window_limited` (flagging peaks whose
reported width is at the limit defined by the integration window).
* Fit diagnostics can now be used for filtering: `filter_peaks` gains `converged` and
`window_limited` arguments, and `filter_peaktable` gains `min_converged`. All three 
error on `fit = "raw"`, which reports no fit diagnostics.
* `get_peaks` now warns when more than a fifth of the peak list is integration
windows too small to fit, which usually means `find_peaks` is detecting noise.
Suppressed for `fit = "raw"` and when `windows` is supplied.


#### Peak tables and plotting

* Added support for plotting spectra from multiple chromatograms in `plot_spectrum` by supplying a vector to `idx` (e.g. `idx = c(1, 2, 4)`). Overlaid traces and spectra are color-coded by sample.
* Allow `plot_chroms` and `plot_chroms_heatmap` to accept peak table as first argument.
* Added new options to `combine_peaks` via the `choose` argument:
`"least_sparse"` retains the peak detected in the most samples, and `"lambda"`
retains the peak measured at a preferred wavelength.
* Added `max_zeros` argument in `filter_peaktable` to exclude features above a threshold proportion of zero values.
* Added support for numeric input to lambda argument in `reshape_chroms`.

#### Other

* Refactored `normalize_data` so it detects metadata and peak table columns
automatically; added an `on_invalid` argument controlling how zero, negative and
`NA` normalization values are handled; and improved performance.
* Added support for new `chrom_list` class introduced in chromConverter v0.9.0.

### Bug fixes and other minor changes

#### Preprocessing and baseline correction

* `preprocess` no longer fails on chromatograms with a few discontinuous wavelengths ---
`smooth.spline` requires at least four unique values, so spectral smoothing is
now skipped below four wavelengths.
* Fixed input validation in `preprocess` so lists containing non-matrix elements no longer pass validation without error.
* Fixed `preprocess` function incorrectly truncating chromatograms when outliers are present and `dim1` is not specified.

#### Peak detection and integration

* `find_peaks(bounds = TRUE)` now returns its full four-column schema when no
peaks are detected, rather than dropping `lower` and `upper`.
* `max_iter` now controls how long the peak fitters run. Previously it was
ignored, leaving `nls.lm`'s default of `100 * (parameters + 1)` evaluations.
`"egh"` and `"bemg"` fits may shift slightly.
* `get_peaks` no longer returns incomplete results when a parallel
worker dies. `mclapply` reports a crash as `NULL`, which used to reach the
returned `peak_list` unflagged and fail later in `get_peaktable`. Lost samples
are now refitted and unrecoverable ones raise an error.
* `get_peaks` no longer errors when a sample contains no peaks; previously a
single empty sample failed the whole call. `fit_peaks` now returns an empty
`data.frame` rather than a one-row matrix of `NA`s.
* `FWHM` is now computed as `2*sqrt(2*log(2))*sigma` rather than the rounded
`2.35*sigma`, and a bug in its calculation for raw peak shapes is fixed.
* Suppressed spurious "essentially perfect fit" warning from `summary.lm` when computing R² in `fit_peaks`.

#### Peak tables and filtering

* Added error when duplicated sample names are provided to `get_peaktable` or `get_peaks` since this will cause downstream issues.
* Fixed `get_peaktable` (with `clust = "sp.rt"`) to no longer take the absolute value of spectral similarity scores, which caused negatively correlated spectra to be treated as similar.
* Fixed `filter_peaktable` and `filter_peaks`: `rts` with a non-zero `tol`
matched every peak within the tolerance rather than the nearest one, so a single
requested time could select most of the table. Also, criteria are now inclusive
(`min_rt = 5` retains `rt == 5`), and peaks with `NA` in some samples are no
longer excluded by `min_value`.
* Fixed `subset.peak_table()`, where `sample_meta` was subset on the wrong
condition.
* Fixed bug in `normalize_data` causing loss of sample names when normalizing chromatograms.
* Fixed error in equation describing spectral clustering in documentation for `get_peaktable`.

#### Plotting

* Fixed `reshape_chroms` so it no longer errors when `lambdas` is left at its
documented `NULL` default. It now returns every wavelength in this case.
* Fixed `reshape_chroms` so it returns a factor for the `sample` column.
* Fixed and standardized plotting details: single chromatograms are now black in
plotly and ggplot; the `base` heatmap colorbar renders correctly;
`plot_spectrum` keeps its plotly axis labels; axis labels and trace titles are
consistent across engines; and `plot_all_spectra` always exports wide format
regardless of the plotting engine.

#### Documentation and packaging

* Made substantial revisions to documentation to improve clarity.
* Formalized requirement for R version 4.1.0.
* Fixed issue with visual tests due to new version of Kaleido.
* Updated the GC-FID vignette to fetch its example data from a mirror on GitHub, since Dryad now requires authentication for all downloads.

# chromatographR 0.7.5

* Fixed bug in `check_chr` when multiple chromatograms are provided.

# chromatographR 0.7.4

* Added new heatmap plot in `plot_chroms` for plotting a list of chromatograms as a heatmap.
* Added labeling of spectral peaks for `ggplot2` and `plotly` plots in `plot_spectrum`.
* Added minimum amplitude for labeling of spectral peaks in `plot_spectrum`, such that peaks are given a text label only if they exceed 5% of the maximum signal intensity.
* Removed `parallel` and `mc.cores` argument from `preprocess` (deprecated in `v0.5.0`).
* Added `title` argument in `plot_chroms`.
* Added option for different plotting engines through `engine` argument in `plot.ptw_list`.
* Added option for heatmap plots through `what` argument in `plot.ptw_list`.
* Changed behavior of `reshape_chroms` such that it can be applied to 2D chromatograms without throwing an error.
* Added index of functions to the package landing page.
* Made minor updates to function documentation.

#### Updates to GC-FID vignette

* Changed alignment plots to heatmaps rather than traces.
* Hid `ggordiplot` function in "folded" code block
* Added alt text to images.
* Updated call to `adonis2` in accordance with new defaults in `vegan v2.7-0`.
* Made other minor changes to text.

#### Updates to chromatographR vignette

* Added figure demonstrating comparison of alignments with new `plot_chroms_heatmap` function.
* Added explanatory figure captions.

# chromatographR 0.7.3

* Updated `reshape_peaktable` function to include wavelength and retention time data.
* Updated `pk_tab` data file to latest peak table format.
* Fixed minor issues with roxygen formatting.
* Added functional "families".
* Fixed bug in `get_peaks` causing occasional errors due to edge cases.
* Added titles for ggplot2 and plotly `plot_spectrum` plots.
* Fixed bug in console output when providing chromatogram index as string to `plot_spectrum`.
* Changed default `ggplot2` theme to `theme_light`.

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
* Deprecated `chr` argument in multiple plotting functions (e.g., `plot_spectrum`) in favor of new `idx` argument.
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
* Added parallel processing through the `pbapply` package for the `correct_rt`, `get_peaks`, and `preprocess` functions by setting the `cl` argument.

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
