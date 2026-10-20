# chromatographR 0.9.0

### Breaking changes

#### Preprocessing

* `preprocess` no longer clamps negative values to zero by default. Set the
new `zero_floor = TRUE` for bit-identical previous output.
* `subtract_blanks` no longer clamps negative values to zero by default. Set
`zero_floor = TRUE` for the previous output; it now warns when it zeroes more
than 5% of a sample.
* `subtract_blanks` now matches blanks to samples by retention time,
interpolating across differing grids. Previously a blank of another length
errored, and one of equal length was subtracted row by row, even when its times
were offset.
* Deprecated `preprocess`'s `remove.time.baseline` argument in favor of
`baseline = "none"`.
* `maxI` now follows `...` in `preprocess`'s signature, so it must be named in
full. Partial matching otherwise captured `m`, which belongs to the baseline
estimators. `maxI` can no longer be given positionally.
* `preprocess` now validates `dim1` and `dim2`, which must be strictly
increasing. A named `by`, `length.out`, or `along.with` element gets an
informative error pointing at `seq()`.

#### Retention time correction

* Inverted the meaning of `correct_rt`'s `penalty` argument for
`alg = "vpdtw"`, so a higher value now allows less warping. The default
changed from `5` to `0.2`; to convert an explicit value, take its reciprocal.
Values above `1` now warn.
* The `trwdth` argument in `correct_rt` now specifies a width in minutes, not a count of timepoints.
Multiply old values by the sampling interval in minutes to convert it; a value
wider than the chromatogram now errors.
* `correct_rt` now passes `trwdth` to `ptw`, which previously ignored it.
The default of `0.2` matches `ptw`'s own `20` timepoints at `dt = 0.01`.
* `max_shift` defaults to `1` and `dilation_span` to `0.3`, so VPdtw alignments
change. See `?correct_rt`.
* `correct_rt(scale = TRUE)` now scales chromatograms only to fit the warp, so
corrected values keep their original units. Scaling divides by the maximum
instead of rescaling to [0, 1]. Use `preprocess(maxI = 1)` to rescale the data.
* `correct_rt` and `correct_rt_group` now error when samples' time axes differ
by a sampling interval or more. Previously they warped such samples by timepoint
index. Interpolate onto a common grid with `preprocess(dim1 = ...)` first.
* `correct_rt` and `correct_rt_group` now error on duplicate sample names, as
`get_peaks` already did.

#### Peak detection and integration

* `fit = "bemg"` now reports `tau_left` and `tau_right` as tail time constants
in `time_unit`, so a larger value means a longer tail. In 0.8.0 they were rates
(1/time), with left and right swapped.
* Corrected an off-by-one error in `find_peaks` that reported every apex one
timepoint early. Retention times from `get_peaks` and `get_peaktable` shift
accordingly, and peaks whose preceding timepoint had been clipped to zero by
`preprocess` are now detected.
* Fitted retention times can now fall between timepoints, so `rt` can shift
by up to one timepoint and `"bemg"` no longer understates `height`. `correct_peaks` no longer rounds `rt.cor` for ptw or lptw models.
Match `rt` to the time axis with a tolerance, not with `==` or `match()`.
* With `baseline = "sloped"`, the fitted local baseline no longer rises above
the signal at the edges of the peak window. Previously it could, underestimating
`"egh"` and `"gaussian"` areas.
* `r.squared` in `get_peaks` and `fit_peaks` output is now `1 - SS_res / SS_tot`,
not the squared correlation between the fit and the data. It is lower for poor
fits and can be negative, so revisit any cutoffs.
* `plot_spectrum` and `plot_all_spectra` now label the correct spectral
maxima, having inherited the same off-by-one error. Reported lambda-max values
shift by one grid step.
* `find_peaks` now rounds the `"mva"` and `"tmva"` smoothing windows up to an
odd width of at least three points, which changes which apices are detected
for these methods.
* `find_peaks` now errors when `smooth_window` is narrower than seven
timepoints under Savitzky-Golay (`"savgol"`) smoothing. Previously it returned
the signal unsmoothed.
* `find_peaks` and `get_peaks` now smooth over 3 timepoints by default. The
old default, `smooth_window = 0.001`, was a proportion of the run, so the
smoothing varied with run length. Set it explicitly to restore the old
behavior; see *Choosing the smoothing window* in `?find_peaks`.
* Renamed `find_peaks`'s `slope_thresh` and `amp_thresh` arguments to
`min_curvature` and `min_amplitude`. `slope_thresh` was a misnomer
because the argument thresholded the second derivative, not the slope.
For now, the old names still work with a deprecation warning.
* `get_peaks` and `fit_peaks` no longer filter on peak width by default.
Previously `sd_max` dropped peaks with sd > 50 timepoints; to filter by
width, call `filter_peaks(max_fwhm = )` on the result.

#### Peak purity

* Overhauled `get_purity`, correcting the purity ratio, so purity values change
for every peak. `cutoff` now applies to the background-corrected trace, as a
proportion of peak height rather than of total absorbance. See `?get_purity`.
* Added arguments to `get_purity`: `background_correction`, `reference`, `metric`,
`min_run` and `max_threshold` control the purity calculation, and
`noise_method`, `noise_quantile` and `noise_window` control noise estimation.
* Deprecated the `weight` argument to `get_purity` and `get_agilent_threshold`
in favor of `noise_inflation`.
* Removed the `noise_threshold` argument from `get_peaks` and `fit_peaks`,
which had no effect. Pass purity options as `purity_args` instead.
* `fit_peaks` now defaults to `estimate_purity = FALSE`, matching `get_peaks`,
so its `purity` column is `NA` unless you set `estimate_purity = TRUE`.
* `plot_purity = TRUE` in `plot.peak_list` now draws dissimilarity and
threshold curves. It also accepts `"curves_local"` and `"ratio"`; see
`?plot.peak_list`.

#### Peak tables and plotting

* Retention times in dimnames are now formatted to a consistent number of
decimal places, so an interpolated time axis is no longer ragged. Lookups by the
old names need updating (`x["10.10", ]` rather than `x["10.1", ]`).
* `plot.peak_list`'s `y.offset` is now a fraction of the plot height, with a
default of `0.03`. Previously it was in absorbance units, with a default of `25`.
Values above 1 warn. Boundary ticks are scaled the same way.
* `get_peaktable`'s `hmax` and `sigma_t`, and `combine_peaks`'s `tol`, are now
in minutes, regardless of the time unit of the data. Results change only for
chromatograms whose times are not in minutes.
* `reshape_chroms` now reads `time_resolution` in minutes, matching
`plot_chroms`. Previously it was in the units of the time axis, so divide old
values by 60 for chromatograms in seconds.
* Removed `get_peaktable`'s `out` argument, so peak tables are always
data.frames. Use `as.matrix(pt$tab)` for a matrix.
* Peak tables built from corrected peaks now name the `pk_meta` rows `rt`,
`start` and `end` rather than `rt.cor`, `start.cor` and `end.cor`, and so does
`reshape_peaktable`. Check `args$use_corrected_rt` to see which times were
used.
* `attach_metadata` now adds to an existing `sample_meta`, replacing only the
columns of the same name. Previously it replaced the whole table. To do that
now, set `peak_table$sample_meta <- NA` first.

#### Example data

* Regenerated the example datasets: `Sa_pr` and `Sa_warp` are no longer clipped at
zero, and most `pk_tab` peak names have shifted; for example, `V13` is now `V11`.

#### Other renamed and deprecated arguments

* Renamed dot-separated argument names to snake_case:
`spec.smooth` to `smooth_spectra`, `min.cor` to `min_cor`,
`n.traces`/`n.zeros` to `n_traces`/`n_zeros`, `sigma.t`/`sigma.r` to
`sigma_t`/`sigma_r`, and `use.cor` to `use_corrected_rt`. The old names still
work, with a warning. `init.coef` and
`max.only` keep their names, since they are passed through to external
functions that define them.
* Deprecated `sd_max` in `get_peaks` and `fit_peaks` in favor of
`filter_peaks(max_fwhm = )`, and `min_sd`/`max_sd` in `filter_peaks` in favor
of `min_fwhm` and `max_fwhm`. A standard deviation is
equivalent to peak width only for `fit = "gaussian"`. See `?filter_peaks`.
* Renamed the Whittaker smoothing penalty passed to `preprocess` and
`correct_baseline` from `lambda` to `stiffness`, reserving `lambda` for
the common usage to specify a wavelength used throughout the package.
* Removed the `lambda` argument from `get_spectra`, which had no effect.
* Renamed `correct_rt`'s `maxshift` to `max_shift`, now in minutes. `maxshift`
still works, read as a count of timepoints as before, with a warning.

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
* Added `match_peaks` function for matching peaks across chromatograms using
spectral correlation and retention time proximity. The resulting
`matched_peak_list` can be passed to `correct_rt` for landmark-based PTW or to
`get_peaktable`. See `?targeted_analysis` for matching to a table of targets.
* Added `get_time_resolution` for reading a chromatogram's sampling interval,
in the units of its time axis by default or in minutes via `unit`.
* Added `plot_chrom_contour` for plotting one chromatogram across retention
time and wavelength, on a compressed color scale by default.

#### New generic functions

* Converted `get_peaktable` to an S3 generic. The default method keeps the original behavior, including for unclassed lists. The new `get_peaktable.matched_peak_list` method builds a peak table directly from matched peak pairs.
* Added a `plot` method for `chrom_list` objects from `read_chroms`, drawing
traces, a heatmap or a contour plot according to `type`.

#### Baseline correction

* Added a `baseline` argument to `preprocess` choosing the time-axis baseline
estimator: `"als"` (the default), `"arpls"`, `"aspls"` or `"snip"`. See
`?baseline_estimators`.
* The smoothing parameter each estimator used is recorded on the result as a
`"baseline_param"` attribute, carrying the estimator's name.
* Added an `attach_baseline` argument to `preprocess`, keeping the estimated
baseline as an attribute for inspection. It is removed by `correct_rt` and
`subtract_blanks`, whose output it would no longer match.
* Added reference-wavelength subtraction to `preprocess`: `reference_wavelengths`
removes an offset shared by every wavelength before the baseline is estimated.
* Baseline estimators now require a complete trace, and give an informative
error when the signal contains missing or non-finite values. Previously a gap
either caused an obscure error or had its two sides spliced together.

#### Retention time correction

* Added landmark parametric time warping (`alg = "lptw"`) to `correct_rt` and
`correct_rt_group`, which fits the warp to matched peak positions rather than
whole signals. See `?correct_rt`.
* Added a `dilation_span` argument to `correct_rt`: the width, in minutes, of
the moving-maximum window used to build the `"vpdtw"` penalty, previously
hardcoded to 350 timepoints. A span wide enough to make the penalty constant
now warns.
* `correct_rt_group` can now return its warping models, one per sample, with
`what = "models"`, for use with `correct_peaks`.
* `correct_peaks` now accepts `"lptw"` and `"vpdtw"` models from `correct_rt`,
not only `"ptw"`.
* Added a `chrom_list_aligned` argument to `correct_peaks`, naming aligned
chromatograms for later spectrum and trace lookups.

#### Peak detection and integration

* Added a `min_prominence_snr` argument to `find_peaks`, filtering apices by
topographic prominence in multiples of the estimated noise. Unlike
`min_curvature` and `min_amplitude` this threshold is scale-free, so it
transfers between datasets. `find_peaks` now also returns a `prominence`
column.
* Added manual integration of fixed retention-time windows to `get_peaks`
via a `windows` argument, quantifying a compound even where no peak is
detected. `get_peaktable` gives one column per window. See
`?get_peaks`.
* Added `targets`, `tol` and `select_by` arguments to `get_peaks` for fitting
only the peaks near target retention times, which is faster than fitting every
peak. `get_peaktable` groups the result by target. See `?targeted_analysis`.

#### Fit diagnostics

* Two fit-diagnostic columns are now reported by `get_peaks`: `converged`
(flagging peaks that didn't converge) and `window_limited` (flagging peaks whose
fitted shape sits at a limit set by the integration window or the sampling
interval).
* Fit diagnostics can now be used for filtering: `filter_peaks` gains `converged` and
`window_limited` arguments, and `filter_peaktable` gains `min_converged`. All three 
error on `fit = "raw"`, which reports no fit diagnostics.
* `get_peaks` now warns when more than a fifth of the peak list is integration
windows too small to fit, which usually means `find_peaks` is detecting noise.
Suppressed for `fit = "raw"` and when `windows` is supplied.

#### Peak purity

* Purity estimation in `get_peaks` is now much faster, so it is no longer a
bottleneck.
* Purity estimation now warns when a chromatogram appears to be clipped at
zero, which makes peaks look less pure. Use `zero_floor = FALSE`, the default
in both `preprocess` and `subtract_blanks`.

#### Peak tables and plotting

* Added `targets`, `tol` and `pad` arguments to `plot.peak_list` for zooming
in on the peaks near given retention times, rescaling both axes to the
selected window.
* `plot.peak_list` now draws axes and titles each plot with its sample,
wavelength and peak model. Pass `xaxt = "n"`, `yaxt = "n"` or `main = ""` to
remove them.
* `plot.peak_list` now accepts `numbers = "rt"`, labeling peaks with their
retention times, which can be passed to `targets` to zoom in on them.
* Added a `what` argument to `plot_chroms` for plotting one element of nested
chromatograms from `read_chroms`, such as `"dad"` or a pump pressure trace.
* Added a `zlim` argument to `plot_chroms_heatmap`, fixing the range of its
color scale so that several plots can share one scale. By default the scale now
spans the plotted window in every engine; previously the `"ggplot"` and
`"plotly"` engines used the whole run when `xlim` was set.
* Added `format = "wide"` to `reshape_peaktable`, returning one row per sample
with metadata and one column per selected peak.
* `annotate_peaks` now accepts retention times through the `loc` argument,
labeling the highest point within `tol` of each supplied retention time
without requiring a peak table. See `?annotate_peaks`.
* `get_peaktable` now fills `sample_meta` with any fields added to the
chromatograms with `chromConverter::add_metadata` (chromConverter 0.11.0 or
later).
* `filter_peaks` and `filter_peaktable` now record each call as code, which
`print` shows for peak lists and `get_peaktable` stores in `args$filters`. See
`?filter_peaks`.
* `get_peaktable` now records how the peaks were fitted (`fit`, `baseline` and
`max_iter`) in `args`, and keeps any `windows` or `targets` as attributes,
including targets found in no sample. See `?get_peaktable`.
* `cluster_spectra` now numbers and colors each box on the dendrogram to match
its spectra panel, with colors set by `col`, and boxes only clusters passing
`min_size` and `max_size`.

### Bug fixes and other minor changes

#### Preprocessing and baseline correction

* `preprocess` no longer fails on chromatograms with a few discontinuous wavelengths ---
`smooth.spline` requires at least four unique values, so spectral smoothing is
now skipped below four wavelengths.
* `preprocess` now rejects a single-valued `dim1` or `dim2` when interpolating
onto it, instead of failing later with an opaque `dimnames` error.
* `preprocess` now reads the time unit after unwrapping nested input, so
`outlier_cutoff` is scaled correctly for chromatograms supplied as
`list(dad = <matrix>)`.
* `preprocess` now ignores missing values when rescaling to `maxI`.
* Spectral smoothing in `preprocess` now accounts for the actual spacing of the
wavelength axis rather than treating the columns as evenly spaced. Results are
unchanged for evenly spaced grids, including all of the bundled data.
* Corrected the documented units of `outlier_cutoff` in `preprocess`, which is
expressed in minutes.
* Corrected the documented order of operations in `preprocess`: the code
interpolates first, then smooths spectrally, then baselines.
* `preprocess`, `correct_rt` and `get_peaks` now stop with an error when a
sample fails in parallel (`cl > 1`). Previously they returned `try-error` or
`NULL` elements, including for samples that had succeeded but were processed
by the failing worker.
* Warnings raised while processing samples in parallel (`cl > 1`) now reach the
user. Each distinct warning is raised once, naming the samples that raised it.

#### Retention time correction

* `correct_rt(alg = "vpdtw")` now reproduces the alignment VPdtw fitted.
Previously samples could remain misaligned when the warp extended past the end
of the reference.
* `correct_rt` no longer misplaces peaks when `n_zeros > 0`.
* `correct_peaks` no longer errors on models fitted by `correct_rt` with
`scale = TRUE` or `n_zeros > 0` when `chrom_list` is omitted, and now gives a
clear error for samples with no warping model and for peak lists built with
`collapse = TRUE`.
* `correct_rt(models = ...)` now matches models to samples by name. Previously
models were applied by position, so reordering `chrom_list` warped each sample
with another sample's model.
* `correct_peaks` now replaces the corrected times of an already corrected peak
list. Previously it added a second set, and `get_peaktable` kept using the first.
* `correct_rt` and `correct_rt_group` now return a `chrom_list` when given one,
so `plot()` and chromConverter's `print` method still apply after alignment.

#### Peak detection and integration

* `get_peaks` no longer inflates the areas of peaks bordering missing values,
such as the padded ends of chromatograms aligned by `correct_rt(alg = "vpdtw")`,
which could be several times too large.
* `get_peaks` now converts `sd`, `FWHM`, `tau` and `area` correctly when the
time axis is not in minutes. `time_unit` now defaults to the chromatogram's own
unit, so default output is unchanged.
* `find_peaks(bounds = TRUE)` now returns its full four-column schema when no
peaks are detected, rather than dropping `lower` and `upper`.
* `max_iter` now controls how long the peak fitters run. Previously it was
ignored, leaving `nls.lm`'s default of `100 * (parameters + 1)` evaluations.
`"egh"` and `"bemg"` fits may shift slightly.
* `get_peaks` no longer errors when a sample contains no peaks; previously a
single empty sample failed the whole call. `fit_peaks` now returns an empty
`data.frame` rather than a one-row matrix of `NA`s.
* Printing a `peak_list` no longer errors when a sample has no peaks at any
wavelength, as in peak lists read by `chromConverter::read_peaklist`.
* `FWHM` is now computed as `2*sqrt(2*log(2))*sigma` rather than the rounded
`2.35*sigma`.
* Requesting a wavelength missing from a single-column chromatogram now errors,
instead of returning the existing column falsely labeled with the requested wavelength.
* `fit = "bemg"` no longer returns its starting values in place of a fit when a
peak's integration window spans about 700 timepoints or more.
* `fit = "bemg"` now rejects fits whose tail grows away from the apex, which
inflated the reported area. Such peaks are refitted, and flagged
`converged = FALSE` if the refit also fails.
* `fit = "gaussian"` and `"egh"` peaks can no longer narrow below one timepoint
at half height. Previously a fit to a noise apex could collapse between two
samples and report a height several times the signal around it.

#### Peak tables and filtering

* `get_peaktable` no longer spreads `pk_meta` across several columns when all
peaks at a wavelength fall in one cluster.
* `get_peaktable` no longer rounds `pk_meta` to two decimal places, which moved
retention times by up to 0.005 minutes and set small `tau` values to zero.
* Peak tables built from `correct_peaks` output now read spectra and traces
from aligned chromatograms, rebuilt from the stored warping models. Previously
they were read from the unaligned data at aligned times. See `?correct_peaks`.
* `filter_peaktable` no longer appends a spurious `rt` row to `pk_meta` when
filtering a table built from corrected peaks.
* `normalize_data` no longer reports `column` as ambiguous when another
column's name merely contains it, such as `"V12"` for `"V1"`.
* `combine_peaks` now removes dropped peaks from the `ref_spectra` element so
that it stays aligned with `tab`.
* `combine_peaks(choose = "least_sparse")` now keeps the peak detected in the
most samples. Previously every candidate tied and the first was kept.
* `filter_peaktable` now errors when `rts` is combined with `min_rt` or `max_rt`.
* `write_peaktable(what = "args")` now writes a CSV that reads back. Previously
arguments holding several values, such as two or more wavelengths, spilled into
extra columns.
* `reshape_peaktable` now keeps `lambda` and `rt` when peaks are renamed with a
named `peaks` vector. Previously both were `NA`.
* `reshape_peaktable` now errors on `metadata` columns it cannot find, which
were previously dropped, and gives clear errors for unattached metadata or
unknown `peaks`.

#### Plotting

* `plot_all_spectra` axes, and the absorbance axis of `plot_spectrum`'s ggplot2
and plotly traces, now show wavelength and absorbance units instead of the time
unit.
* `plot_all_spectra` now names the returned spectra when `idx` gives sample
names, and errors on names not in `chrom_list`.
* `plot.peak_table(idx = "all")` no longer errors when `chrom_list` is omitted,
and `attach_ref_spectra`, `mirror_plot` and `get_peaktable(clust = "sp")` now
treat `chrom_list = NULL` as omitted. Omitting `lambdas` for multi-wavelength
data in `mirror_plot` now gives a clear error.
* `plot_all_spectra(what = "rt")` and `what = "idx"` no longer require a
`peak_table`, and a `loc` outside every chromatogram now gives a clear error.
* `plot_chroms(engine = "base")` now extends the y axis below zero when the
traces do. Previously it always started at 0, cutting off negative values.
* The base-engine legend in `plot_chroms` no longer covers the traces at a top
or bottom `legend_position`, including the default `"topright"`; unless `ylim`
is set, the y axis extends to make room.
* `plot_chroms` now fits the y axis to the traces within `xlim` when `ylim` is
not given, under the `"base"` and `"plotly"` engines. Previously it spanned the
whole run.
* `plot_spectrum_inset` now draws spectra from the supplied `chrom_list` when
`idx` covers several samples. Previously it read the peak table's own
chromatograms instead.
* `plot.peak_list(numbers = TRUE)` now draws each number just above its peak.
Previously the numbers were drawn at arbitrary heights.
* `plot.peak_list` now draws shaded peaks at the right width for peak lists
built with `get_peaks(time_unit = "s")` or `"ms"`. Previously they were drawn
too wide.
* `plot.peak_list` now accepts `xlab`, `ylab` and `type`. Previously passing any
of them gave an error.
* `plot.peak_list` no longer fails with "subscript out of bounds" when the
wavelength given to `get_peaks` is written differently from the chromatogram's
column name, such as `"220"` for `"220.00000"`.
* `plot_chroms` no longer errors on single-wavelength chromatograms under the
`"ggplot"` and `"plotly"` engines. A wavelength given in `lambdas` for such data
must now match the column name. Previously it was ignored.
* `plot_chroms` and `reshape_chroms` no longer fail with
`invalid '(to - from)/by'` when the sampling interval is more than twice the
requested `time_resolution`.
* `reshape_chroms` no longer errors when `lambdas` is left at its `NULL`
default, which now returns every wavelength, or when `sample_var` is not
`"sample"`.
* `plot_chroms` now draws each wavelength as a separate line under the
`"ggplot"` and `"plotly"` engines, instead of joining them into one.
* `plot_chroms_heatmap(engine = "base")` now labels its color bar with the
range of the plotted colors. Previously the labels came from the wrong rows and
were rounded to the nearest 100, so they could all read 0.
* The `plot.ptw_list` heatmap now plots the wavelength in `lambdas` and puts
both panels on one color scale. Previously it always used the first wavelength,
the base engine showed only the warped panel, and the other engines' legend
matched only the raw panel.
* `scan_chrom(export_spectrum = TRUE)` now returns the spectrum instead of the
scan index.
* `plot.peak_table` now warns when several `loc` are given but only the first
is returned.

#### Documentation and packaging

* Updated the GC-FID vignette to fetch its example data from a mirror on GitHub, since Dryad now requires authentication for all downloads.

# chromatographR 0.8.0

### Breaking changes

* Added support for bidirectional exponentially modified gaussian (BEMG) peak
shape (`fit = "bemg"` in `get_peaks`/`fit_peaks`). The default peak shape
remains `"egh"` for now but this may change in a future release.
* Added `summarize_by` argument to `get_peaktable` controlling how peak metadata
are summarized within each cluster. The new default, `"weighted.mean"`, weights
peaks by intensity; `"mean"`, `"median"` and `"max"` are also available. Set
`summarize_by = "mean"` to reproduce previous behavior.
* Added automatic thinning of chromatograms in `plot_chroms` to speed up
plotting of raw data, adjustable with the new `time_resolution` argument.
* Changed behavior of `plot_chroms` so the legend is no longer displayed by
default.
* When `export_spectrum == TRUE`, `plot_spectrum` and `plot_all_spectra` now
return spectra invisibly.
* Corrected `get_peaks` to return model-fitted peak heights rather than raw peak
heights. This changes numerical output for workflows relying on the previous
behavior.
* Added `fill_zeros` to `correct_rt`, controlling whether out-of-bounds regions
produced by warping are filled with zeros or returned as `NA`. Defaults to
`FALSE`, which changes VPdtw behavior.
* Changed the name of the `a` parameter in `plot.peak_list` to `alpha`.
* Renamed the `chrom_list` argument to `x` in `plot_chroms_heatmap`, which now
also accepts a peak table in the first position.
* Removed `index` and `chr` arguments across plotting functions, previously
deprecated in v0.7.0. Use `idx` instead.
* Corrected `merge_peaks` so the result no longer depends on the order in which
`peaks` are supplied.
* `get_purity` now decides which timepoints belong to a peak from the
wavelength of maximum absorbance at the apex, rather than the first wavelength
in the chromatogram.
* Deprecated arguments in `get_peaks` still using dot notation (e.g., `sd.max`).
Use the snake_case equivalents instead.
* Deprecated `peak` argument in `plot_all_spectra` in favor of `loc`.
* Deprecated the `save` argument in `cluster_spectra` in favor of `outfile`.
* Removed deprecated `peak_no` argument from `cluster_spectra`.

### New features

#### New functions

* Added `subtract_blanks` function for subtracting blank chromatograms from a
`chrom_list`.
* Added `correct_rt_group` function for group-wise warping function estimation
using within-group averages. This improves alignment when samples share
batch-level shifts or contain group-specific peaks.
* Added `plot_spectrum_inset` function for plotting UV spectra overlaid onto a
chromatogram.
* Added `annotate_peaks` function for adding text labels to peaks in plots
created by `plot_chroms`.
* Added `get_spectra` function for extracting spectra from chromatograms. The
`export_spectrum` argument in `plot_spectrum` and `plot_all_spectra` is
deprecated in favor of `get_spectra`.

#### New generic functions

* Added `summary` method for `peak_table` objects, reporting dimensions, retention time range, intensity statistics, and metadata availability.
* Added `print.peak_table` function to return just the peak data.
* Added `[` method for `peak_table` objects for simultaneous subsetting of `tab`, `pk_meta`, `sample_meta`, and `ref_spectra` elements.
* Added `[` method for `peak_list` objects to preserve metadata attributes when subsetting.
* Added `print` method for `peak_list` objects, summarizing samples, wavelengths, fit method, and peak counts.

#### Enhancements

* Added support for fitting a local baseline beneath each peak through the
`baseline` argument in `get_peaks`: `"none"` (default), `"flat"` or `"sloped"`.
Peak heights and areas are computed above the fitted baseline.
* Added support for plotting spectra from multiple chromatograms in `plot_spectrum` by supplying a vector to `idx` (e.g. `idx = c(1, 2, 4)`). Overlaid traces and spectra are color-coded by sample.
* Added new options to `combine_peaks` via the `choose` argument:
`"least_sparse"` retains the peak detected in the most samples, and `"lambda"`
retains the peak measured at a preferred wavelength.
* Preprocess can now take nested lists containing a `dad` element in each list. This change allows better compatibility with chromConverter.
* Added more informative warnings and errors to `preprocess` function.
* Added support for numeric input to lambda argument in `reshape_chroms`.
* Refactored `normalize_data` so it detects metadata and peak table columns
automatically; added an `on_invalid` argument controlling how zero, negative and
`NA` normalization values are handled; and improved performance.
* Refactored `ptw` warping in `correct_rt` to interpolate the whole matrix at once instead of column by column.
* Added `max_zeros` argument in `filter_peaktable` to exclude features above a threshold proportion of zero values.
* Added support for new `chrom_list` class introduced in chromConverter v0.9.0.
* Allow `plot_chroms` and `plot_chroms_heatmap` to accept peak table as first argument.

### Bug fixes and other minor changes

* Corrected two minor bugs in `get_purity`: one which counted a pure timepoint
as impure when the estimated noise variance was zero, and another which returned
`NA` without warning when only one timepoint fell below the noise threshold.
* Fixed `preprocess` function incorrectly truncating chromatograms when outliers are present and `dim1` is not specified.
* Fixed input validation in `preprocess` so lists containing non-matrix elements no longer pass validation without error.
* Fixed `filter_peaktable` and `filter_peaks`: `rts` with a non-zero `tol`
matched every peak within the tolerance rather than the nearest one, so a single
requested time could select most of the table. Also, criteria are now inclusive
(`min_rt = 5` retains `rt == 5`), and peaks with `NA` in some samples are no
longer excluded by `min_value`.
* Added error when duplicated sample names are provided to `get_peaktable` or `get_peaks` since this will cause downstream issues.
* Fixed issue with visual tests due to new version of Kaleido.
* Added additional arguments as placeholder for `subset.peak_table`.
* Fixed bug in `normalize_data` causing loss of sample names when normalizing chromatograms.
* Formalized requirement for R version 4.1.0.
* Fixed `plot_spectrum` plots with plotly engine so they always maintain their axis labels.
* Standardized axis labels across `plot_spectrum` plots.
* Standardized titles of trace plots across engines so they always use unicode lambda character.
* Fixed `reshape_chroms` so it returns a factor for the `sample` column.
* Fixed `plot_chroms` so chromatograms are black when plotting a single chromatogram with plotly or ggplot.
* Fixed `plot_all_spectra` to consistently export spectra in wide format regardless of plotting engine.
* Fixed error in equation describing spectral clustering in documentation for `get_peaktable`.
* Fixed `get_peaktable` (with `clust = "sp.rt"`) to no longer take the absolute value of spectral similarity scores, which caused negatively correlated spectra to be treated as similar.
* Fixed bug in FWHM calculation in `get_peaktable` for raw peak shape.
* Suppressed spurious "essentially perfect fit" warning from `summary.lm` when computing R² in `fit_peaks`.
* Fixed `subset.peak_table()`, where `sample_meta` was subset on the wrong
condition.
* Made substantial revisions to documentation to improve clarity.

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
