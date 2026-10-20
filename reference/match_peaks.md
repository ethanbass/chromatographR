# Match peaks across chromatograms

Matches peaks across chromatograms by spectral correlation and retention
time proximity. Each reference peak is matched to at most one sample
peak, taking candidate pairs in order of score, highest first. The
resulting `matched_peak_list` object can be passed to
[`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md)
to perform landmark parametric time warping (`alg = "lptw"`).

## Usage

``` r
match_peaks(
  peak_list,
  chrom_list = NULL,
  reference = NULL,
  targets = NULL,
  lambda = NULL,
  rt_window = 0.3,
  spectral_cutoff = 0.95,
  background_correction = c("none", "linear"),
  rt_weight = 0.15,
  use_corrected_rt = NULL,
  spectral_test = c("cutoff", "noise"),
  alpha = 0.01,
  plot_it = FALSE,
  verbose = getOption("verbose")
)
```

## Arguments

- peak_list:

  A peak list object from
  [`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md).

- chrom_list:

  List of chromatograms in matrix format. If `NULL`, will be retrieved
  from the `peak_list` object.

- reference:

  Index or name of the reference sample, or `"best"` to select the
  reference algorithmically using
  [`bestref`](https://rdrr.io/pkg/ptw/man/bestref.html). Defaults to
  `NULL`, which uses the first sample, or, if `targets` is given, a
  sample chosen from the targets.

- targets:

  Optional table of targets, or a numeric vector of their retention
  times (see
  [`targeted_analysis`](https://ethanbass.github.io/chromatographR/reference/targeted_analysis.md)).
  Used to choose the reference sample, unless `reference` is given, and
  to name the targets' columns in the peak table. See *Targets and the
  reference sample*. Cannot be combined with `reference = "best"`.
  Defaults to `NULL`.

- lambda:

  A character or numeric value specifying the wavelength to use for peak
  matching. Required if peaks were detected at multiple wavelengths,
  unless `targets` has a `lambda` column, which sets the wavelength per
  target; `lambda` then fills in targets without one.

- rt_window:

  Numeric. Maximum allowed difference in retention time, in minutes,
  between a reference and sample peak to be considered a candidate
  match. An `rt_window` column in `targets` overrides it for those
  targets. Defaults to `0.3`.

- spectral_cutoff:

  Numeric. Minimum spectral correlation for a candidate peak pair to be
  retained. A `spectral_cutoff` column in `targets` overrides it for
  those targets. Ignored, with a warning, for one-dimensional
  chromatograms. Defaults to `0.95`.

- background_correction:

  How to correct the apex spectrum of each peak before correlating:
  `"none"` (the default) uses it as is, and `"linear"` subtracts the
  background under the peak, interpolated between the spectra at its
  `start` and `end`. See *Background correction*.

- rt_weight:

  Numeric in \[0, 1\]. Penalizes retention time offset when ranking
  candidate pairs: each pair's spectral correlation is multiplied by
  `1 - rt_weight * d / rt_window`, where `d` is the pair's retention
  time difference. At `rt_weight = 1`, a pair at the edge of the window
  scores zero. Defaults to `0.15`.

- use_corrected_rt:

  Logical. Whether to compare corrected retention times (`rt.cor`) from
  [`correct_peaks`](https://ethanbass.github.io/chromatographR/reference/correct_peaks.md)
  rather than raw ones (`rt`), both between peaks and against `targets`.
  Defaults to `NULL`, which uses corrected times when every sample has
  them.

- spectral_test:

  How a candidate pair's spectra are compared: `"cutoff"` (the default)
  requires their correlation to be at least `spectral_cutoff`; `"noise"`
  tests whether they could be the same spectrum given each sample's
  noise, so that a small, noisy peak is not rejected for its noise
  alone. The `"noise"` test is experimental; see *Testing spectra
  against noise*.

- alpha:

  Significance level for `spectral_test = "noise"`. A pair is rejected
  if its p-value is below `alpha`. Pairs whose spectra are the same
  apart from noise are rejected at most about `alpha` of the time,
  except when the reference peak is small (up to roughly 20 times the
  noise), where they are rejected more often. Pairs whose spectra truly
  correlate at `spectral_cutoff` are rejected more often too, up to
  about half the time for peaks far above the noise, as with a
  correlation threshold. Defaults to `0.01`.

- plot_it:

  Logical. Whether to plot the matches with
  [`plot.matched_peak_list`](https://ethanbass.github.io/chromatographR/reference/plot.matched_peak_list.md).
  Defaults to `FALSE`.

- verbose:

  Whether to print messages to console.

## Value

A `matched_peak_list` object: a list with one element per sample
containing the matched peak pairs, per-pair spectral correlation scores,
and peak data in index space for use with
[`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md).
A sample with no peaks or no matches is `NULL`, with a warning. When
`targets` is given, the targets with `rt` replaced by the reference's
own retention times are stored as the `"reference_targets"` attribute.
When the targets span several wavelengths, peaks are matched separately
at each, and the result is a `matched_peak_list_set`: a list of
`matched_peak_list` objects named by wavelength, which
[`get_peaktable`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)
combines into one table. Unless `reference` is given, every wavelength
uses the same reference sample, chosen from all the targets; if no
sample contains every target, each wavelength gets its own, with a
warning.

## Targets and the reference sample

With `targets`, a sample's peak for a target is its tallest peak within
`tol` of it that, when fitted, has `r.squared` of at least 0.5; a peak
within `tol` of several targets counts only for the nearest. Unless
`reference` is given, the reference is chosen among samples with a peak
for every target, where a peak counts if it is at least a tenth of that
target's median height across samples. The sample whose weakest target
is largest relative to its median is chosen, then the one whose farthest
peak lies closest to its target time.

## Background correction

With `background_correction = "linear"`, the apex spectrum is averaged
over 2 scans either side and each boundary spectrum over 5 scans either
side, to limit the noise the subtraction adds. The subtraction removes a
background that varies linearly across the peak, such as a sloping
baseline or a neighbor's tail. It also adds noise, and where a boundary
lies on a neighboring peak it subtracts part of that peak's spectrum.

## Testing spectra against noise

With `spectral_test = "noise"`, each candidate pair is tested rather
than thresholded. The sample spectrum is fitted as an offset plus a
multiple of the reference spectrum by least squares, each wavelength
weighted by the reciprocal of the sample's noise there. The hypothesis
is that, without noise, the two spectra would correlate at
`spectral_cutoff`; the residual that this mismatch would leave is
subtracted, and what remains is compared with the residuals left under
the same fit by pairs of noise spectra, one from the sample and one from
the reference scaled by the fitted multiple. For a peak far above the
noise, this reduces to requiring a noise-weighted correlation of at
least `spectral_cutoff`; for a small peak, it allows for the noise, so a
real peak is not rejected because its spectrum is noisy. Each matched
sample's `spectral_p` holds the p-values of its matched pairs.

This test is experimental. It has had little testing on real data, and
the way its p-values are calculated may change in future versions.

The noise spectra come from each chromatogram's own peak-free
timepoints: those outside every peak in `peak_list` taller than 5 times
its trace's noise, at any wavelength in `peak_list`, after dropping
timepoints whose spectra still stand out from the rest. A trace's noise
is the scaled median absolute deviation of its timepoints outside every
peak taller than 3 times the noise, found by iterating from the estimate
given by its point-to-point differences, `mad(diff(y)) / sqrt(2)` (von
Neumann 1941; Donoho and Johnstone 1994).

The test cannot tell compounds apart near the noise level, where every
spectrum is consistent with every other; decide whether a small peak is
present by its height instead. A peak overlapped by a neighbor or
sitting on structured background may be rejected because its apex
spectrum really is a mixture. The test needs
`background_correction = "none"` and, in each chromatogram, at least 30
peak-free timepoints and at least `1 / alpha`, since with fewer it could
never reject a match. On short traces with strongly autocorrelated
noise, the peak-free timepoints are fewer independent spectra than their
number, so a `spectral_cutoff` of `1` rejects identical spectra up to
about twice as often as `alpha`.

## References

- Donoho, D. L., and Johnstone, I. M. 1994. Ideal spatial adaptation by
  wavelet shrinkage. *Biometrika* 81:425-455.
  [doi:10.1093/biomet/81.3.425](https://doi.org/10.1093/biomet/81.3.425)
  .

- von Neumann, J. 1941. Distribution of the ratio of the mean square
  successive difference to the variance. *Ann. Math. Stat.* 12:367-395.
  [doi:10.1214/aoms/1177731677](https://doi.org/10.1214/aoms/1177731677)
  .

## See also

[`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md),
[`get_peaks`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md),
[`targeted_analysis`](https://ethanbass.github.io/chromatographR/reference/targeted_analysis.md)

## Author

Ethan Bass

## Examples

``` r
data(Sa_pr)
pks <- get_peaks(Sa_pr, lambdas = 210, fit = "egh")
match_peaks(pks, lambda = 210)
#> A matched_peak_list with 4 samples, matched at 210 nm against reference sample 119
#> Reference peaks: 29 | spectral_cutoff: 0.95 | background correction: none
#> Matched peaks per sample: 121 12, 122 18, 458 17
```
