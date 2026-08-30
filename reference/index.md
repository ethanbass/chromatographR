# Package index

## Package overview

- [`chromatographR`](https://ethanbass.github.io/chromatographR/reference/chromatographR-package.md)
  [`chromatographR-package`](https://ethanbass.github.io/chromatographR/reference/chromatographR-package.md)
  : chromatographR: Chromatographic Data Analysis Toolset

## Analysis functions

- [`reexports`](https://ethanbass.github.io/chromatographR/reference/reexports.md)
  [`read_chroms`](https://ethanbass.github.io/chromatographR/reference/reexports.md)
  : Read chromatograms.

- [`preprocess()`](https://ethanbass.github.io/chromatographR/reference/preprocess.md)
  : Preprocess time/wavelength data

- [`subtract_blanks()`](https://ethanbass.github.io/chromatographR/reference/subtract_blanks.md)
  : Subtract blank chromatograms

- [`correct_rt()`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md)
  : Retention time alignment via time warping

- [`correct_rt_group()`](https://ethanbass.github.io/chromatographR/reference/correct_rt_group.md)
  : Correct retention time with group-based warping

- [`correct_peaks()`](https://ethanbass.github.io/chromatographR/reference/correct_peaks.md)
  : Correct peak positions according to a PTW warping model

- [`get_peaks()`](https://ethanbass.github.io/chromatographR/reference/get_peaks.md)
  : Get peak list.

- [`fit_peaks()`](https://ethanbass.github.io/chromatographR/reference/fit_peaks.md)
  : Fit chromatographic peaks

- [`get_peaktable()`](https://ethanbass.github.io/chromatographR/reference/get_peaktable.md)
  : Convert peak lists into an aligned peak table.

- [`get_spectra()`](https://ethanbass.github.io/chromatographR/reference/get_spectra.md)
  : Extract a spectrum from a chromatogram

- [`attach_metadata()`](https://ethanbass.github.io/chromatographR/reference/attach_metadata.md)
  :

  Attaches sample metadata to a `peak_table` object by matching sample
  names.

- [`attach_ref_spectra()`](https://ethanbass.github.io/chromatographR/reference/attach_ref_spectra.md)
  :

  Attach reference spectra to a `peak_table` object.

- [`normalize_data()`](https://ethanbass.github.io/chromatographR/reference/normalize_data.md)
  : Normalize peak table or chromatograms

- [`cluster_spectra()`](https://ethanbass.github.io/chromatographR/reference/cluster_spectra.md)
  : Cluster spectra

## Visualization functions

- [`plot_all_spectra()`](https://ethanbass.github.io/chromatographR/reference/plot_all_spectra.md)
  : Plot all spectra for chosen peak.
- [`plot_chroms()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms.md)
  : Plot traces from list of chromatograms.
- [`plot_chroms_heatmap()`](https://ethanbass.github.io/chromatographR/reference/plot_chroms_heatmap.md)
  : Plot chromatograms as heatmap
- [`plot_spectrum()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum.md)
  : Plot spectrum from peak table
- [`plot_spectrum_inset()`](https://ethanbass.github.io/chromatographR/reference/plot_spectrum_inset.md)
  : Plot a chromatogram or chromatograms with an inset UV spectrum
- [`annotate_peaks()`](https://ethanbass.github.io/chromatographR/reference/annotate_peaks.md)
  : Annotate peaks on a chromatogram
- [`boxplot(`*`<peak_table>`*`)`](https://ethanbass.github.io/chromatographR/reference/boxplot.peak_table.md)
  : Make boxplot from peak table.
- [`mirror_plot()`](https://ethanbass.github.io/chromatographR/reference/mirror_plot.md)
  : Make mirror plot from peak table.
- [`scan_chrom()`](https://ethanbass.github.io/chromatographR/reference/scan_chrom.md)
  : Plot spectra by clicking on the chromatogram.
- [`plot(`*`<peak_list>`*`)`](https://ethanbass.github.io/chromatographR/reference/plot.peak_list.md)
  : Overlay fitted peak shapes on chromatograms
- [`plot(`*`<peak_table>`*`)`](https://ethanbass.github.io/chromatographR/reference/plot.peak_table.md)
  : Plot spectrum from peak table
- [`plot(`*`<ptw_list>`*`)`](https://ethanbass.github.io/chromatographR/reference/plot.ptw_list.md)
  : Plot PTW alignments

## Utility functions

- [`combine_peaks()`](https://ethanbass.github.io/chromatographR/reference/combine_peaks.md)
  : Combine duplicate peaks
- [`filter_peaks()`](https://ethanbass.github.io/chromatographR/reference/filter_peaks.md)
  : Filter peak lists
- [`filter_peaktable()`](https://ethanbass.github.io/chromatographR/reference/filter_peaktable.md)
  : Filter peak table
- [`merge_peaks()`](https://ethanbass.github.io/chromatographR/reference/merge_peaks.md)
  : Merge split peaks
- [`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md)
  : Get retention times
- [`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md)
  : Get lambdas
- [`reshape_chroms()`](https://ethanbass.github.io/chromatographR/reference/reshape_chroms.md)
  : Reshape chromatograms
- [`reshape_peaktable()`](https://ethanbass.github.io/chromatographR/reference/reshape_peaktable.md)
  : Reshape peaktable
- [`write_peaktable()`](https://ethanbass.github.io/chromatographR/reference/write_peaktable.md)
  : Export peak table

## Data structures

- [`peak_table-class`](https://ethanbass.github.io/chromatographR/reference/peak_table-class.md)
  [`peak_table`](https://ethanbass.github.io/chromatographR/reference/peak_table-class.md)
  : Peak Table Object

## Methods

- [`print(`*`<peak_table>`*`)`](https://ethanbass.github.io/chromatographR/reference/print.peak_table.md)
  : Print peak table

- [`print(`*`<peak_list>`*`)`](https://ethanbass.github.io/chromatographR/reference/print.peak_list.md)
  : Print peak list

- [`summary(`*`<peak_table>`*`)`](https://ethanbass.github.io/chromatographR/reference/summary.peak_table.md)
  : Summarize Peak Table

- [`` `[`( ``*`<peak_table>`*`)`](https://ethanbass.github.io/chromatographR/reference/sub-.peak_table.md)
  : Subset peak table

- [`` `[`( ``*`<peak_list>`*`)`](https://ethanbass.github.io/chromatographR/reference/sub-.peak_list.md)
  : Subset peak list

- [`dim(`*`<peak_table>`*`)`](https://ethanbass.github.io/chromatographR/reference/dim.peak_table.md)
  :

  Return dimensions of a `peak_table` object.

- [`head(`*`<peak_table>`*`)`](https://ethanbass.github.io/chromatographR/reference/head.peak_table.md)
  [`tail(`*`<peak_table>`*`)`](https://ethanbass.github.io/chromatographR/reference/head.peak_table.md)
  :

  Return first or last parts of a `peak_table`.

- [`row.names(`*`<peak_table>`*`)`](https://ethanbass.github.io/chromatographR/reference/row.names.peak_table.md)
  :

  Return row names from a `peak_table` object.

- [`subset(`*`<peak_table>`*`)`](https://ethanbass.github.io/chromatographR/reference/subset.peak_table.md)
  : Subset peak table

## Data objects

- [`Sa`](https://ethanbass.github.io/chromatographR/reference/Sa.md) :
  Raw goldenrod root chromatograms
- [`Sa_pr`](https://ethanbass.github.io/chromatographR/reference/Sa_pr.md)
  : Preprocessed goldenrod root chromatograms
- [`Sa_warp`](https://ethanbass.github.io/chromatographR/reference/Sa_warp.md)
  : Warped goldenrod root chromatograms.
- [`pk_tab`](https://ethanbass.github.io/chromatographR/reference/pk_tab.md)
  : Goldenrod peak table
