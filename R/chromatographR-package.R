#' chromatographR: Chromatographic Data Analysis Toolset
#' 
#' Tools for high-throughput analysis of HPLC-DAD/UV chromatographic data. 
#' The package provides functionality for preprocessing, alignment, peak
#' detection and fitting, peak table construction, and visualization of
#' chromatographic data.
#' 
#' A typical workflow includes signal preprocessing and peak table construction
#' following general approaches described in Wehrens et al. (2015). 
#' 
#' Chromatogram alignment is implemented using parametric time warping (PTW) or
#' variable penalty dynamic time warping (VPdtw).
#' 
#' Peak detection is based on zero-crossings in smoothed derivatives. Peaks are
#' then modeled using Gaussian, exponential-Gaussian hybrid, or bidirectional 
#' exponentially modified Gaussian peak functions via nonlinear least squares.
#' 
#' Further details and example workflows are provided in the package vignettes.
#' 
#' @section Core analysis functions:
#' - `read_chroms()`: Import chromatograms from a variety of vendor formats.
#' - `preprocess()`: Preprocess chromatographic matrices.
#' - `correct_rt()`: Align chromatograms.
#' - `get_peaks()`: Find and integrate peaks.
#' - `get_peaktable()`: Assemble peak table.
#' - `attach_metadata()`: Attach metadata to peak table.
#' - `attach_ref_spectra()`: Attach reference spectra to peak table.
#' - `normalize_data()`: Normalize `peak_table` or `chrom_list`.
#' @section Peak refinement:
#' - `combine_peaks()`: Combine duplicate peaks based on retention time and spectral similarity.
#' - `merge_peaks()`: Merge split peaks into a single feature in a peak table.
#' - `filter_peaktable()`: Filter peak features in a peak table.
#' - `filter_peaks()`: Filter peaks in peak lists.
#' @section Visualization functions:
#' ### Chromatogram plots
#' - `plot_chroms()`: Plot chromatograms as traces.
#' - `plot_chroms_heatmap()`: Plot chromatograms as heatmap.
#' - `annotate_peaks()`: Add peak labels to chromatogram plots.
#' - `mirror_plot()`: Plot chromatograms as mirror plots.
#' - `plot.ptw_list()`: Plot parametric time warping (PTW) alignment object.
#' ### Spectral plots
#' - `plot_spectrum()`: Plot spectrum and/or trace of specified peak.
#' - `plot_spectrum_inset()`: Plot spectrum over chromatogram.
#' - `plot_all_spectra()`: Plot all spectra for specified peak.
#' - `scan_chrom()`: Interactively extract spectra from a chromatogram.
#' ### Peak plots
#' - `plot.peak_list()`: Plot fitted peaks over chromatographic trace.
#' - `plot.peak_table()`: Plot chromatograms and/or spectra from a peak table.
#' - `boxplot.peak_table()`: Create boxplot from a peak table object.
#' @section Data utilities and I/O:
#' - `get_times()`: Return retention times from a peak table or a list of chromatograms.
#' - `get_lambdas()`: Return wavelengths from a peak table or a list of chromatograms.
#' - `write_peaktable()`: Export peak table to CSV or XLSX format.
#' - `reshape_chroms()`: Reshape a list of chromatograms to long format.
#' - `reshape_peaktable()`: Reshape a peak table to long format.
#' @section Example data:
#' - `Sa`: A list of four goldenrod root chromatograms.
#' - `Sa_pr`: Preprocessed goldenrod root chromatograms.
#' - `Sa_warp`: Preprocessed and aligned goldenrod root chromatograms.
#' - `pk_tab`: Peak table from aligned goldenrod root chromatograms.
#' @references 
#' - Clifford, D., Stone, G., Montoliu, I., Rezzi, S., Martin, F. P., Guy, P.,
#' Bruce, S., & Kochhar, S. 2009. Alignment using variable penalty dynamic time
#' warping. *Analytical chemistry*, 81(3):1000-1007. \doi{10.1021/ac802041e}.
#'
#' - Clifford, D., & Stone, G. 2012. Variable Penalty Dynamic Time Warping Code
#' for Aligning Mass Spectrometry Chromatograms in R. *Journal of
#' Statistical Software*, 47(8):1-17. \doi{10.18637/jss.v047.i08}.
#' 
#' - Eilers, P.H.C. 2004. Parametric Time Warping.
#' *Analytical Chemistry*, 76:404-411. \doi{10.1021/ac034800e}.
#' 
#' - Wehrens, R., Bloemberg, T.G., and Eilers P.H.C. 2015. Fast
#' parametric time warping of peak lists. *Bioinformatics*, 31:3063-3065.
#' \doi{10.1093/bioinformatics/btv299}.
#' 
#' - Wehrens, R., Carvalho, E., Fraser, P.D. 2015. Metabolite profiling in
#' LC–DAD using multivariate curve resolution: the alsace package for R.
#' *Metabolomics*, 11:143-154. \doi{10.1007/s11306-014-0683-5}.
#' @author Ethan Bass
#' @md
"_PACKAGE"

#' Raw goldenrod root chromatograms
#' 
#' A list of four HPLC-DAD data matrices of *Solidago altissima* roots
#' extracted in 90\% methanol. Retention times are stored in rows and 
#' wavelengths are stored in columns. Data were collected on a Agilent 1100 HPLC.
#' 
#' @name Sa
#' @docType data
#' @keywords data
#' @usage data(Sa)
#' @format  A list of four matrices (1301 times x 60 wavelengths).
#' @family data objects
NULL

#' Preprocessed goldenrod root chromatograms
#' 
#' A list of four pre-processed HPLC-DAD chromatograms derived from the raw data
#' stored in [`Sa`]. Retention times are stored in rows and wavelengths
#' are stored in columns. The time axis is compressed to save space and 
#' processing time so the data are a little choppy. 
#'
#' @name Sa_pr
#' @keywords data
#' @docType data
#' @usage data(Sa_pr)
#' @format  A list of four pre-processed matrices (434 retention times x 60 wavelengths).
#' @family data objects
NULL

#' Warped goldenrod root chromatograms.
#' 
#' A list of four pre-processed and warped goldenrod root chromatograms derived
#' from the raw data stored in [`Sa`].
#'
#' @name Sa_warp
#' @docType data
#' @keywords data
#' @usage data(Sa_warp)
#' @format  A list of four pre-processed and warped matrices (434 times x 
#' 60 wavelengths).
#' @family data objects
NULL

#' Goldenrod peak table
#' 
#' A peak table generated from the exemplary goldenrod root extracts stored in
#' [`Sa`].
#'
#' @name pk_tab
#' @docType data
#' @keywords data
#' @usage data(pk_tab)
#' @format  A `peak_table` object.
#' @family data objects
NULL
