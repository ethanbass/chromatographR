# chromatographR 0.4.0

* Added support for variable penalty dynamic time warping (VPdtw) through `correct_rt` function.
* Fixed bug in get_peaks function.
* Allow preprocessing without interpolation.
* Fixed bug so preprocess can work on windows (without parallel processing).
* Allow use of raw data for peak integration in `get_peaks`.
* Added `verbose` option to `correct_rt` to print reference chromatogram.

# chromatographR 0.3

* Added a `NEWS.md` file to track changes to the package.
