# chromatographR

`chromatographR` is a package for the reproducible analysis of HPLC-DAD chromatographic data in R.


<!-- badges: start -->
  [![R-CMD-check](https://github.com/ethanbass/chromatographR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ethanbass/chromatographR/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->
  
### Installation

You can install `chromatographR` from GitHub using the devtools package:
```
install.packages("devtools")
devtools::install_github("https://github.com/ethanbass/chromatographR/")
```

To build the vignette, include the argument `build_vignettes=TRUE` (**Note:** this will take considerably longer than loading the package without building the vignette).

### Usage
Please see the vignette included with the package for details on the application of the `chromatographR` package for the analysis of HPLC data. Your data will need to be in `csv` format first. Unfortunately, this can be considerably harder than it sounds, since most HPLC software is not designed to export raw data which is often stored in proprietary, binary formats.
