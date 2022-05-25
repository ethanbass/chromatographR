# chromatographR <a href='https://ethanbass.github.io/chromatographR/'><img src='man/figures/logo.png' align="right" height="160" /></a>

<!-- badges: start -->
  [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/chromatographR)](https://cran.r-project.org/package=chromatographR)
  [![metacran downloads](https://cranlogs.r-pkg.org/badges/last-week/chromatographR)](https://cran.r-project.org/package=chromatographR)
  [![R-CMD-check](https://github.com/ethanbass/chromatographR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ethanbass/chromatographR/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

## Overview
`chromatographR` is a package for the reproducible analysis of HPLC-DAD chromatographic data in R.
  
## Installation

`chromatographR` can now be installed from CRAN:

```
install.packages("chromatographR")
```

You can also install the latest development version of `chromatographR` from GitHub using the devtools package:

```
install.packages("devtools")
devtools::install_github("https://github.com/ethanbass/chromatographR/")
```

To build the vignette, include the argument `build_vignettes=TRUE` (**Note:** this will take considerably longer than building the package without the vignette).

## Usage
Please see the [vignette](https://ethanbass.github.io/chromatographR/articles/chromatographR.html) included with the package for details on the application of `chromatographR` for the analysis of HPLC data. `chromatographR` can now import ChemStation and MassHunter file formats (using the parsers included with the [Aston](https://github.com/bovee/aston) package for Python) as well as regular `csv` files. A second vignette with a suggested workflow for the analysis of GC-FID data will be forthcoming soon. 

## Contributing

Contributions to the package are very welcome. Please get in touch (preferable by opening a GitHub [issue](https://github.com/ethanbass/chromatographR/issues)) to discuss any suggestions or to file a bug report. If you get a cryptic error message that you can't understand, that would also be a good reason to file an issue. You can also file an issue if you have a file format you'd like to read that isn't currently supported by chromatographR. (If you do this, please make sure to include a link to an example file!)
