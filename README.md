# chromatographR 
<!-- badges: start -->
  [![R-CMD-check](https://github.com/ethanbass/chromatographR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ethanbass/chromatographR/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

`chromatographR` is a package for the reproducible analysis of HPLC-DAD chromatographic data in R.
  
### Installation

You can install `chromatographR` from GitHub using the devtools package:
```
install.packages("devtools")
devtools::install_github("https://github.com/ethanbass/chromatographR/")
```

To build the vignette, include the argument `build_vignettes=TRUE` (**Note:** this will take considerably longer than loading the package without building the vignette).

### Usage
Please see the vignette included with the package for details on the application of the `chromatographR` package for the analysis of HPLC data. `chromatographR` can now import chemstation and masshunter file formats (using the parsers included with the [Aston](https://github.com/bovee/aston) package for Python) as well as regular `csv` files.
