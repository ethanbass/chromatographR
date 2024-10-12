# chromatographR <a href='https://ethanbass.github.io/chromatographR/'><img src='man/figures/logo.png' align="right" height="160" /></a>

<!-- badges: start -->
  [![R-CMD-check](https://github.com/ethanbass/chromatographR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ethanbass/chromatographR/actions/workflows/R-CMD-check.yaml)
  [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/chromatographR)](https://cran.r-project.org/package=chromatographR)
  [![chromatographR status badge](https://ethanbass.r-universe.dev/badges/chromatographR)](https://ethanbass.r-universe.dev)

  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7016988.svg)](https://doi.org/10.5281/zenodo.7016988)
  [![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/chromatographR)](https://cran.r-project.org/package=chromatographR)
  [![metacran downloads](https://cranlogs.r-pkg.org/badges/last-month/chromatographR)](https://cran.r-project.org/package=chromatographR)
   <!-- badges: end -->

## Overview
`chromatographR` is a package for the reproducible analysis of HPLC-DAD chromatographic data in R. It can also be used to analyze other "simple" chromatographic data like GC-FID, HPLC-UV, or HPLC-FD.
  
## Installation

chromatographR can now be installed from CRAN:

```
install.packages("chromatographR")
```

However, it's recommended to install the latest development version of chromatographR from GitHub using the devtools package:

```
install.packages("remotes")
remotes::install_github("https://github.com/ethanbass/chromatographR/")
```

or from R Universe:

```
install.packages("chromatographR", repos="https://ethanbass.r-universe.dev/", type="source")
```

## Usage

#### Importing data
chromatographR can import a variety of vendor formats, including 'Agilent ChemStation' and 'MassHunter' files. This is accomplished using parsers from my chromConverter package. See the [chromConverter page](https://ethanbass.github.io/chromConverter/) for a detailed list of supported formats. Alternatively, chromatographR can also be used with regular `csv` files.

#### Analysis
Please see the [vignette](https://ethanbass.github.io/chromatographR/articles/chromatographR.html) included with the package for details on the application of chromatographR for the analysis of HPLC data. A second vignette with a suggested workflow for the analysis of GC-FID data will be forthcoming soon. 

## Contributing

Contributions are always welcome. Please get in touch (preferable by opening a GitHub [issue](https://github.com/ethanbass/chromatographR/issues)) to discuss any suggestions or to file a bug report. Some good reasons to file an issue:

- You've found an actual bug.  
- You're getting a cryptic error message that you don't understand.  
- You have a file format you'd like to read that isn't currently supported by chromatographR.  (If you do this, please make sure to include a link to an example file!)  
- You have a new feature you'd like to see implemented.  

Also see the [contributing.md](https://github.com/ethanbass/chromatographR/blob/master/.github/contributing.md) page for more details.

(**Note:** Please post questions about file conversions to the [chromConverter](https://github.com/ethanbass/chromConverter/issues) page).

## Citation:

If you use chromatographR in published work, please cite it as follows:

Bass, E. (2023). chromatographR: Chromatographic Data Analysis Toolset (version 0.7.2). http://doi.org/10.5281/zenodo.6944334
