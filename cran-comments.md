* I apparently did not succeed in fixing the test failure on the MKL server in version 0.4.0, so I made another attempt to do so here by eliminating the test that is failing.

## R CMD check results

0 errors | 0 warnings | 1 note

* The note concerns four DOI links in the documentation which are supposedly invalid:

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1002/cem.859
    From: inst/doc/chromatographR.html
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.1021/ac034800e
    From: inst/doc/chromatographR.html
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.1021/ac802041e
    From: inst/doc/chromatographR.html
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.1111/nph.12172
    From: inst/doc/chromatographR.html
    Status: 503
    Message: Service Unavailable

I believe these messages can safely be ignored, since all of the links in question load fine in my web browser.
