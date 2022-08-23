## R CMD check results

There was an issue with some of the tests in my last submission that was only showing up when tested on the MKL machine (https://www.stats.ox.ac.uk/pub/bdr/Rblas/MKL/chromatographR.out). I believe I have remedied these issues in this new release. I apologize for the double release!

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

*There's an additional note when I used rhub to build the package on the Windows Server 2022, R-devel, 64 bit.

checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
  
I don't really understand why this is being generated, but I think it can likewise be ignored.
