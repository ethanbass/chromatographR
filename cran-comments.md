## R CMD check results

0 errors | 0 warnings | 4 notes

*There is a note about the possible mispelling of 3 words in the description: "HPLC (3:27, 9:52)", "chromatograms (10:5)" and "preprocessing (10:61)". None of these words are actually mispelled.

* Two of the notes concern the package VPdtw which is included as a suggested dependency. VPdtw was formerly available on CRAN but has since been removed. The notes are generated because the package is not available (on CRAN) for checking. One of the notes concerns the Rd references to VPdtw: "Package unavailable to check Rd xrefs: 'VPdtw'". I think it appropriate to include these links to the VPdtw documentation, but I am open to alternative suggestions for how to deal with this. If VPdtw is not installed locally, the links do return an "internal server error". For example, I could use weblinks to the documentation on github rather than assuming that it is present locally. I have also been in touch with the package maintainer about returning the package to CRAN, but have not yet received a response.

* Another note concerns three DOI links in the documentation which are supposedly invalid:

    URL: https://doi.org/10.1021/ac034800e
    From: man/correct_rt.Rd
          inst/doc/chromatographR.html
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.1021/ac501530d
    From: man/get_peaktable.Rd
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.1021/ac802041e
    From: man/correct_rt.Rd
    Status: 503
    Message: Service Unavailable
      
I believe these messages can safely be ignored, since all links load fine in my browser. I'm not sure why the error only affects these three DOIs. Hadley Wickam has stated in this (tweet)[https://twitter.com/hadleywickham/status/1358170607314235392] that this type of message can safely be ignored.
