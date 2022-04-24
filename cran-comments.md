## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

* The note concerns two DOI links in the documentation which are supposedly invalid:

      URL: https://doi.org/10.1021/ac034800e
      From: man/correct_rt.Rd
            inst/doc/chromatographR.html
      Status: 503
      Message: Service Unavailable
      
      URL: https://doi.org/10.1021/ac501530d
      From: man/get_peaktable.Rd
      Status: 503
      Message: Service Unavailable
      
I believe these messages can safely be ignored, since both links load fine in my browser. I'm not sure why the error only affects these two DOIs. Hadley Wickam has stated in this (tweet)[https://twitter.com/hadleywickham/status/1358170607314235392] that this type of message can safely be ignored.
