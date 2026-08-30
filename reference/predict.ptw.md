# Predict PTW

Predict method for objects of class `ptw`. Reproduced from the
[`ptw`](https://rdrr.io/pkg/ptw/man/ptw.html) package as the original is
not exported.

## Usage

``` r
# S3 method for class 'ptw'
predict(object, newdata, what = c("response", "time"), RTref = NULL, ...)
```

## Arguments

- object:

  an object of class `ptw`.

- newdata:

  an optional matrix of new data to predict. If missing, returns the
  warped sample.

- what:

  character string specifying whether to return the warped `"response"`
  or corrected `"time"`. Defaults to `"response"`.

- RTref:

  optional numeric vector of reference retention times.

- ...:

  additional arguments (currently ignored).

## Value

a matrix of warped responses or corrected retention times.

## Author

Ron Wehrens
