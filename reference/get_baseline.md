# Retrieve estimated baseline

Retrieves the baseline retained by
[`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md)
or
[`correct_baseline`](https://ethanbass.github.io/chromatographR/reference/correct_baseline.md)
when called with `attach_baseline = TRUE`. This function does not
estimate a baseline itself.

## Usage

``` r
get_baseline(x)
```

## Arguments

- x:

  A chromatogram, or a list of chromatograms, returned by
  [`preprocess`](https://ethanbass.github.io/chromatographR/reference/preprocess.md)
  or
  [`correct_baseline`](https://ethanbass.github.io/chromatographR/reference/correct_baseline.md)
  with `attach_baseline = TRUE`.

## Value

A matrix of the same dimensions as `x` containing the estimated
baseline, or a list of such matrices if `x` is a list.

## Details

The baseline is dropped by
[`correct_rt`](https://ethanbass.github.io/chromatographR/reference/correct_rt.md),
which warps the signal along the time axis, and by
[`subtract_blanks`](https://ethanbass.github.io/chromatographR/reference/subtract_blanks.md),
after which it no longer describes the signal. Retrieve it before either
step.
[`crop_chroms`](https://ethanbass.github.io/chromatographR/reference/crop_chroms.md)
keeps it, cropped to match.

## See also

Other accessors:
[`get_lambdas()`](https://ethanbass.github.io/chromatographR/reference/get_lambdas.md),
[`get_time_resolution()`](https://ethanbass.github.io/chromatographR/reference/get_time_resolution.md),
[`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md)

## Examples

``` r
data(Sa)
pr <- preprocess(Sa[[1]], dim1 = seq(10, 18.66, by = .01),
                 dim2 = seq(200, 318, by = 2), attach_baseline = TRUE)
bl <- get_baseline(pr)
```
