# Get lambdas

Get wavelengths from a chromatogram, a list of chromatograms, or a
`peak_table` object.

## Usage

``` r
get_lambdas(x)
```

## Arguments

- x:

  A chromatogram, a list of chromatograms, or a `peak_table` object.

## Value

A numeric vector of wavelengths, taken from the first chromatogram if
`x` is a list.

## See also

Other accessors:
[`get_baseline()`](https://ethanbass.github.io/chromatographR/reference/get_baseline.md),
[`get_time_resolution()`](https://ethanbass.github.io/chromatographR/reference/get_time_resolution.md),
[`get_times()`](https://ethanbass.github.io/chromatographR/reference/get_times.md)

## Examples

``` r
data(Sa_warp)
get_lambdas(Sa_warp)
#>  [1] 200 202 204 206 208 210 212 214 216 218 220 222 224 226 228 230 232 234 236
#> [20] 238 240 242 244 246 248 250 252 254 256 258 260 262 264 266 268 270 272 274
#> [39] 276 278 280 282 284 286 288 290 292 294 296 298 300 302 304 306 308 310 312
#> [58] 314 316 318
```
